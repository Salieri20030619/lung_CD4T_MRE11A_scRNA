#!/usr/bin/env bash
set -Eeuo pipefail
ROOT="${ROOT:-/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA}"
mkdir -p "$ROOT/logs" "$ROOT/tmp" "$ROOT/02_analysis"
OUT=$(mktemp -d "$ROOT/02_analysis/CD4_DEG_collect_$(date +%Y%m%d_%H%M%S)_XXXXXX")
export ROOT OUT
trap 'rc=$?; printf "CD4_DEG_COLLECT_FAIL EXIT_CODE=%s OUTPUT=%s\n" "$rc" "$OUT"; exit "$rc"' ERR
printf 'OUTPUT=%s\n' "$OUT"
cat > "$OUT/collect.py" <<'PY'
import csv
import gzip
import json
import os
import re
import shutil
import tarfile
from pathlib import Path

root = Path(os.environ['ROOT']).resolve()
out = Path(os.environ['OUT']).resolve()
limit = 25 * 1024 * 1024
effects = {'logfc', 'log2fc', 'log2foldchange', 'avglog2fc', 'avglogfc'}
pvals = {'pvalue', 'pval', 'pvaladj', 'padj', 'fdr', 'p', 'pvalues', 'pvals', 'pvalsadj'}
def norm(s):
    return re.sub(r'[^a-z0-9]', '', s.lower())

def header(p):
    opener = gzip.open if p.name.lower().endswith('.gz') else open
    with opener(p, 'rt', encoding='utf-8-sig', errors='strict') as f:
        line = f.readline(131072)
    delimiter = '\t' if '\t' in line else ','
    return next(csv.reader([line], delimiter=delimiter))

def write_table(name, rows, fields):
    with (out / name).open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(rows)

paths = []
issues = []
for dirname in ['02_analysis', '07_DEG', '10_tables', 'scripts']:
    base = root / dirname
    if not base.exists():
        continue
    def walk_error(e):
        issues.append({'path': str(e.filename), 'reason': str(e)})
    for folder, dirs, files in os.walk(base, onerror=walk_error, followlinks=False):
        dirs[:] = sorted(d for d in dirs if not d.startswith('CD4_DEG_collect_') and d not in {'__pycache__', '.git'})
        for name in sorted(files):
            p = Path(folder) / name
            if p.is_symlink():
                continue
            paths.append(p)

inventory = []
selected = {}
deg_parents = set()
context_pattern = re.compile(r'design|sample|metadata|provenance|method|manifest|summary|status|readme|contrast', re.I)
for p in paths:
    rel = p.relative_to(root).as_posix()
    low = p.name.lower()
    is_table = low.endswith(('.tsv', '.csv', '.txt', '.tsv.gz', '.csv.gz'))
    if not is_table:
        continue
    try:
        size = p.stat().st_size
        if size > limit:
            if re.search(r'deg|differ|pseudobulk|marker', rel, re.I):
                issues.append({'path': rel, 'reason': 'TABLE_OVER_25_MIB_NOT_READ'})
            continue
        cols = header(p)
        keys = {norm(c) for c in cols}
        eff = keys & effects
        pv = keys & pvals
        if eff and pv:
            selected[p] = 'DE_TABLE'
            deg_parents.add(p.parent)
            inventory.append({'path': rel, 'bytes': size,
                'mtime': p.stat().st_mtime,
                'effect_columns': ';'.join(c for c in cols if norm(c) in effects),
                'p_columns': ';'.join(c for c in cols if norm(c) in pvals),
                'columns': ';'.join(cols),
                'scope_hint': 'CD4_PATH' if 'cd4' in rel.lower() else 'REVIEW_REQUIRED',
                'comparison_status': 'REVIEW_SCRIPT_AND_DESIGN'})
    except Exception as e:
        issues.append({'path': rel, 'reason': f'{type(e).__name__}: {e}'})

# Include nearby code and sample design so marker tests, time contrasts and
# disease contrasts can be distinguished before ranking any gene.
near = set(deg_parents)
for p in deg_parents:
    if p.parent != root:
        near.add(p.parent)
for p in paths:
    rel = p.relative_to(root).as_posix()
    name = p.name.lower()
    code = name.endswith(('.r', '.py', '.sh'))
    context = context_pattern.search(p.name) and name.endswith(('.tsv', '.csv', '.txt', '.json', '.md'))
    nearby = p.parent in near or p.parent.parent in near
    project_code = p.parent == root / 'scripts' and re.search(r'cd4|gse314943|gse274372|gse183682|gse174725|gse288887', name)
    try:
        size = p.stat().st_size
        if size <= 4 * 1024 * 1024 and ((nearby and (code or context)) or (project_code and code)):
            selected.setdefault(p, 'CONTEXT')
    except OSError as e:
        issues.append({'path': rel, 'reason': str(e)})

manifest = []
for p, role in sorted(selected.items(), key=lambda x: str(x[0])):
    rel = p.relative_to(root)
    dest = out / 'source' / rel
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(p, dest)
    manifest.append({'path': rel.as_posix(), 'role': role, 'bytes': dest.stat().st_size})
write_table('DE_table_inventory.tsv', inventory, ['path', 'bytes', 'mtime', 'effect_columns', 'p_columns', 'columns', 'scope_hint', 'comparison_status'])
write_table('collected_files.tsv', manifest, ['path', 'role', 'bytes'])
write_table('read_issues.tsv', issues, ['path', 'reason'])
status = {'DE_tables': len(inventory), 'collected_files': len(manifest), 'read_issues': len(issues),
          'purpose': 'Existing result collection; no differential test or gene ranking performed',
          'scope': 'Each file retains its original project-relative path; versions are not merged',
          'next': 'Confirm CD4 identity, biological replication, contrasts and source correction before ranking'}
(out / 'status.json').write_text(json.dumps(status, ensure_ascii=False, indent=2) + '\n')
print(json.dumps(status, ensure_ascii=False), flush=True)
print('ANALYSIS_LOGIC_PASS=COLLECTION_ONLY_SOURCE_PATHS_PRESERVED', flush=True)
archive = out.with_suffix('.tar.gz')
with tarfile.open(archive, 'w:gz') as t:
    for p in sorted(out.rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts:
            t.add(p, arcname=p.relative_to(out).as_posix(), recursive=False)
with tarfile.open(archive, 'r:gz') as t:
    names = set(t.getnames())
    assert {'status.json', 'DE_table_inventory.tsv', 'collected_files.tsv'} <= names
stable = root / 'tmp' / 'lung_CD4_DEG_collect_return.tar.gz'
temp = stable.with_name(stable.name + '.' + str(os.getpid()) + '.part')
shutil.copyfile(archive, temp)
os.replace(temp, stable)
print('RETURN_PACKAGE=' + str(stable), flush=True)
print('DE_TABLES_FOUND=' + str(len(inventory)), flush=True)
print('READ_ISSUES=' + str(len(issues)), flush=True)
print('CD4_DEG_COLLECT_PASS' if inventory else 'CD4_DEG_COLLECT_NO_TABLES', flush=True)
PY
bash -n "$0"
conda run -n batch_bench python -m py_compile "$OUT/collect.py"
echo 'SYNTAX_CHECK=PASS'
conda run --no-capture-output -n batch_bench python -u "$OUT/collect.py"
