import csv
import math
import os
import re
import shutil
import tarfile
from pathlib import Path

root, snap, out = [Path(os.environ[k]) for k in ['ROOT', 'SNAP', 'OUT']]
def read(p):
    with p.open(encoding='utf-8-sig', newline='') as f:
        r = csv.DictReader(f, delimiter='\t')
        return r.fieldnames, list(r)
def write(name, rows, fields):
    with (out / name).open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)
def number(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return math.nan

spec = []
mouse = '02_analysis/GSE314943/stage6_20260905_231807_L44XdT/tables'
for scope in ['Core', 'ActCD4', 'Naive', 'Treg']:
    for arm in ['Full', 'NoBG']:
        spec.append(('GSE314943', scope, arm, f'{mouse}/{scope}_{arm}_DEG.tsv'))
timebase = '02_analysis/GSE274372/stage10_20260908_141459_a8IPaq/tables'
for scope in ['Author', 'NoMix', 'TLike']:
    for day in ['D7', 'D28', 'D56']:
        spec.append(('GSE274372', scope, f'{day}_D0', f'{timebase}/{scope}_{day}_D0.tsv'))
human = '07_DEG/CD4_pseudobulk_v1/source_adjusted'
for scope in ['Core', 'Primary']:
    for arm in ['SrcAdj', 'G227', 'Unadj']:
        spec.append(('Human_integrated', scope, arm, f'{human}/{scope}_{arm}_edgeR.tsv'))

loaded = []
for dataset, scope, arm, rel in spec:
    fields, rows = read(snap / 'source' / rel)
    assert {'logFC', 'PValue', 'FDR'} <= set(fields), rel
    assert 'gene' in fields or 'gene_id' in fields, rel
    loaded.append((dataset, scope, arm, rel, rows))
assert len(loaded) == 23
print('ANALYSIS_LOGIC_PASS=23_EXPLICIT_GENE_TABLES_SCOPES_SEPARATE', flush=True)

summary, selected, previews, targets, invalid = [], [], [], [], []
target_set = {'MRE11', 'MRE11A', 'IL18', 'IL18R1', 'IL18RAP', 'CXCL16', 'CXCR6'}
for dataset, scope, arm, rel, rows in loaded:
    clean = []
    for i, r in enumerate(rows):
        gene = r.get('gene_id') or r.get('gene') or ''
        symbol = r.get('symbol') or gene
        fc, p, fdr = [number(r.get(k)) for k in ['logFC', 'PValue', 'FDR']]
        item = dict(dataset=dataset, scope=scope, arm=arm, source=rel,
                    gene_id=gene, symbol=symbol, logFC=r.get('logFC', ''),
                    logCPM=r.get('logCPM', ''), PValue=r.get('PValue', ''),
                    FDR=r.get('FDR', ''), background=r.get('background', ''),
                    original_status=r.get('status', ''),
                    original_direction=r.get('direction', ''),
                    direction='POS' if fc > 0 else 'NEG' if fc < 0 else 'ZERO',
                    housekeeping_flag=bool(re.match(r'^(RPL\d|RPS\d|MRPL\d|MRPS\d|MT-)', symbol, re.I)),
                    selected_reason='', within_direction_rank='')
        if symbol.upper() in target_set or gene.upper() in target_set:
            targets.append(item.copy())
        if not (math.isfinite(fc) and math.isfinite(p) and 0 <= p <= 1):
            invalid.append(dict(source=rel, row=i+2, gene=gene, reason='INVALID_FC_OR_P'))
            continue
        if math.isfinite(fdr) and not 0 <= fdr <= 1:
            raise ValueError('FDR outside [0,1]: ' + rel)
        clean.append((p, -abs(fc), gene, item, fc, fdr))
    clean.sort(key=lambda x: x[:3])
    eligible = [x for x in clean if x[0] < 0.05 and abs(x[4]) >= 0.5]
    for direction in ['POS', 'NEG']:
        keep = [x for x in eligible if x[3]['direction'] == direction][:20]
        for rank, x in enumerate(keep, 1):
            item = x[3].copy()
            item.update(selected_reason='P_LT_0.05_ABS_LOG2FC_GE_0.5', within_direction_rank=rank)
            selected.append(item)
    for x in clean[:10]:
        item = x[3].copy()
        item['selected_reason'] = 'RANK_PREVIEW_NOT_SIGNIFICANCE_CLAIM'
        previews.append(item)
    summary.append(dict(dataset=dataset, scope=scope, arm=arm, source=rel,
                        rows=len(rows), valid=len(clean), nominal_p05=sum(x[0]<0.05 for x in clean),
                        fdr05=sum(math.isfinite(x[5]) and x[5]<0.05 for x in clean),
                        eligible=len(eligible), eligible_pos=sum(x[4]>0 for x in eligible),
                        eligible_neg=sum(x[4]<0 for x in eligible),
                        symbols_duplicated=len(clean)-len({x[3]['symbol'] for x in clean})))

columns = ['dataset', 'scope', 'arm', 'symbol', 'gene_id', 'direction', 'logFC', 'logCPM',
           'PValue', 'FDR', 'background', 'housekeeping_flag', 'original_status',
           'original_direction', 'selected_reason', 'within_direction_rank', 'source']
write('candidate_top20_each_direction.tsv', selected, columns)
write('rank_preview_top10.tsv', previews, columns)
write('previous_targets.tsv', targets, columns)
write('comparison_summary.tsv', summary, list(summary[0]))
write('invalid_rows.tsv', invalid, ['source', 'row', 'gene', 'reason'])

# Small original scripts/design tables establish contrast direction and model.
# Oversized cell metadata is not included in this return package.
_, manifest = read(snap / 'collected_files.tsv')
context_log = []
for r in manifest:
    if r['role'] != 'CONTEXT':
        continue
    p = snap / 'source' / r['path']
    if not p.is_file() or p.stat().st_size > 128 * 1024:
        context_log.append(dict(path=r['path'], action='OMITTED_OVER_128_KIB_OR_MISSING'))
        continue
    dest = out / 'context' / r['path']
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(p, dest)
    context_log.append(dict(path=r['path'], action='COPIED'))
write('context_manifest.tsv', context_log, ['path', 'action'])
(out / 'methods.txt').write_text(
    'Existing statistics only; no new hypothesis tests. P<0.05 and abs(log2FC)>=0.5 are exploratory screening thresholds.\n'
    'Top 20 positive and negative coefficients per table ranked by nominal P, then abs(logFC). FDR is retained unchanged.\n'
    'POS/NEG are coefficient signs. Confirm biological contrast direction using original code before interpretation.\n'
    'Full/NoBG, Core/Primary, Author/NoMix/TLike and adjusted/unadjusted arms are not independent replications.\n'
    'GSE274372 time comparisons and human integrated comparisons retain their distinct scope.\n'
    'Housekeeping flags are a limited name heuristic, not an exclusion rule.\n'
    'Independent sample consistency and detection rates are not computed by this script.\n'
    'GSE183682 and GSE174725 genome-wide DE tables were not present in the supplied inventory.\n', encoding='utf-8')
print('CANDIDATE_ROWS=' + str(len(selected)), flush=True)
print('INVALID_ROWS=' + str(len(invalid)), flush=True)
archive = out.with_suffix('.tar.gz')
with tarfile.open(archive, 'w:gz') as t:
    for p in sorted(out.rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts:
            t.add(p, arcname=p.relative_to(out), recursive=False)
assert archive.stat().st_size < 3 * 1024 * 1024, 'Return archive exceeds 3 MiB: ' + str(archive)
dest = root / 'tmp/lung_CD4_candidate_slim_return.tar.gz'
temp = dest.with_name(dest.name + f'.{os.getpid()}.part')
shutil.copyfile(archive, temp)
os.replace(temp, dest)
print('RETURN_BYTES=' + str(dest.stat().st_size), flush=True)
print('RETURN_PACKAGE=' + str(dest), flush=True)
print('CANDIDATE_SLIM_PASS', flush=True)
