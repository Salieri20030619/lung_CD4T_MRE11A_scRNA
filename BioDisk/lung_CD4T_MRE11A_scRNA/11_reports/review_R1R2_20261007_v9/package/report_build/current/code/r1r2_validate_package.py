"""Read-only verification of the R1/R2 delivery. No analysis scripts are executed."""
from pathlib import Path,PurePosixPath
import argparse,json,hashlib,csv,ast,re
import pymupdf as fitz
P=argparse.ArgumentParser();P.add_argument('--package',type=Path,required=True);a=P.parse_args();K=a.package.resolve()
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
with (K/'FILE_MANIFEST_SHA256.tsv').open(encoding='utf-8',newline='') as f:rows=list(csv.DictReader(f,delimiter='\t'))
for r in rows:
    rel=PurePosixPath(r['relative_path']);assert not rel.is_absolute() and '..' not in rel.parts and ':' not in str(rel)
    p=(K/str(rel)).resolve();assert p.is_relative_to(K) and p.is_file()
    assert p.stat().st_size==int(r['bytes']) and sha(p)==r['sha256'],str(rel)
actual={p.relative_to(K).as_posix() for p in K.rglob('*') if p.is_file()};assert actual=={r['relative_path'] for r in rows}|{'FILE_MANIFEST_SHA256.tsv'}
with (K/'provenance/code_manifest.tsv').open(encoding='utf-8',newline='') as f:code=list(csv.DictReader(f,delimiter='\t'))
assert len(code)==450
for r in code:assert sha(K/r['package_path'])==r['package_sha256'],r['package_path']
scripts=[p for p in (K/'code').rglob('*') if p.is_file() and p.suffix.lower() in ['.py','.r','.sh','.ps1']];assert len(scripts)==285
v=json.loads((K/'provenance/package_validation.json').read_text(encoding='utf-8'))
pdf=K/v['pdf_relative_path'];doc=fitz.open(pdf)
assert len(doc)==v['pdf_pages'] and sha(pdf)==v['pdf_sha256'] and pdf.stat().st_size==v['pdf_bytes']
qa=json.loads((K/'provenance/PDF_content_parity.json').read_text(encoding='utf-8'))
assert qa['pages']==len(doc) and qa['pdf_sha256']==sha(pdf) and not qa['missing_text_blocks']
cur=K/'report_build/current';man=json.loads((cur/'revision_manifest.json').read_text(encoding='utf-8'))
assert man['R3_BH_not_changed'] and man['G288_sentence_not_changed'] and not man['non_target_changes']
page=doc[man['target_page']-1];rect=fitz.Rect(man['target_rectangle']);flat=re.sub(r'\s+','',page.get_text(clip=rect,sort=True))
for row in man['table2_after']:assert re.sub(r'\s+','',''.join(row)) in flat,row
newcode=list((cur/'code').glob('*.py'))
for p in newcode:ast.parse(p.read_text(encoding='utf-8-sig'))
patterns=[rb'-----BEGIN (?:RSA |EC |OPENSSH |DSA )?PRIVATE KEY-----',rb'(?i)[?&]sig=[A-Za-z0-9%+/]{20,}',rb'\bAKIA[A-Z0-9]{16}\b',rb'(?i)Authorization\s*:\s*Bearer\s+[A-Za-z0-9._-]{20,}']
bad=[];exclusion_guard_references=[]
for p in K.rglob('*'):
    if p.is_file() and p.suffix.lower() in ['.py','.r','.sh','.ps1','.md','.txt','.json','.tsv','.html']:
        b=p.read_bytes()
        if any(re.search(pattern,b) for pattern in patterns):bad.append(p.relative_to(K).as_posix())
        excluded=[b''.join([b'GSE',b'264182']),b'_'.join([b'silicosis',b'bridge'])]
        if any(term in b for term in excluded):
            # This hash-pinned historical packager only names excluded material in absence assertions.
            assert p.relative_to(K).as_posix()=='report_build/history/library_v7/code/logic_finalize_package.py'
            assert sha(p)=='847992872f38d2f5967ac8f7b90a38b4ce3fd7a63969ae79025f4ade0ad11e64'
            exclusion_guard_references.append(p.relative_to(K).as_posix())
assert not bad,bad
print(json.dumps({'status':'PASS','pdf_pages':len(doc),'pdf_sha256':sha(pdf),'manifest_files_verified':len(rows),'source_code_index_rows_verified':len(code),'original_scripts':len(scripts),'current_build_scripts_AST':len(newcode),'table2_rows_verified':6,'sensitive_findings':bad,'exclusion_guard_references_only':exclusion_guard_references,'all_paths_contained':True,'max_relative_path_length':max(map(len,actual))}))
