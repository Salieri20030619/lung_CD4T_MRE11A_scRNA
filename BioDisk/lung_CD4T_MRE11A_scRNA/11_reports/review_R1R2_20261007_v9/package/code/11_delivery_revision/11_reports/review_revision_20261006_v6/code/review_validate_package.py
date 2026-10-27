"""Validate and archive the curated package; never execute scientific analyses."""
from pathlib import Path,PurePosixPath
import csv,json,hashlib,zipfile,re,collections,ast,os
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
O=P/'11_reports/review_revision_20261006_v6';C=O/'package'
manifest=list(csv.DictReader((C/'provenance/code_manifest.tsv').open(encoding='utf-8'),delimiter='\t'))
expected=set()
for r in manifest:
 f=C/r['package_path'];expected.add(f.resolve());assert f.is_file() and not f.is_symlink()
 assert hashlib.sha256(f.read_bytes()).hexdigest()==r['package_sha256']==r['original_sha256']
 assert hashlib.sha256(Path(r['original_path']).read_bytes()).hexdigest()==r['original_sha256']
actual={x.resolve() for x in (C/'code').rglob('*') if x.is_file()}
assert actual==expected,{'unlisted_code':list(map(str,actual-expected)),'missing':list(map(str,expected-actual))}
checks=list(csv.DictReader((C/'provenance/static_checks.tsv').open(encoding='utf-8'),delimiter='\t'))
rs=list(csv.DictReader((C/'provenance/r_parse_checks.tsv').open(encoding='utf-8'),delimiter='\t'))
assert all(r['status']=='PASS' for r in checks+rs)
ps=json.loads((C/'provenance/powershell_parse_checks.json').read_text(encoding='utf-8-sig'))
assert len(ps)==4 and all(x['status']=='PASS' for x in ps)
qa=json.loads((C/'provenance/QA_result.json').read_text(encoding='utf-8'))
parity=json.loads((C/'provenance/PDF_content_parity.json').read_text(encoding='utf-8'))
rev=json.loads((C/'provenance/revision_manifest.json').read_text(encoding='utf-8'))
num=json.loads((C/'provenance/numeric_revision_checks.json').read_text(encoding='utf-8'))
assert qa['pages']==29 and qa['tables']==18 and qa['figures']==14 and not qa['out_of_bounds'] and qa['visual_review']['all_pages_reviewed']==29
assert not parity['missing_text_blocks'] and parity['docx_text_blocks_verified']==1643
assert rev['original_tables_preserved']==15 and rev['figure_bytes_identical'] and rev['original_decimal_values_preserved']
assert num['aggregate_rows_checked']==19 and num['baseline_endpoint_checked'] and num['old_word_unchanged']
pdfs=list((C/'report').glob('*.pdf'));assert len(pdfs)==1
pdf=pdfs[0];sha=hashlib.sha256(pdf.read_bytes()).hexdigest();assert sha==qa['sha256']['review_report.pdf']
for f in C.rglob('*'):
 if not f.is_file():continue
 assert not f.is_symlink() and f.resolve().is_relative_to(C.resolve())
 assert f.suffix.lower() not in ['.rds','.h5ad','.loom','.h5','.hdf5','.mtx','.gz','.docx','.doc','.log','.zip','.feather','.rdata']
 assert not any(x in f.name.lower() for x in ['id_rsa','id_ed25519','credentials','azcopy','upload_session'])
 assert f.stat().st_size<1000000 or f==pdf
# Reuse exactly the scanner definitions, without importing/running the packaging script.
module=ast.parse((O/'code/review_package_code.py').read_text(encoding='utf-8'))
node=next(x for x in module.body if isinstance(x,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='patterns' for t in x.targets))
env={'re':re};exec(compile(ast.Module(body=[node],type_ignores=[]),'<scanner patterns>','exec'),env)
findings=[]
for f in C.rglob('*'):
 if not f.is_file() or f.suffix=='.pdf':continue
 s=f.read_text(encoding='utf-8-sig')
 for kind,rx in env['patterns'].items():
  for m in rx.finditer(s):
   if kind=='literal_secret' and any(x in m.group(1).lower() for x in ['environ','placeholder','your_','getenv','none','example','replace']):continue
   findings.append({'file':str(f.relative_to(C)),'kind':kind,'line':s[:m.start()].count('\n')+1})
assert not findings,findings
summary=json.loads((C/'provenance/package_summary.json').read_text(encoding='utf-8'))
out={'pdf_pages':29,'pdf_sha256':sha,'pdf_bytes':pdf.stat().st_size,'code_files':summary['unique_code_files'],'deduplicated_aliases':summary['duplicate_aliases'],'category_code_counts':summary['category_code_counts'],'original_sources_hash_verified':len(manifest),'python_shell_static_checks':len(checks),'R_parse_checks':len(rs),'sensitive_pattern_findings':findings,'individual_data_tables_or_large_objects_included':False,'old_word_unchanged':True,'new_scientific_models':0,'all_pdf_pages_reviewed':29,'remaining_reproduction_gaps':['GSE314943 early stage1-5 standalone execution scripts not located in project source inventory','External raw data, frozen objects and model databases not bundled','Animal-library mapping not closed'],'delivery_blockers':[]}
out['PowerShell_parse_checks']=len(ps)
(C/'provenance/package_validation.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
files=sorted(f for f in C.rglob('*') if f.is_file() and f.name!='FILE_MANIFEST_SHA256.tsv')
with (C/'FILE_MANIFEST_SHA256.tsv').open('w',encoding='utf-8',newline='') as h:
 w=csv.writer(h,delimiter='\t');w.writerow(['relative_path','bytes','sha256'])
 for f in files:w.writerow([str(f.relative_to(C)),f.stat().st_size,hashlib.sha256(f.read_bytes()).hexdigest()])
z=O/'CWP_Silica_CD4_IKZF2_STAT1_20261006_review_PDF_code.zip'
with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as h:
 for f in sorted(C.rglob('*')):
  if f.is_file():h.write(f,str(f.relative_to(C)))
with zipfile.ZipFile(z) as h:
 assert h.testzip() is None
 names=h.namelist();assert len(names)==len(set(n.casefold() for n in names))
 for name in names:
  pp=PurePosixPath(name);assert not pp.is_absolute() and '..' not in pp.parts and '\\' not in name and ':' not in name
 assert hashlib.sha256(h.read(str(pdf.relative_to(C)))).hexdigest()==sha
 out.update(zip_files=len(names),zip_max_relative_path=max(map(len,names)),zip_bytes=z.stat().st_size,zip_sha256=hashlib.sha256(z.read_bytes()).hexdigest(),zip_path=str(z))
(O/'FINAL_PACKAGE_STATUS.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(out,ensure_ascii=True))
