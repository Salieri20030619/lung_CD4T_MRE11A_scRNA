from pathlib import Path
import json,hashlib,zipfile,shutil,subprocess
import pandas as pd
O=Path(__file__).parent;K=O/'review_bundle';d=json.loads((O/'delivery.json').read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
q=json.loads((O/'pdf_QA_v19.json').read_text(encoding='utf-8'));q['visual_contact_sheet_checked']=True;q['visual_result']='no clipping, overlaps, missing glyphs or pagination defects observed';(O/'pdf_QA_v19.json').write_text(json.dumps(q,ensure_ascii=False,indent=2),encoding='utf-8')
shutil.copy2(O/'pdf_QA_v19.json',K/'pdf_QA_v19.json');shutil.copy2(O/'verify_delivery.py',K/'code/v19/verify_delivery.py')
# Check model file pointers supplied by the DE trace.
t=pd.read_csv(K/'R01_candidate_DE_all_tested_records.tsv',sep='\t');col='sample_design_file'
if col in t:
 for n in t[col].dropna().unique():assert (K/n).exists(),n
files=[p for p in K.rglob('*') if p.is_file() and p.name!='MANIFEST_SHA256.tsv'];pd.DataFrame([{'path':p.relative_to(K).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)} for p in files]).to_csv(K/'MANIFEST_SHA256.tsv',sep='\t',index=False,encoding='utf-8')
z=Path(d['zip'])
with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as f:
 for p in K.rglob('*'):
  if p.is_file():f.write(p,p.relative_to(K).as_posix())
with zipfile.ZipFile(z) as f:
 assert f.testzip() is None
 for p in files:assert hashlib.sha256(f.read(p.relative_to(K).as_posix())).hexdigest()==sha(p)
d.update(zip_sha256=sha(z),files=len(files)+1,zip_MB=z.stat().st_size/1e6);(O/'delivery.json').write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')
P='/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/feedback_revision_20261017_v19';kw={'check':True,'creationflags':subprocess.CREATE_NO_WINDOW}
subprocess.run(['ssh','-T','sgs10','mkdir -p '+P+'/delivery'],**kw)
subprocess.run(['scp','-q','-r',str(O/'tables'),'sgs10:'+P+'/'],**kw)
subprocess.run(['scp','-q',d['pdf'],d['zip'],'sgs10:'+P+'/delivery/'],**kw)
for p in [O/'QA_v19.json',O/'pdf_QA_v19.json',O/'delivery.json',O/'reconstruction_QA.json',O/'scope.json']:
 subprocess.run(['scp','-q',str(p),'sgs10:'+P+'/'],**kw)
subprocess.run(['scp','-q',*[str(p) for p in O.glob('*.py')],'sgs10:'+P+'/code/'],**kw)
r=subprocess.run(['ssh','-T','sgs10','sha256sum '+P+'/delivery/*'],stdout=subprocess.PIPE,**kw).stdout.decode('utf-8');assert d['pdf_sha256'] in r and d['zip_sha256'] in r
record={'status':'targeted revision completed; independent review pending','delivery':d,'remote_checksums_match':True,'remote_directory':P};(O/'completion.json').write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8');subprocess.run(['scp','-q',str(O/'completion.json'),'sgs10:'+P+'/'],**kw);print(json.dumps(record,ensure_ascii=False))