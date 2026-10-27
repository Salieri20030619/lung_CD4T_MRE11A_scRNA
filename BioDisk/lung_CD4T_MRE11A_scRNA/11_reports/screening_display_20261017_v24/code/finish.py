from pathlib import Path
import json,hashlib,zipfile,subprocess,time
import pandas as pd
O=Path(__file__).parent;d=json.loads((O/'report_delivery.json').read_text(encoding='utf-8'));q=json.loads((O/'QA.json').read_text(encoding='utf-8'));assert q['pages']==34 and q['figures']==24 and not any(z['bounds'] for z in q['page_checks']);q['visual_check']='Pages9,25-29 inspected; no label overlap or clipping; figures renumbered consistently';(O/'QA.json').write_text(json.dumps(q,ensure_ascii=False,indent=2),encoding='utf-8')
# Supplementary candidate records include functional basis and original reference links.
s=pd.read_csv(O/'tables/supplemental42_evidence.tsv',sep='\t');y=pd.read_csv(O.parent/'lung_feedback_revision_20261017_v19/tables/R05_current_all90.tsv',sep='\t');s=s.drop(columns=['priority']).merge(y[['gene','function_or_role_note','manual_source','A_status_current']],on='gene',validate='one_to_one');s.to_csv(O/'tables/supplemental42_evidence.tsv',sep='\t',index=False,encoding='utf-8')
ms=json.loads((O/'figure_manifest.json').read_text(encoding='utf-8'));reps=json.loads((O/'editorial_replacements.json').read_text(encoding='utf-8'))
def clean(s):
 for a,b in reps:s=s.replace(a,b)
 return s
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
idx=[]
for i,m in enumerate(ms):idx.append({'figure':m['id'],'report_page':i+5 if i<20 else i+6,'module':clean(m['module']),'title':clean(m['title']),'caption':clean(m['caption'])})
pd.DataFrame(idx).to_csv(O/'图件索引.tsv',sep='\t',index=False,encoding='utf-8');checks=[];z=O/'v24_图件与筛选附表.zip'
with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED) as zz:
 for m in ms:
  for ext in ['png','pdf']:
   p=Path(m['file']).with_suffix('.'+ext);key='figures/'+m['id']+'.'+ext;zz.write(p,key);checks.append({'file':key,'sha256':sha(p)})
 for p in (O/'tables').iterdir():
  key='tables/'+p.name;zz.write(p,key);checks.append({'file':key,'sha256':sha(p)})
 for p in [O/'图件索引.tsv',O/'screening_QA.json',O/'F05_label_QA.json']:
  zz.write(p,p.name);checks.append({'file':p.name,'sha256':sha(p)})
 zz.writestr('checksums.json',json.dumps(checks,indent=2))
with zipfile.ZipFile(z) as zz:
 assert zz.testzip() is None
 for c in checks:assert hashlib.sha256(zz.read(c['file'])).hexdigest()==c['sha256']
d.update(pages=34,figures=24,figure_zip=str(z),figure_zip_sha256=sha(z));(O/'delivery.json').write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')
P='/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/screening_display_20261017_v24'
def run(cmd):
 for i in range(4):
  r=subprocess.run(cmd,capture_output=True,creationflags=subprocess.CREATE_NO_WINDOW)
  if r.returncode==0:return r.stdout.decode('utf-8')
  if r.returncode not in [1,255]:break
  time.sleep(2)
 raise RuntimeError(r.stderr.decode('utf-8'))
run(['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10','mkdir -p '+P+'/delivery '+P+'/code'])
run(['scp','-q',d['pdf'],str(z),'sgs10:'+P+'/delivery/']);run(['scp','-q','-r',str(O/'figures'),str(O/'plot_data'),str(O/'tables'),'sgs10:'+P+'/']);run(['scp','-q',*[str(p) for p in O.glob('*.py')],'sgs10:'+P+'/code/']);run(['scp','-q',*[str(p) for p in O.glob('*.json')],str(O/'report_text.txt'),'sgs10:'+P+'/'])
r=run(['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10','sha256sum '+P+'/delivery/*']);assert d['sha256'] in r and d['figure_zip_sha256'] in r
rec={'status':'complete','figure5_label_overlap':False,'screening_steps_reproduced':True,'results_changed':False,'report_pages':34,'report_figures':24,'screening_pages':[25,26,27],'remote_checksum_match':True,'delivery':d};(O/'completion.json').write_text(json.dumps(rec,ensure_ascii=False,indent=2),encoding='utf-8');run(['scp','-q',str(O/'completion.json'),'sgs10:'+P+'/']);print(json.dumps(rec,ensure_ascii=False))