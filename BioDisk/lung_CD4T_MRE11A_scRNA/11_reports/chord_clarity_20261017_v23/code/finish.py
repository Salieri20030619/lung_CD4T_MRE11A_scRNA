from pathlib import Path
import json,hashlib,zipfile,subprocess,time
import pandas as pd
O=Path(__file__).parent;d=json.loads((O/'report_delivery.json').read_text(encoding='utf-8'));q=json.loads((O/'QA.json').read_text());assert q['pages']==31 and q['figures']==22 and not q['bounds'];q['visual_check']='Summary, methods, chords, trajectory and matrix pages inspected; no occlusion/overflow observed';(O/'QA.json').write_text(json.dumps(q,ensure_ascii=False,indent=2),encoding='utf-8')
ms=json.loads((O/'figure_manifest.json').read_text(encoding='utf-8'));reps=json.loads((O/'editorial_replacements.json').read_text(encoding='utf-8'))
def clean(s):
 for a,b in reps:s=s.replace(a,b)
 return s
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
idx=pd.DataFrame([{'figure':m['id'],'page':i+5,'title':clean(m['title']),'caption':clean(m['caption'])} for i,m in enumerate(ms)]);idx.to_csv(O/'图件索引.tsv',sep='\t',index=False,encoding='utf-8');z=O/'v23_高清图件.zip';checks=[]
with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED) as zz:
 for m in ms:
  for ext in ['png','pdf']:
   p=Path(m['file']).with_suffix('.'+ext);key='figures/'+m['id']+'.'+ext;zz.write(p,key);checks.append({'file':key,'sha256':sha(p)})
 zz.write(O/'图件索引.tsv','图件索引.tsv');zz.writestr('checksums.json',json.dumps(checks,indent=2))
with zipfile.ZipFile(z) as zz:
 assert zz.testzip() is None
 for c in checks:assert hashlib.sha256(zz.read(c['file'])).hexdigest()==c['sha256']
d.update(figure_zip=str(z),figure_zip_sha256=sha(z));(O/'delivery.json').write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8');P='/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/chord_clarity_20261017_v23'
def run(cmd):
 for i in range(4):
  r=subprocess.run(cmd,capture_output=True,creationflags=subprocess.CREATE_NO_WINDOW)
  if r.returncode==0:return r.stdout.decode('utf-8')
  if r.returncode not in [1,255]:break
  time.sleep(2)
 raise RuntimeError(r.stderr.decode('utf-8'))
run(['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10','mkdir -p '+P+'/delivery'])
run(['scp','-q',d['pdf'],str(z),'sgs10:'+P+'/delivery/']);run(['scp','-q','-r',str(O/'figures'),str(O/'plot_data'),'sgs10:'+P+'/'])
run(['scp','-q',*[str(p) for p in O.glob('*.py')],str(O/'chord_clarity.R'),'sgs10:'+P+'/code/']);run(['scp','-q',*[str(p) for p in O.glob('*.json')],str(O/'report_text.txt'),'sgs10:'+P+'/'])
s=run(['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10','sha256sum '+P+'/delivery/*']);assert d['sha256'] in s and d['figure_zip_sha256'] in s
record={'status':'complete','pages':31,'figures':22,'updates':['summary rewritten with purpose/data/methods/findings/validation','colored directional chords','DPT neighborhood backbone','NicheNet evidence matrix'],'original_statistical_results_unchanged':True,'remote_checksum_match':True,'delivery':d};(O/'completion.json').write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8');run(['scp','-q',str(O/'completion.json'),'sgs10:'+P+'/']);print(json.dumps(record,ensure_ascii=False))