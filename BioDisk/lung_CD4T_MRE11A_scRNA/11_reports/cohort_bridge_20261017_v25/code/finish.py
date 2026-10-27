from pathlib import Path
import json,hashlib,zipfile,subprocess,time,pandas as pd
O=Path(__file__).parent;B=O.parent/'lung_screening_display_20261017_v24';q=json.loads((O/'QA.json').read_text(encoding='utf-8'));assert q['pages']==37 and q['total_images']==27 and not any(c['bounds'] for c in q['page_checks']);q['visual_check']='Modified pages2,3,4,7,9,11,12,13 visually checked: no clipping/overlapping legends; source cohorts identified before joint analysis';(O/'QA.json').write_text(json.dumps(q,ensure_ascii=False,indent=2),encoding='utf-8')
d=json.loads((O/'report_delivery.json').read_text(encoding='utf-8'));ms=json.loads((O/'figure_manifest.json').read_text(encoding='utf-8'));extra=json.loads((O/'extra_manifest.json').read_text(encoding='utf-8'));heads=json.loads((O/'report_page_index.json').read_text(encoding='utf-8'));reps=json.loads((O/'editorial_replacements.json').read_text(encoding='utf-8'))
def clean(s):
 for a,b in reps:s=s.replace(a,b)
 return s
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
index=[];files=[]
for m in ms:
 title=clean(m['module']+'｜'+m['title']);matches=[h['page'] for h in heads if h['heading']==title];assert len(matches)==1
 index.append({'id':m['id'],'page':matches[0],'title':clean(m['title']),'caption':clean(m['caption'])})
 for ext in ['png','pdf']:files.append((Path(m['file']).with_suffix('.'+ext),'figures/'+m['id']+'.'+ext))
for e,n in zip(extra,[3,4,11]):
 index.append({'id':e['id'],'page':n,'title':e['title'],'caption':e['role']})
 for ext in ['png','pdf']:files.append((O/'figures'/(e['id']+'.'+ext),'figures/'+e['id']+'.'+ext))
pd.DataFrame(index).sort_values('page').to_csv(O/'图件索引.tsv',sep='\t',index=False,encoding='utf-8');files.append((O/'图件索引.tsv','图件索引.tsv'))
for folder,prefix in [(O/'tables','cohort_tables'),(O/'plot_data','new_plot_data'),(B/'tables','screening_tables')]:
 for p in folder.iterdir():
  if p.is_file():files.append((p,prefix+'/'+p.name))
z=O/'v25_图件与队列筛选附表.zip';checks=[]
with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED) as zz:
 for p,key in files:zz.write(p,key);checks.append({'file':key,'sha256':sha(p)})
 zz.writestr('checksums.json',json.dumps(checks,indent=2))
with zipfile.ZipFile(z) as zz:
 assert zz.testzip() is None
 for r in checks:assert hashlib.sha256(zz.read(r['file'])).hexdigest()==r['sha256']
d.update(pages=37,figure_zip=str(z),figure_zip_sha256=sha(z));(O/'delivery.json').write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')
P='/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/cohort_bridge_20261017_v25'
def run(cmd):
 for i in range(4):
  r=subprocess.run(cmd,capture_output=True,creationflags=subprocess.CREATE_NO_WINDOW)
  if r.returncode==0:return r.stdout.decode('utf-8')
  if r.returncode not in [1,255]:break
  time.sleep(2)
 raise RuntimeError(r.stderr.decode('utf-8'))
run(['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10','mkdir -p '+P+'/delivery'])
run(['scp','-q',d['pdf'],str(z),'sgs10:'+P+'/delivery/']);run(['scp','-q','-r',str(O/'figures'),str(O/'plot_data'),'sgs10:'+P+'/']);run(['scp','-q',*[str(p) for p in O.glob('*.py')],'sgs10:'+P+'/code/']);run(['scp','-q',*[str(p) for p in O.glob('*.json')],str(O/'report_text.txt'),'sgs10:'+P+'/'])
s=run(['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10','sha256sum '+P+'/delivery/*']);assert d['sha256'] in s and d['figure_zip_sha256'] in s
c={'status':'complete','report_pages':37,'numbered_figures':24,'front_displays':3,'joint10':{'GSE314943':6,'GSE183682':4},'supplement_CD4_cells':1413,'front_pages':[3,4,11],'statistics_unchanged':True,'remote_checksum_match':True,'delivery':d};(O/'completion.json').write_text(json.dumps(c,ensure_ascii=False,indent=2),encoding='utf-8');run(['scp','-q',str(O/'completion.json'),'sgs10:'+P+'/']);print(json.dumps(c,ensure_ascii=False))