from pathlib import Path
import json,hashlib,zipfile,subprocess,time
import pandas as pd,pymupdf as fitz
O=Path(__file__).parent;d=json.loads((O/'report_delivery.json').read_text(encoding='utf-8'));q=json.loads((O/'layout_QA.json').read_text(encoding='utf-8'))
assert q['all_portrait'] and q['pages']==31 and q['figure_images']==22 and not any(p['bounds'] for p in q['page_checks'])
assert json.loads((O/'R_server_exit.json').read_text())['exit_code']==0
for id in ['F07','F08']:
 doc=fitz.open(O/'figures'/(id+'.pdf'));t=''.join(p.get_text() for p in doc);assert t.count('P=')==20 and 'q=' not in t and 'FDR=' not in t
q['GO_KEGG_labels']='20 raw-P annotations per figure; no q/FDR annotation';q['visual_check']='all31pages inspected; legends outside plotting areas; full names and tables visible; final page25 line-wrap repaired';(O/'layout_QA.json').write_text(json.dumps(q,ensure_ascii=False,indent=2),encoding='utf-8')
reps=json.loads((O/'editorial_replacements.json').read_text(encoding='utf-8'))
def clean(t):
 for a,b in reps:t=t.replace(a,b)
 return t
ms=json.loads((O/'figure_manifest.json').read_text(encoding='utf-8'));idx=pd.DataFrame([{'figure':m['id'],'report_page':i+5,'module':clean(m['module']),'title':clean(m['title']),'caption':clean(m['caption'])} for i,m in enumerate(ms)]);idx.to_csv(O/'图件索引.tsv',sep='\t',index=False,encoding='utf-8')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
files=sorted((O/'figures').glob('*.png'))+sorted((O/'figures').glob('*.pdf'));assert len(files)==44
manifest=[{'path':'figures/'+p.name,'sha256':sha(p)} for p in files];(O/'figure_checksums.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
z=O/'v22_高清图件.zip'
with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED) as zz:
 for p in files:zz.write(p,'figures/'+p.name)
 zz.write(O/'图件索引.tsv','图件索引.tsv');zz.write(O/'figure_checksums.json','figure_checksums.json')
with zipfile.ZipFile(z) as zz:
 assert zz.testzip() is None
 for r in manifest:assert hashlib.sha256(zz.read(r['path'])).hexdigest()==r['sha256']
d.update(figure_zip=str(z),figure_zip_sha256=sha(z),report_pages=31,all_portrait=True)
(O/'delivery.json').write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')
P='/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/client_layout_20261017_v22'
def run(cmd):
 for i in range(4):
  r=subprocess.run(cmd,capture_output=True,creationflags=subprocess.CREATE_NO_WINDOW)
  if r.returncode==0:return r.stdout.decode('utf-8')
  if r.returncode not in [1,255]:break
  time.sleep(2)
 raise RuntimeError(r.stderr.decode('utf-8'))
run(['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10','mkdir -p '+P+'/delivery'])
run(['scp','-q',d['pdf'],str(z),'sgs10:'+P+'/delivery/'])
run(['scp','-q','-r',str(O/'figures'),str(O/'plot_data'),'sgs10:'+P+'/'])
run(['scp','-q',*[str(p) for p in O.glob('*.py')],*[str(p) for p in O.glob('*.R')],'sgs10:'+P+'/code/'])
run(['scp','-q',*[str(p) for p in O.glob('*.json')],str(O/'图件索引.tsv'),'sgs10:'+P+'/'])
s=run(['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10','sha256sum '+P+'/delivery/*']);assert d['sha256'] in s and d['figure_zip_sha256'] in s
c={'status':'completed','report_pages':31,'representative_figures':22,'all_portrait':True,'source_statistics_preserved':True,'GSEA_reference_plotting_applied':True,'GO_annotations':'raw P','communication':'chord diagrams','remote_checksums_match':True,'delivery':d};(O/'completion.json').write_text(json.dumps(c,ensure_ascii=False,indent=2),encoding='utf-8');run(['scp','-q',str(O/'completion.json'),'sgs10:'+P+'/']);print(json.dumps(c,ensure_ascii=False))