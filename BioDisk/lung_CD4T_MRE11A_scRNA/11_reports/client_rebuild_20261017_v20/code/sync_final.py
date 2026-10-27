from pathlib import Path
import subprocess,json,hashlib,time
O=Path(__file__).parent;d=json.loads((O/'delivery.json').read_text(encoding='utf-8'));P='/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/client_rebuild_20261017_v20'
def run(cmd):
 for i in range(4):
  r=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,creationflags=subprocess.CREATE_NO_WINDOW)
  if r.returncode==0:return r.stdout
  if r.returncode not in [1,255]:break
  time.sleep(2)
 raise RuntimeError(r.stderr.decode('utf-8'))
run(['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10','mkdir -p '+P+'/delivery'])
run(['scp','-q',d['pdf'],d['support_zip'],'sgs10:'+P+'/delivery/'])
run(['scp','-q',*[str(p) for p in O.glob('*.py')],str(O/'enrichment_v20.R'),'sgs10:'+P+'/code/'])
run(['scp','-q','-r',str(O/'tables'),str(O/'plot_data'),'sgs10:'+P+'/'])
run(['scp','-q',*[str(O/n) for n in ['plan.json','review_scope.json','delivery.json','report_QA.json','report_page_index.json','figure_manifest.json']],'sgs10:'+P+'/'])
r=run(['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10','sha256sum '+P+'/delivery/*']).decode('utf-8');assert d['sha256'] in r and d['support_sha256'] in r
complete={'status':'completed','report_pages':31,'representative_figures':22,'modules':11,'predelivery_checks':'completed within documented scope','independent_review':'v19 previously reviewed; v20 added analyses not independently reviewed','remote_checksums_match':True,'remote_directory':P,'delivery':d};(O/'completion.json').write_text(json.dumps(complete,ensure_ascii=False,indent=2),encoding='utf-8');run(['scp','-q',str(O/'completion.json'),'sgs10:'+P+'/']);print(json.dumps(complete,ensure_ascii=False))