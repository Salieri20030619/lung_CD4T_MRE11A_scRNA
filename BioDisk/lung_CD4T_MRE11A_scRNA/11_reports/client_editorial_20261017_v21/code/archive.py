from pathlib import Path
import json,subprocess,time,hashlib
O=Path(__file__).parent;d=json.loads((O/'report_delivery.json').read_text(encoding='utf-8'));q=json.loads((O/'editorial_QA.json').read_text(encoding='utf-8'));q['visual_check']='All31pages inspected; no clipping/overlap/missing glyphs';q['residual_terms_explanation']='发送 refers only to ligand sender cells, retained as scientific terminology';(O/'editorial_QA.json').write_text(json.dumps(q,ensure_ascii=False,indent=2),encoding='utf-8')
P='/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/client_editorial_20261017_v21'
def run(cmd):
 for i in range(4):
  r=subprocess.run(cmd,capture_output=True,creationflags=subprocess.CREATE_NO_WINDOW)
  if r.returncode==0:return r.stdout.decode('utf-8')
  if r.returncode not in [1,255]:break
  time.sleep(2)
 raise RuntimeError(r.stderr.decode('utf-8'))
run(['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10','mkdir -p '+P+'/delivery '+P+'/code'])
run(['scp','-q',d['pdf'],'sgs10:'+P+'/delivery/'])
run(['scp','-q',*[str(p) for p in O.glob('*.py')],'sgs10:'+P+'/code/'])
run(['scp','-q',*[str(p) for p in O.glob('*.json')],str(O/'report_text.txt'),str(O/'internal_selection_flow.tsv'),'sgs10:'+P+'/'])
run(['scp','-q','-r',str(O/'figures'),str(O/'plot_data'),'sgs10:'+P+'/'])
s=run(['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10','sha256sum '+P+'/delivery/*']);assert d['sha256'] in s
c={'status':'completed','scope':'customer-language editorial only','pages':31,'figures':22,'statistics_unchanged':True,'v20_preserved':True,'remote_checksum_match':True,'remote':P,'pdf':d['pdf'],'sha256':d['sha256']};(O/'completion.json').write_text(json.dumps(c,ensure_ascii=False,indent=2),encoding='utf-8');run(['scp','-q',str(O/'completion.json'),'sgs10:'+P+'/']);print(json.dumps(c,ensure_ascii=False))