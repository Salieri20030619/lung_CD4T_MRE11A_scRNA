from pathlib import Path
import json,tarfile,hashlib,shutil,datetime
import pandas as pd
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/death_priority_20261018_v2');N=O/'joint_network'
j=json.loads((N/'QA_complete.json').read_text(encoding='utf-8'));assert j['status']=='success'
for name in ['numerical_QA_complete.json','module_member_QA_complete.json','communication_export_complete.json']:assert json.loads((O/'logs'/name).read_text(encoding='utf-8'))['status']=='success'
inv=pd.read_csv(O/'tables/completion_inventory.tsv',sep='\t',dtype=str)
extra=pd.DataFrame([
 {'id':'16','analysis':'联合10样本重新推断TF网络','sample_scope':'5+5','status':'新增完成','details':f'All/BG1 each3seeds; stable scored regulons {j["All_scored"]}/{j["BG1_scored"]}','evidence':'joint_network_*'},
 {'id':'17','analysis':'评分成员逐一移除敏感性','sample_scope':'5+5','status':'新增完成','details':'5relationships;660tests;20families','evidence':'module_member_sensitivity*.tsv'}
])
inv=inv[~inv.id.isin(['16','17'])];pd.concat([inv,extra],ignore_index=True).to_csv(O/'tables/completion_inventory.tsv',sep='\t',index=False)
s=json.loads((O/'summary.json').read_text(encoding='utf-8'));s['joint_network']=j;s['module_member_sensitivity']={'tests':660,'families':20};s['analysis_complete']=True
(O/'summary.json').write_text(json.dumps(s,ensure_ascii=False,indent=2),encoding='utf-8')
old=O.parent.parent/'02_analysis/scenic1_20260929_103513_xRoEnr/pipeline.py';shutil.copy2(old,N/'original_pipeline_source.py')
records=[]
for p in list(N.glob('Mouse/*/seed*/adj.tsv'))+list(N.glob('Mouse/*/seed*/regulons.csv'))+list(N.glob('Mouse/*/score.loom'))+list(N.glob('Mouse/*/clean_stable_cell_AUC.tsv.gz'))+[O/'objects/GSE183682_whole_clustered.h5ad']:
 h=hashlib.sha256()
 with p.open('rb') as f:
  for ch in iter(lambda:f.read(8*1024*1024),b''):h.update(ch)
 records.append({'path':str(p),'bytes':p.stat().st_size,'sha256':h.hexdigest()})
(N/'large_output_manifest.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
(O/'analysis_complete.json').write_text(json.dumps({'status':'complete','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'joint_network_runs':6,'numerical_checks_original':52,'module_sensitivity_families':20,'joint_network_QA':j,'report_delivery':'pending final report integration and visual QA'},ensure_ascii=False,indent=2),encoding='utf-8')
dest=O/'review_data_v2_final.tar.gz'
with tarfile.open(dest,'w:gz') as tar:
 for folder in ['code','tables','inputs','logs','figures']:
  for p in (O/folder).rglob('*'):
   if p.is_file():tar.add(p,arcname=str(p.relative_to(O)))
 for p in N.rglob('*'):
  if not p.is_file() or p.suffix=='.loom' or p.name in ['adj.tsv','clean_stable_cell_AUC.tsv.gz']:continue
  tar.add(p,arcname=str(p.relative_to(O)))
 for f in ['protocol.json','summary.json','analysis_complete.json']:tar.add(O/f,arcname=f)
print('FINAL_DATA_BUNDLE',dest.stat().st_size,hashlib.sha256(dest.read_bytes()).hexdigest())
