from pathlib import Path
import subprocess,json,tarfile,hashlib
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/death_priority_20261018_v3')
if not (O/'network_loo/statistics_complete.json').exists():
 pids=[]
 for p in Path('/proc').iterdir():
  if not p.name.isdigit():continue
  try:
   if any(a.endswith(b'/code/loo_stats.py') for a in (p/'cmdline').read_bytes().split(b'\0')):pids.append(int(p.name))
  except OSError:pass
 assert len(pids)==1,pids
 subprocess.run(['tail','--pid='+str(pids[0]),'-f','/dev/null'],check=True)
for name in ['complete.json','scoring_complete.json','statistics_complete.json']:
 assert json.loads((O/'network_loo'/name).read_text())['status']=='success'
folds=list((O/'network_loo').glob('*/*/grn_exit.json'));assert len(folds)==20
for grn in folds:
 for stage in ['grn','ctx']:
  f=grn.parent/(stage+'_exit.json');assert json.loads(f.read_text())['returncode']==0
 assert (grn.parent/'regulons.csv').stat().st_size>0
for name in ['external_DE_complete.txt','mouse_HIF1A_complete.txt','human_validation_complete.txt','human_DE_complete.txt','BAL_member_validation_complete.txt']:
 assert (O/'logs'/name).read_text().strip()=='success'
out=O/'v3_review_materials.tar.gz'
with tarfile.open(out,'w:gz') as tf:
 for folder in ['tables','logs','code','figures']:tf.add(O/folder,arcname=folder)
 for f in O.glob('*.json'):tf.add(f,arcname=f.name)
 for f in (O/'network_loo').glob('*.json'):tf.add(f,arcname='network_loo/'+f.name)
 for f in (O/'network_loo').glob('*/*/*_exit.json'):tf.add(f,arcname=str(f.relative_to(O)))
print('ALL_ANALYSES_COMPLETE',out.stat().st_size,hashlib.sha256(out.read_bytes()).hexdigest(),flush=True)
