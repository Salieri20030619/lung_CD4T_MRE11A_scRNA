from pathlib import Path
import os,json,time,shutil,psutil
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/death_priority_20261018_v3');p=psutil.Process(1705790)
assert p.uids().real==os.getuid() and any(a.endswith('/code/network_loo.py') for a in p.cmdline())
children=p.children(recursive=True);assert all(q.uids().real==os.getuid() for q in children)
snapshot=[{'pid':q.pid,'cmd':q.cmdline()} for q in [p]+children]
(O/'rescale_process_snapshot.json').write_text(json.dumps(snapshot,indent=2),encoding='utf-8')
# Stop submission first, then only this task's descendants.
p.terminate()
for q in reversed(children):
 try:q.terminate()
 except psutil.NoSuchProcess:pass
gone,alive=psutil.wait_procs([p]+children,timeout=15)
for q in alive:
 try:q.kill()
 except psutil.NoSuchProcess:pass
psutil.wait_procs(alive,timeout=10)
complete=0;interrupted=[]
for d in (O/'network_loo').glob('*/*'):
 if not d.is_dir():continue
 status={}
 for stage in ['grn','ctx']:
  f=d/(stage+'_exit.json');status[stage]=f.exists() and json.loads(f.read_text()).get('returncode')==0
 if all(status.values()):complete+=1;continue
 for stage in ['grn','ctx']:
  if status[stage]:continue
  for name in [stage+'.log',stage+'_exit.json',('adj.tsv' if stage=='grn' else 'regulons.csv')]:
   f=d/name
   if f.exists():f.rename(d/(name+'.initial8workers'))
 interrupted.append({'fold':str(d.relative_to(O)),'completed_stages':status})
shutil.copy2(O/'protocol_v3.json',O/'protocol_v3.initial8workers.json')
(O/'rescale_audit.json').write_text(json.dumps({'reason':'reuse completed checkpoints; use 8 concurrent jobs x10 workers within current256CPU shared host load and268GiB available memory; explicit BioData temporary directory','completed_folds_preserved':complete,'unfinished':interrupted,'old_parent_pid':1705790,'terminated_only_own_process_tree':True,'timestamp':time.time()},indent=2),encoding='utf-8')
assert complete>=6
print('RESCALE_CHECKPOINTS_PRESERVED',complete,flush=True)
