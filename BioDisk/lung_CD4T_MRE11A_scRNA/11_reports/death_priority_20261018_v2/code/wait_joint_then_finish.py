import os,sys,json,select,subprocess
from pathlib import Path
import psutil
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/death_priority_20261018_v2');N=O/'joint_network';target=str(O/'code/joint_network.py')
roots=[p for p in psutil.process_iter(['pid','cmdline','uids']) if p.info['uids'].real==os.getuid() and p.info['cmdline'] and target in p.info['cmdline']]
if roots:
 assert len(roots)==1
 pid=roots[0].pid
 print('NATIVE_COMPLETION_WAIT',pid,flush=True)
 if hasattr(os,'pidfd_open'):
  fd=os.pidfd_open(pid);select.select([fd],[],[]);os.close(fd)
 else:
  subprocess.run(['tail','--pid='+str(pid),'-f','/dev/null'],check=True)
assert json.loads((N/'network_inference_complete.json').read_text())['status']=='success','network inference not complete'
commands=[
 ('joint_network_score.py','/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/software/pyscenic_0.12.1_py310/bin/python'),
 ('joint_network_stats.py','/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python'),
 ('joint_network_finish.py','/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python')]
for name,exe in commands:
 print('FOLLOWUP_START',name,flush=True)
 with (N/(name+'.log')).open('wb') as log:r=subprocess.run([exe,'-u',str(O/'code'/name)],stdout=log,stderr=subprocess.STDOUT)
 (N/(name+'.exit.json')).write_text(json.dumps({'returncode':r.returncode}),encoding='utf-8')
 print('FOLLOWUP_EXIT',name,r.returncode,flush=True)
 if r.returncode:raise RuntimeError(name+' failed; inspect stage log')
assert json.loads((N/'QA_complete.json').read_text())['status']=='success'
print('JOINT_NETWORK_AND_DOWNSTREAM_COMPLETE',flush=True)
