from pathlib import Path
import subprocess,time
O=Path(__file__).parent;P='/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/chord_clarity_20261017_v23'
def run(cmd):
 for i in range(4):
  r=subprocess.run(cmd,capture_output=True,creationflags=subprocess.CREATE_NO_WINDOW)
  if r.returncode==0:return r
  if r.returncode!=255:break
  time.sleep(2)
 print(r.stderr.decode('utf-8'));r.check_returncode()
run(['scp','-q',str(O/'trajectory_backbone.py'),'sgs10:'+P+'/code/']);r=run(['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10','/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python -u '+P+'/code/trajectory_backbone.py']);(O/'trajectory.log').write_bytes(r.stdout+r.stderr);run(['scp','-q','sgs10:'+P+'/plot_data/trajectory_backbone*',str(O/'plot_data')+'/']);run(['scp','-q','sgs10:'+P+'/trajectory_backbone_QA.json',str(O)+'/']);print(r.stdout.decode('utf-8'))