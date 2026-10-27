from pathlib import Path
import subprocess,sys
O=Path(__file__).parent;P='/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/client_rebuild_20261017_v20';f=O/sys.argv[1]
def run(cmd):
 for i in range(4):
  r=subprocess.run(cmd,capture_output=True,creationflags=subprocess.CREATE_NO_WINDOW)
  if r.returncode==0:return r
  if r.returncode!=255:sys.stdout.buffer.write(r.stdout+r.stderr);r.check_returncode()
 sys.stdout.buffer.write(r.stdout+r.stderr);r.check_returncode()
run(['scp','-q',str(f),'sgs10:'+P+'/code/'+f.name]);exe='/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript --vanilla' if f.suffix=='.R' else '/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python -u';r=run(['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10',exe+' '+P+'/code/'+f.name]);(O/(f.stem+'.log')).write_bytes(r.stdout+r.stderr);sys.stdout.buffer.write(r.stdout+r.stderr)