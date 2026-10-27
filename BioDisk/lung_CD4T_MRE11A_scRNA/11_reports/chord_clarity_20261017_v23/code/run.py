from pathlib import Path
import subprocess,time,json
O=Path(__file__).parent;P='/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/chord_clarity_20261017_v23'
def run(cmd):
 for i in range(4):
  r=subprocess.run(cmd,capture_output=True,creationflags=subprocess.CREATE_NO_WINDOW)
  if r.returncode==0:return r
  if r.returncode!=255:break
  time.sleep(2)
 print(r.stderr.decode('utf-8'));r.check_returncode()
run(['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10','mkdir -p '+P+'/code'])
run(['scp','-q',str(O/'chord_clarity.R'),'sgs10:'+P+'/code/'])
r=run(['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10','/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript --vanilla '+P+'/code/chord_clarity.R']);(O/'R.log').write_bytes(r.stdout+r.stderr)
run(['scp','-q','-r','sgs10:'+P+'/figures','sgs10:'+P+'/plot_data',str(O)+'/']);(O/'R_exit.json').write_text(json.dumps({'exit_code':r.returncode}),encoding='utf-8');print('CHORDS_COMPLETE')