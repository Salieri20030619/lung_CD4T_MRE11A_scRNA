from pathlib import Path
import subprocess,json,time
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/four_axis_screen_20261017_v15');t=time.time()
with (O/'provenance/screen_run.log').open('w',encoding='utf-8') as log:
 p=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python','-u',str(O/'code/resume_screen.py')],stdout=log,stderr=subprocess.STDOUT)
s={'exit_code':p.returncode,'elapsed_seconds':time.time()-t};(O/'provenance/screen_exit.json').write_text(json.dumps(s,indent=2),encoding='utf-8');print(json.dumps(s));print((O/'provenance/screen_run.log').read_text()[-9000:])