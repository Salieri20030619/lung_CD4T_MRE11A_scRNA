from pathlib import Path
import subprocess
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/four_axis_screen_20261017_v15');T=O/'inputs/decontx_diagnostic'
print('INPUT_FILES',[(f.name,f.stat().st_size) for f in T.iterdir()]);p=O/'provenance/decontx_diagnostic.log';print(p.read_text()[-1500:] if p.exists() else 'input export not yet completed')
r=subprocess.run(['ps','-u','wangmingyuan','-o','pid,etime,pcpu,pmem,args'],capture_output=True,text=True);print('\n'.join(l for l in r.stdout.splitlines() if 'decontx' in l or ('python -' in l and 'grep' not in l)))