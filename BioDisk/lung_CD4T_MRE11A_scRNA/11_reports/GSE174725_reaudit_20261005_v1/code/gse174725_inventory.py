from pathlib import Path
import subprocess,shutil,json
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
print('HOST',subprocess.run(['hostname'],capture_output=True,text=True).stdout.strip())
p=subprocess.run(['ps','-eo','pid,etime,args'],capture_output=True,text=True)
print('MATCHING_PROCESSES','\n'.join(x for x in p.stdout.splitlines() if ('174725' in x or 'human_gap' in x) and 'base64' not in x and 'ps -eo' not in x))
if shutil.which('rg'):
 r=subprocess.run(['rg','--files','--hidden',str(P)],capture_output=True,text=True)
 files=r.stdout.splitlines()
else:files=[str(p) for p in P.rglob('*') if p.is_file()]
hits=[p for p in files if '174725' in p or '532486' in p]
for p in hits:
 q=Path(p);print('FILE',q.stat().st_size,p)
print('TEXT_SEARCH')
if shutil.which('rg'):
 r=subprocess.run(['rg','-l','-i','--glob','*.py','--glob','*.R','--glob','*.log','--glob','*.md','--glob','*.txt','--glob','*.json','174725',str(P/'11_reports'),str(P/'02_analysis')],capture_output=True,text=True);print(r.stdout[:30000])
