from pathlib import Path
import json,subprocess
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
print('PROCESSES','\n'.join(x for x in subprocess.run(['ps','-eo','pid,etime,args'],capture_output=True,text=True).stdout.splitlines() if ('g174_' in x or 'gse174725_' in x) and 'base64' not in x))
for folder in ['02_analysis','11_reports']:
 for p in (P/folder).rglob('*'):
  if p.is_file() and p.suffix in ['.py','.R','.md','.txt','.json'] and p.stat().st_size<1000000:
   t=p.read_text(errors='replace')
   if ('GSE174725' in t or 'GSM5324866' in t) and any(w in t.lower() for w in ['macroph','myeloid','annotation','annotat','doublet','ambient']):
    print('HIT',p)
    if p.name in ['annotation.json','README.txt','methods.txt']:print(t[:3500])
for n in ['02_analysis/screen2_TNK_20260911_103239_nfW0yu/analysis.py','02_analysis/screen2_Trefine_nuzVQu/refine.py','02_analysis/screen2_CD4_target_wz1lfx/README.txt']:
 p=P/n
 if p.exists():print('CODE',n,p.read_text()[:12000])
