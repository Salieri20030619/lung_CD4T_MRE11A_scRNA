from pathlib import Path
import re
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');T=P/'11_reports/technical_delivery_20261002_ikzf2_stat1_v1/revision_v2/final_release'
for f in (T/'provenance').glob('*'):
 if f.is_file() and f.suffix.lower() in ['.txt','.md']:
  ss=f.read_text(errors='replace')
  if re.search(r'BH|FDR|correlation',ss,re.I):print(f.name,ss[:15000])
for f in (P/'02_analysis').glob('scenic2*/**/*.py'):
 ss=f.read_text(errors='replace');ls=ss.splitlines();ix=set()
 for i,l in enumerate(ls):
  if re.search(r'groupby|FDR|bh|multipletest',l,re.I):ix.update(range(max(0,i-2),min(len(ls),i+4)))
 print(str(f), '\n'.join(f'{i+1}: {ls[i]}' for i in sorted(ix))[-18000:])
f=P/'logs/17_GSE288887_CD4_pathway_sensitivity.log';print('LOGEND',f.read_text(errors='replace')[-6000:])
