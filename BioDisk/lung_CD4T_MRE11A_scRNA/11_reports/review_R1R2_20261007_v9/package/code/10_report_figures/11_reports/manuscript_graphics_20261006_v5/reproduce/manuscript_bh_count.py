from pathlib import Path
import pandas as pd,json
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');R=next((P/'02_analysis').glob('scenic2_*'));O=P/'11_reports/manuscript_integrated_20261006_v3/provenance'
ls=(R/'analyze.py').read_text().splitlines();print('\n'.join(ls[96:122]))
for f in (R/'tables').glob('*'):
 if 'association' in f.name:
  print('FILE',f.name)
  d=pd.read_csv(f,sep='\t');print(list(d))
  if set(['species','filter','seed','variant','scope','model','p','FDR']).issubset(d):
   z=d[(d.species=='Human')&(d.scope=='Core')&(d.model=='SourceDep')];a=z.groupby(['filter','seed','variant']).agg(rows=('p','size'),finite_tests=('p','count'),TFs=('regulon','nunique')).reset_index();print(a.to_string(index=False));a.to_json(O/'TF_program_BH_sizes.json',orient='records',indent=2)
print('GSEA source filename done')
