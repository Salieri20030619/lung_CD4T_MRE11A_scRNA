from pathlib import Path
import pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');B=P/'11_reports';S=B/'independent_candidates_revision_20261005_v1/results/source_evidence'
d=pd.read_csv(next(S.glob('old_TF_program.tsv*')),sep='\t')
print(d[(d.species=='Human')&(d.tf=='STAT1')&(d['filter']=='All')&(d.variant=='DropOverlap')&(d.scope.isin(['GZMA','Core']))&(d.model=='SourceDep')&(d.feature=='score_Inflam')].to_string(index=False))
T=B/'technical_delivery_20261002_ikzf2_stat1_v1/revision_v2/final_release'
for root in [T,P/'10_tables/CD4_final/pathway_v1']:
 print('ROOT',root)
 for x in root.rglob('*'):
  if x.suffix.lower() in ['.py','.r','.md','.json','.txt']:
   print('CODE',x)
   txt=x.read_text(errors='replace')
   lines=txt.splitlines()
   for i,l in enumerate(lines):
    if any(w in l for w in ['SourceDep','S2_TF_associations','source_adjusted_full_DEG','score_Inflam','multipletests','p.adjust']):
     print('\n'.join(f'{j+1}: {lines[j]}' for j in range(max(0,i-2),min(len(lines),i+4))))
