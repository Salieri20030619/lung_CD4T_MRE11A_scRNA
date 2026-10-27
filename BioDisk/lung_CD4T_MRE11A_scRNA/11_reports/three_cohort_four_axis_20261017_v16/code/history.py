from pathlib import Path
import pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');s=P/'11_reports/methods_complete_20261016_v13/supplement';o=P/'11_reports/three_cohort_four_axis_20261017_v16'
records=[]
for f in s.glob('S17*.tsv'):
 d=pd.read_csv(f,sep='\t');print(repr(f.name),d.shape,d.columns.tolist());print(d.head(1).to_string(index=False))
 for col in ['gene','canonical','symbol','候选基因','基因']:
  if col in d:
   for v in d[col].dropna().astype(str).unique():records.append({'gene':v.upper(),'file':str(f),'field':col})
pd.DataFrame(records).to_csv(o/'inputs/historical_candidates_extracted.tsv',sep='\t',index=False)
print('LOGS',[f.name for f in (o/'logs').iterdir()])