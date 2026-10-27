from pathlib import Path
import pandas as pd
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/three_cohort_four_axis_20261017_v16');d=pd.read_csv(O/'tables/DE_all_models.tsv',sep='\t');f=pd.read_csv(O/'tables/function_membership.tsv',sep='\t');a=d[d.model.isin(['314_disease','joint_stratified','274_time']) & (d.PValue<.05) & d.gene.isin(f.gene)]
print('DE_FUNCTION',a.gene.nunique(),';'.join(sorted(a.gene.unique())))
for x in ['S17_逐基因说明.tsv','S17_候选基准.tsv']:
 p=O.parent/'methods_complete_20261016_v13/supplement'/x;z=pd.read_csv(p,sep='\t');print(x,z.columns.tolist(),z.head(2).to_string(index=False))