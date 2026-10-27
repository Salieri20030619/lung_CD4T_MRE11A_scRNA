from pathlib import Path
import pandas as pd,json,hashlib
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');Q=P/'11_reports/review_revision_20261017_v14/supplement_v14'
a=pd.read_csv(Q/'S24_complex_codetection_raw.tsv',sep='\t');print('CASE COMPLEX',a[(a.group=='Silica')&a.scope.isin(['ActCD4','AMlike','AllMac'])][['scope','sample','complex','positive_cells','cells','pct']].to_string(index=False))
a=pd.read_csv(Q/'S26_whole_lung_competing_sources.tsv',sep='\t');print('SOURCE',a.groupby(['gene','group','source']).agg(mean_CP10k=('mean_CP10k','mean'),mean_pct=('positive_pct','mean'),UMI=('total_UMI','sum'),samples=('sample','nunique')).reset_index().sort_values(['gene','group','mean_CP10k'],ascending=[True,True,False]).groupby(['gene','group']).head(7).to_string(index=False))
a=pd.read_csv(Q/'S21_complete_STAT1_fixed.tsv',sep='\t');print('FIXED',a[a.y=='GSDMD'].to_string(index=False))
a=pd.read_csv(Q/'S23_targeted_proximal_associations.tsv',sep='\t');print('HITS',a[a.P<.05].to_string(index=False))