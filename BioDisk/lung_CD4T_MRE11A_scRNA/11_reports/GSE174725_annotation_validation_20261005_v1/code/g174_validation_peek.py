from pathlib import Path
import pandas as pd
R=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/GSE174725_annotation_validation_20261005_v1/results')
p=R/'mature_nonAM_macrophage_marker_audit.tsv'
if p.exists():
 d=pd.read_csv(p,sep='\t');print(d[d.cluster.isin([0,1,6,9,17,18,19])].pivot(index='gene',columns='cluster',values='pct').round(3).to_string())
 p=pd.read_csv(R/'all_fixed_validation_tests.tsv',sep='\t');q=p[(p.cell_filter=='All')&(p.min_cells==1)&p.family.isin(['receiver_RNA','receiver_program'])];print(q[['scope','endpoint','n_exposed','n_silicosis','effect','P_exact','FDR_BH']].to_string(index=False));print('STATUS',p.groupby(['cell_filter','min_cells','family','status']).size().to_string())
else:print('NOT_COMPLETE')
