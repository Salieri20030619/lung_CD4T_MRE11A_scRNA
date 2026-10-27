from pathlib import Path
import pandas as pd
R=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/GSE296682_patient_localization_20261005_v1/results')
d=pd.read_csv(R/'all_fixed_and_marker_localization.tsv.gz',sep='\t');c=pd.read_csv(R/'same_cell_codetection_and_marker_context.tsv',sep='\t');p=pd.read_csv(R/'fixed_gene_positive_marker_context.tsv.gz',sep='\t');sc=['Alveolar macrophage','Monocyte-like macrophage','Interstitial-like macrophage','Proliferating alveolar macrophage','CD4 effector T cell','CD4 helper T cell','Regulatory T cell']
print('CO',c[c.scope.isin(sc)&c.pair.isin(['IL15_IL15RA','IL2RB_IL2RG','TCR2','CD4_detected','CD8AB','NKcyto2','MatureMac3','Mono2','AM2'])][['sample','scope','pair','N','n_joint','pct']].to_string(index=False))
markers='C1QA C1QC APOE APOC1 MRC1 TREM2 GPNMB LGMN SPP1 FABP4 PPARG FCN1 VCAN S100A8 CD3D TRAC CD4 CD8A CD8B NKG7 GNLY FOXP3 IL2RA'.split()
for sample in ['B2_RZ725','B6_RZ725']:
 q=d[(d['sample']==sample)&d.scope.isin(sc)&d.gene.isin(markers)];print('MARKERS',sample,q.pivot(index='gene',columns='scope',values='pct').round(3).to_string())
print('FIXED12',d[d.scope.isin(['Alveolar macrophage','Monocyte-like macrophage','CD4 effector T cell'])&d.gene.isin('IL15 TNF CLSTN1 ADAM15 PLAU DSC2 CXCL2 ICAM2 VTN VEGFB GSTP1 RTN4'.split())][['sample','scope','gene','N','n_positive','pct']].to_string(index=False))
print('POSCONTEXT',p[p.scope.isin(['Alveolar macrophage','Monocyte-like macrophage','CD4 effector T cell'])&p.gene.isin(['IL15','IKZF2','IL2RB','IL2RG'])&p.context.isin(['TCR2','CD8AB','NKcyto2','MatureMac3'])].to_string(index=False))
old=pd.read_csv(R/'prior_main_fixed_expression_REUSED_SERVER_ONLY.tsv.gz',sep='\t');print('OLD_KEYS',old.key.unique());print('OLD_GROUPS',old.group.unique());print(old.head().to_string(index=False))
