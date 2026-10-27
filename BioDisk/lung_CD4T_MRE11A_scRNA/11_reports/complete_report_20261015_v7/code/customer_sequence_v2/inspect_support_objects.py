from pathlib import Path
import pandas as pd,anndata as ad,sys
O=Path(sys.argv[1]);P=O.parents[1]
r=pd.read_csv(O/'03A_state_pseudotime/genomewide_trajectory_tests.tsv.gz',sep='\t');print(r[r.feature.isin(['AHNAK','MDFIC','TTC39C','ST6GALNAC3'])][['state','run','feature','mean_slope','P','BH','same_sign_fraction']].to_string(index=False))
H=P/'11_reports/GSE174725_annotation_validation_20261005_v1/results'
for name in ['frozen_scope_donor_counts.tsv','frozen_cluster_annotation_review.tsv']:
 print(name);print((H/name).read_text(encoding='utf-8')[:4000])
a=ad.read_h5ad(H/'broad_annotation_review_SERVER_ONLY.h5ad',backed='r');print('HUMAN OBJECT',a.shape,'OBS',list(a.obs.columns),'VAR',list(a.var.columns),'LAYERS',list(a.layers.keys()),'RAW',None if a.raw is None else a.raw.shape);print(a.obs.head(2).to_string(index=False));print(a.var.head(3).to_string());a.file.close()
for n in ['all_fixed_and_marker_localization.tsv.gz','fixed_feature_mapping.tsv']:
 a=pd.read_csv(P/'11_reports/GSE296682_patient_localization_20261005_v1/results'/n,sep='\t');print('LOCALIZE',n,a.shape,list(a.columns));print(a.head(2).to_string(index=False))