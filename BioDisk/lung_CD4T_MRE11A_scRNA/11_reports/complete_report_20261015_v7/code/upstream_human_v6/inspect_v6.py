from pathlib import Path
import sys,pandas as pd,numpy as np
O=Path(sys.argv[1]);P=O.parents[1];H=P/'11_reports/GSE174725_annotation_validation_20261005_v1/results';m=pd.read_csv(H/'frozen_annotation_cells_SERVER_ONLY.tsv.gz',sep='\t',index_col=0);a=pd.read_csv(P/'02_analysis/scenic2_20261001_092655_ZKNWD6/objects/Human_All_17_DropOverlap_auc.csv.gz',index_col=0);print('AUC',a.shape,a.index[:3].tolist(),'STAT',[c for c in a if c.startswith('STAT1')]);print('MD',m.shape,m.index[:3].tolist(),list(m.columns));print('OVERLAP',len(set(a.index)&set(m.index)),'core',int(m.old_CD4.ne('NotInCore').sum()))
J=P/'11_reports/GSE296682_patient_localization_20261005_v1/results'
for n in ['fixed_feature_mapping.tsv','official_silicosis_metadata_SERVER_ONLY.tsv.gz','all_fixed_and_marker_localization.tsv.gz']:
 d=pd.read_csv(J/n,sep='\t',nrows=4);print(n,list(d.columns));print(d.to_string(index=False)[:2000])
f=np.load(J/'official_fixed_expression_SERVER_ONLY.npz');print('NPZ',f.files);print([(k,f[k].shape) for k in f.files])