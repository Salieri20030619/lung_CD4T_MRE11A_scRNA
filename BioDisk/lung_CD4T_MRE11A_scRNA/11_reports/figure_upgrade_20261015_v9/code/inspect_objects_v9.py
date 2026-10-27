from pathlib import Path
import sys,pandas as pd,anndata as ad
O=Path(sys.argv[1]);P=O.parents[1]
for name in ['stage5_20260905_222457_a2dApT/objects/CD4_core.h5ad','stage11_20260906_100314_pALwAB/objects/macrophages_core.h5ad']:
 a=ad.read_h5ad(P/'02_analysis/GSE314943'/name,backed='r');print('OBJECT',name,a.shape);print('OBS',list(a.obs.columns),'OBSM',list(a.obsm.keys()));print(a.obs.head(2).to_string()[:1800]);a.file.close()
for name in ['02_analysis/scenic1_20260929_103513_xRoEnr/Mouse/All/regulon_edges.tsv','02_analysis/scenic1_20260929_103513_xRoEnr/Mouse/All/edge_stability.tsv','02_analysis/scenic2_20261001_092655_ZKNWD6/tables/candidate_edge_audit.tsv']:
 a=pd.read_csv(P/name,sep='\t');print('TABLE',name,a.shape,list(a.columns));print(a.head(2).to_string(index=False))