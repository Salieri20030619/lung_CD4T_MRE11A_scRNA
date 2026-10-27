from pathlib import Path
import sys,pandas as pd,numpy as np
from scipy.io import mmread
O=Path(sys.argv[1]);V=O.parent/'integrated_macrophage_CD4_20261015_v4';P=O.parents[1];M=P/'02_analysis/mechanism_step1_20260927_152218_s4V5eF';md=pd.read_csv(M/'comm_input/CD4_cells.tsv',sep='\t');gn=pd.read_csv(M/'comm_input/CD4_genes.tsv',sep='\t');x=mmread(M/'comm_input/CD4_counts.mtx.gz').T.tocsr();gi=pd.Index(gn.canonical);rows=[]
for sample in sorted(md['sample'].unique()):
 mask=md['sample'].eq(sample).to_numpy()
 for g in ['CD3D','CD3E','TRAC','CD247','CSF1R','LST1','TYROBP','FCER1G','LYZ2','EPCAM','SFTPA1','SFTPB','SFTPC','SFTPD','SCGB1A1','GSDMD']:
  j=gi.get_indexer([g])[0]
  if j<0:continue
  v=x[mask,j].toarray().ravel();rows.append(dict(sample=sample,gene=g,cells=int(mask.sum()),positive_cells=int((v>0).sum()),median_UMI=float(np.median(v)),max_UMI=float(v.max()),fraction_total_UMI=float(v.sum()/md.loc[mask,'library'].sum())))
pd.DataFrame(rows).to_csv(O/'tables/identity_marker_counts.tsv',sep='\t',index=False)
a=pd.read_csv(O/'tables/CD4_executor_identity_audit.tsv',sep='\t');print('IDENTITY_GSDMD');print(a[(a.scope=='Core')&(a.gene=='GSDMD')].to_string(index=False));print('SALINE2');print(pd.DataFrame(rows).query("sample=='Saline_2'").to_string(index=False));up=pd.read_csv(O/'tables/all_nominal_upstream_alternatives.tsv',sep='\t');print('ALTERNATIVES');print(up[up.x.isin(['C3','CXCL16','IL12A'])][['scope','sender','x','rho','P','BH','condition_P']].to_string(index=False));lr=pd.read_csv(V/'02_sender/ligand_scope_eligibility.tsv',sep='\t');lr[lr.ligand.isin(['C3','CXCL16','IL12A','IL27','IFNG'])].to_csv(O/'tables/alternative_receptor_coverage.tsv',sep='\t',index=False)