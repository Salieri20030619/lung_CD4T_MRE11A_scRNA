from pathlib import Path
import os,json,hashlib,importlib.util
os.environ['OPENBLAS_NUM_THREADS']='1'
import pandas as pd,numpy as np,h5py
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/IKZF2_IL7R_target_deletion_20261005_v1';V=P/'11_reports/human_alignment_20261005_v1/results';G=P/'11_reports/human_gap_diagnostics_20261005_v1/results';N=P/'02_analysis/scenic1_20260929_103513_xRoEnr/Human/All';T=P/'11_reports/technical_delivery_20261002_ikzf2_stat1_v1/revision_v2/final_release'
for s in ['code','results','provenance','figures','deliverables']:(O/s).mkdir(parents=True,exist_ok=True)
m=pd.read_csv(N/'cells.tsv',sep='\t');print('CELL_METADATA',list(m.columns));counts=m.groupby(['sample','group','series']).size().rename('Core_cells').reset_index();ic=m[m.subtype=='CD4-IL7R'].groupby('sample').size();counts['IL7R_cells']=counts['sample'].map(ic).fillna(0).astype(int);counts['included16']=counts.IL7R_cells>=20;print('ROSTER',counts.to_string(index=False));counts.to_csv(O/'results/roster_17_to_16_audit.tsv',sep='\t',index=False)
a=pd.read_csv(V/'TF_readout_donors.tsv.gz',sep='\t');a=a[(a.readout=='AUC')&(a.scope=='CD4-IL7R')&(a.tf=='IKZF2')&(a['filter']=='All')&(a.variant=='DropOverlap')];print('AUC_ROSTER',a.groupby('seed')['sample'].agg(list).to_string());print('REF',pd.read_csv(G/'IL7R_fixed_reference_contrast_donors.tsv',sep='\t').head(2).to_string(index=False))
keep=sorted(counts.loc[counts.included16,'sample']);assert len(keep)==16
for seed in [17,29,43]:assert set(a[a.seed==seed]['sample'])==set(keep)
net=pd.read_csv(T/'data/TF_network_all_edges.tsv',sep='\t');net=net[(net.species=='Human')&(net['filter']=='All')&(net.tf=='IKZF2')&net.retained_DropOverlap];print('NET',net.groupby('seed').size().to_string());assert not(net.target=='IKZF2').any()
with h5py.File(N/'score.loom','r') as h:print('LOOM',h['matrix'].shape,'row_attrs',list(h['row_attrs']),'col_attrs',list(h['col_attrs']))
print('ALL_NETWORK',pd.read_csv(N/'regulon_edges.tsv',sep='\t').head(2).to_string(index=False))
for mod in ['pyscenic','ctxcore','loompy']:print('PACKAGE',mod,importlib.util.find_spec(mod))
if importlib.util.find_spec('pyscenic'):
 import inspect,pyscenic.aucell as au
 print('AUC_SIGNATURE',inspect.signature(au.aucell));print('AUC4R_SIGNATURE',inspect.signature(au.aucell4r));print('CREATE_RANKINGS',inspect.getsource(au.create_rankings));print('AUC4R_SOURCE',inspect.getsource(au.aucell4r));print('AUC_SOURCE',inspect.getsource(au.aucell))
print('PREFLIGHT_COMPLETE',flush=True)
