from pathlib import Path
import h5py,numpy as np,pandas as pd,json
from anndata.io import read_elem
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/cohort_bridge_20261017_v25';O.mkdir(exist_ok=True);(O/'tables').mkdir(exist_ok=True)
f=P/'02_analysis/screen2_CD4_target_wz1lfx/GSE183682/CD4_core.h5ad'
with h5py.File(f,'r') as h:
 m=read_elem(h['obs']);keys=list(h['obsm'].keys());print('SHAPE',len(m),'OBS',list(m.columns),'EMBEDDINGS',keys,flush=True)
 emb='X_umap_Harmony' if 'X_umap_Harmony' in keys else 'X_umap' if 'X_umap' in keys else None
 if emb is None:raise RuntimeError('No stored UMAP; do not invent combined embedding')
 xy=h['obsm'][emb][:];z=pd.DataFrame({'cell':m.index,'sample':m['sample'].astype(str).to_numpy(),'state':m['core_label'].astype(str).replace({'NaiveCD4':'Naive'}).to_numpy(),'UMAP1':xy[:,0],'UMAP2':xy[:,1]})
 mp={'GSM5567501':('Saline',7),'GSM5567502':('Silica',7),'GSM5567503':('Saline',56),'GSM5567504':('Silica',56)};z['group']=z['sample'].map({k:v[0] for k,v in mp.items()});z['day']=z['sample'].map({k:v[1] for k,v in mp.items()});assert not z[['group','day']].isna().any().any()
 z.to_csv(O/'tables/GSE183682_CD4_coordinates.tsv.gz',sep='\t',index=False);meta={'source':str(f),'embedding':emb,'cells':len(z),'samples':z['sample'].nunique(),'states':z.state.value_counts().to_dict(),'note':'existing cohort-specific embedding; no integration/refit'};(O/'supplement_embedding_source.json').write_text(json.dumps(meta,indent=2),encoding='utf-8');print(json.dumps(meta),flush=True)