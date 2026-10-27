from pathlib import Path
import os,json
os.environ['OPENBLAS_NUM_THREADS']='1'
import h5py,numpy as np,pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/GSE296682_patient_localization_20261005_v1';f=P/'00_raw/GSE296682/GSE296682_BAL_FINAL.h5ad'
with h5py.File(f,'r') as h:
 x=h['X'];v=x['data'][:100000];info=dict(source='ARCHIVED_GEO_NOT_CURRENT_ZENODO',X_encoding=dict(x.attrs),value_sample_min=float(v.min()),value_sample_max=float(v.max()),value_sample_integer=float(np.mean(np.isclose(v,np.round(v)))),layers=list(h['layers'].keys()),raw='raw' in h)
 def cv(v):
  if isinstance(v,np.ndarray):return v.tolist()
  if isinstance(v,bytes):return v.decode()
  if isinstance(v,np.generic):return v.item()
  return str(v)
 print(json.dumps(info,default=cv));(O/'provenance/archived_GEO_layer_schema.json').write_text(json.dumps(info,default=cv,indent=2))
m=pd.read_csv(P/'00_raw/GSE296682/GSE296682_BAL_FINAL_metadata.txt.gz',sep='\t',index_col=0);q=m[m.disease_final=='Silicosis'];print(pd.crosstab(q.Cell_Subtype,q.Sample).to_string())
