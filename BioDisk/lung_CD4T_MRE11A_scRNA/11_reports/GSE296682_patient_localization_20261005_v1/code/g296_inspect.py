from pathlib import Path
import os,json
os.environ['OPENBLAS_NUM_THREADS']='1'
import pandas as pd,numpy as np,h5py,anndata as ad
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/GSE296682_patient_localization_20261005_v1';R=O/'results';C=O/'public_cache'
assert (O/'provenance/official_download_verified.json').exists()
m=pd.read_csv(C/'BAL_FINAL_META.csv').set_index('barcode');print('META',m.shape,m.columns.tolist(),flush=True);q=m[m.disease_final=='Silicosis'];assert set(q.Sample)=={'B2_RZ725','B6_RZ725'};print('SIL_COUNTS',pd.crosstab(q.Cell_Subtype,q.Sample).to_string(),flush=True);q.to_csv(R/'official_silicosis_metadata_SERVER_ONLY.tsv.gz',sep='\t')
f=C/'BAL_FINAL.h5ad';structure=[]
with h5py.File(f,'r') as h:
 def visitor(name,o):
  if name.count('/')<3:structure.append(dict(path=name,kind=type(o).__name__,shape=str(getattr(o,'shape','')),dtype=str(getattr(o,'dtype','')),attrs={k:str(v) for k,v in o.attrs.items()}))
 h.visititems(visitor)
(O/'provenance/h5ad_structure.json').write_text(json.dumps(structure,indent=2));print('STRUCTURE',json.dumps([x for x in structure if x['path'].split('/')[0] in ['X','layers','raw','var','uns']],indent=2)[:18000],flush=True)
a=ad.read_h5ad(f,backed='r');assert a.obs_names.is_unique and a.var_names.is_unique;print('OBJECT',a.shape,'layers',list(a.layers.keys()),'raw',None if a.raw is None else a.raw.shape,'var',a.var.head().to_string(),'obs',a.obs.columns.tolist(),'uns',str(a.uns)[:3000],flush=True)
common=q.index.intersection(a.obs_names);print('MATCHED_SIL',len(common),'of',len(q),'objecttotal',a.n_obs,'metadata_total',len(m),flush=True);assert len(common)==len(q);b=a[common].to_memory();a.file.close();b.obs=m.loc[b.obs_names].copy();b.write_h5ad(R/'two_silicosis_official_SERVER_ONLY.h5ad',compression='gzip')
rows=[]
for key,v in [('X',b.X)]+[(k,b.layers[k]) for k in b.layers.keys()]+([('raw.X',b.raw.X)] if b.raw is not None else []):
 data=v.data if hasattr(v,'tocsr') else np.asarray(v).ravel();tot=np.asarray(v.sum(axis=1)).ravel();stored=np.expm1(pd.to_numeric(b.obs.log1p_total_counts).to_numpy());row=dict(layer=key,shape=str(v.shape),dtype=str(v.dtype),nonzero=int(np.count_nonzero(data)),min=float(data.min()),max=float(data.max()),fraction_integer=float(np.mean(np.isclose(data,np.round(data),atol=1e-7))),rowtotal_min=float(tot.min()),rowtotal_median=float(np.median(tot)),rowtotal_max=float(tot.max()),relative_error_to_metadata_total_median=float(np.median(np.abs(tot-stored)/stored)));rows.append(row)
pd.DataFrame(rows).to_csv(R/'matrix_layer_value_audit.tsv',sep='\t',index=False);print('LAYERS',json.dumps(rows,indent=2),flush=True)
# Public release difference: metadata identity only, do not modify either version.
old=pd.read_csv(P/'00_raw/GSE296682/GSE296682_BAL_FINAL_metadata.txt.gz',sep='\t',index_col=0);removed=old.loc[~old.index.isin(m.index)];added=m.loc[~m.index.isin(old.index)];diff=dict(old_cells=len(old),new_cells=len(m),old_donors=old.Sample.nunique(),new_donors=m.Sample.nunique(),removed=len(removed),added=len(added),removed_by_sample=removed.Sample.value_counts().to_dict(),added_by_sample=added.Sample.value_counts().to_dict(),silicosis_cell_ids_equal=set(old[old.disease_final=='Silicosis'].index)==set(q.index));(R/'public_release_metadata_difference.json').write_text(json.dumps(diff,indent=2));print('RELEASE_DIFFERENCE',json.dumps(diff),flush=True)
