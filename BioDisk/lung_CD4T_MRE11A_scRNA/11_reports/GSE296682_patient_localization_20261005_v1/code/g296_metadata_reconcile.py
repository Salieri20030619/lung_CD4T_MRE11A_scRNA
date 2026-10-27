from pathlib import Path
import json,hashlib
import pandas as pd,numpy as np
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/GSE296682_patient_localization_20261005_v1';C=O/'public_cache';R=O/'results';f=C/'BAL_FINAL_META.csv';part=C/'BAL_FINAL_META.csv.part'
if not f.exists():
 assert part.stat().st_size==28044314,'Metadata still incomplete'
 assert hashlib.md5(part.read_bytes()).hexdigest()=='89a6e220d3ea1da5418027c6fd71f8b0';part.rename(f)
assert hashlib.md5(f.read_bytes()).hexdigest()=='89a6e220d3ea1da5418027c6fd71f8b0'
assert hashlib.sha256(f.read_bytes()).hexdigest()=='c34e6e8727d21c1c77816e280a237e391865bbb46d11b294c70f3a26f028f32b'
m=pd.read_csv(f).set_index('barcode');old=pd.read_csv(P/'00_raw/GSE296682/GSE296682_BAL_FINAL_metadata.txt.gz',sep='\t',index_col=0);q=m[m.disease_final=='Silicosis'];qo=old[old.disease_final=='Silicosis'];assert set(q.Sample)=={'B2_RZ725','B6_RZ725'};assert len(q)==4579
removed=old.loc[~old.index.isin(m.index)];added=m.loc[~m.index.isin(old.index)];d=dict(old_cells=len(old),new_cells=len(m),old_donors=old.Sample.nunique(),new_donors=m.Sample.nunique(),removed=len(removed),added=len(added),removed_by_sample=removed.Sample.value_counts().to_dict(),added_by_sample=added.Sample.value_counts().to_dict(),silicosis_cell_ids_equal=set(qo.index)==set(q.index));rows=[]
for col in q.columns.intersection(qo.columns):
 x=q[col];y=qo.loc[q.index,col]
 if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):eq=np.isclose(x,y,equal_nan=True,rtol=1e-6,atol=1e-6)
 else:eq=x.fillna('NA').astype(str).to_numpy()==y.fillna('NA').astype(str).to_numpy()
 rows.append(dict(column=col,total=len(x),equal=int(eq.sum()),different=int((~eq).sum())))
pd.DataFrame(rows).to_csv(R/'target_metadata_version_comparison.tsv',sep='\t',index=False);(R/'public_release_metadata_difference.json').write_text(json.dumps(d,indent=2));q.to_csv(R/'official_silicosis_metadata_SERVER_ONLY.tsv.gz',sep='\t');pd.crosstab(q.Cell_Subtype,q.Sample).to_csv(R/'official_target_subtype_counts.tsv',sep='\t');(O/'provenance/metadata_verified.json').write_text(json.dumps(dict(bytes=f.stat().st_size,md5=hashlib.md5(f.read_bytes()).hexdigest(),sha256=hashlib.sha256(f.read_bytes()).hexdigest(),url='https://zenodo.org/records/18626573/files/BAL_FINAL_META.csv?download=1'),indent=2));print(json.dumps(d,indent=2));print(pd.DataFrame(rows).to_string(index=False));print('META_VERIFIED')
