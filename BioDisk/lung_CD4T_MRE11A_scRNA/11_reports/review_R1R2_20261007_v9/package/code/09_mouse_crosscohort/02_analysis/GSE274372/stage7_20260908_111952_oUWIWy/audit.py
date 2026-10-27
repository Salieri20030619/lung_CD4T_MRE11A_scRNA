import os
from pathlib import Path
import numpy as np,pandas as pd,scanpy as sc
root=Path(os.environ['ROOT']);run=Path(os.environ['RUN']);out=run/'tables'
a=sc.read_h5ad(root/'02_analysis/GSE274372/stage5_20260907_224503_cUAf52/objects/T_refined.h5ad')
md=pd.read_csv(out/'author_metadata.tsv.gz',sep='\t',dtype=str,keep_default_na=False)
markers=['Cd4','Cd8a','Cd8b1','Zbtb7b','Runx3','Cd3d','Cd3e','Cd3g','Bcl11b','Themis','Skap1','Tcf7','Lef1','Ccr7','Sell','Il7r','Foxp3','Il2ra','Ctla4','Ikzf2','Icos','Nkg7','Ccl5','Prf1','Mre11a']
sym=a.var['symbol'].astype('string').fillna('').to_numpy()
e=a.obs.copy();e.index.name='cell_id'
for g in markers:
 idx=np.flatnonzero(sym==g)
 e[g]=np.asarray(a.layers['counts'][:,idx].sum(axis=1)).ravel() if len(idx) else np.nan
pd.DataFrame({'gene':markers,'feature_matches':[int((sym==g).sum()) for g in markers]}).to_csv(out/'marker_availability.tsv',sep='\t',index=False)
e.to_csv(out/'T_marker_evidence.tsv.gz',sep='\t')
e.loc[e['K20R08'].astype(str).eq('0')].to_csv(out/'Naive_marker_evidence.tsv.gz',sep='\t')
# Audit exact keys only. Repeated sci-RNA barcodes cannot establish sample identity.
rows=[]
for col in md.columns:
 vals=md[col].astype(str)
 for key,target in [('cell_id',pd.Series(a.obs_names)),('barcode',a.obs['barcode'].astype(str))]:
  matched=target.isin(set(vals))
  rows.append([col,key,int(matched.sum()),bool(vals.is_unique),bool(target.is_unique)])
pd.DataFrame(rows,columns=['author_column','local_key','matched_cells','author_unique','local_unique']).to_csv(out/'exact_key_audit.tsv',sep='\t',index=False)
ids=md['author_cell_id']; local=pd.Index(a.obs_names)
if ids.is_unique and local.isin(ids).any():
 joined=e.join(md.set_index('author_cell_id'),how='left',rsuffix='_author')
 joined.to_csv(out/'exact_id_join.tsv.gz',sep='\t')
(out/'notes.txt').write_text('Author labels are independent reference evidence, not ground truth. No fuzzy barcode matching or automatic relabeling. Barcode-only matches across samples are not sufficient. No expression imputation, new clustering, or time trend testing performed.\n')
print('GSE274372_STAGE7_ANALYSIS_PASS')
