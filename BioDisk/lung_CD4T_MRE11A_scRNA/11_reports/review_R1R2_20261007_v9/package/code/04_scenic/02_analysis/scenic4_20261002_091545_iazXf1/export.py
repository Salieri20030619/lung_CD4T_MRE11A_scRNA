import os,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy import sparse
import anndata as ad
root=Path(os.environ['ROOT']);p1=Path(os.environ['PREV1']);out=Path(os.environ['OUT'])
source=root/'02_analysis/screen2_CD4_target_wz1lfx/GSE183682/CD4_core.h5ad'
a=ad.read_h5ad(source);assert a.n_obs==1413 and a.obs_names.is_unique and 'counts' in a.layers
x=sparse.csr_matrix(a.layers['counts'],dtype=float);x.sum_duplicates();x.eliminate_zeros()
assert np.isfinite(x.data).all() and (x.data>=0).all() and np.allclose(x.data,np.rint(x.data))
m=pd.read_csv(p1/'objects/GSE183682_cell_scores.tsv.gz',sep='\t');assert m.cell.is_unique and set(m.cell)==set(a.obs_names)
idx=pd.Index(a.obs_names).get_indexer(m.cell);x=x[idx];assert np.allclose(np.asarray(x.sum(axis=1)).ravel(),m.library)
assert m['sample'].nunique()==4 and set(m.day)=={7,56} and set(m.group)=={'NS','Silica'}
assert m[['sample','group','day','subtype','library','bg']].notna().all().all() and (m.library>0).all()
design={'GSM5567501':('NS',7),'GSM5567502':('Silica',7),'GSM5567503':('NS',56),'GSM5567504':('Silica',56)}
assert set(m['sample'])==set(design)
assert all((r.group,int(r.day))==design[r.sample] for r in m.itertuples())
assert np.array_equal(a.obs['sample'].astype(str).to_numpy()[idx],m['sample'].astype(str).to_numpy())
original=a.var['symbol'].astype('string').fillna('').astype(str).to_numpy() if 'symbol' in a.var else np.array(['']*a.n_vars)
ids=pd.Index(a.var_names.astype(str)).str.replace(r'\.\d+$','',regex=True)
alias={'MRE11A':'MRE11','ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','PARK2':'PRKN'}
genes=[]
for sy,gi in zip(original,ids):
 if not sy.strip() or sy.upper().startswith('ENS'):sy='' if gi.upper().startswith('ENS') else gi
 sy=sy.strip().upper();genes.append(alias.get(sy,sy))
genes=np.asarray(genes);valid=(genes!='')&~np.char.startswith(genes,'ENS');unique=sorted(set(genes[valid]))
pr=sparse.csr_matrix((np.ones(valid.sum()),(np.flatnonzero(valid),pd.Index(unique).get_indexer(genes[valid]))),shape=(len(genes),len(unique)))
x=(x@pr).tocsr();assert len(unique)>5000
sparse.save_npz(out/'input/validation_counts.npz',x)
pd.DataFrame({'gene':unique}).to_csv(out/'input/genes.tsv',sep='\t',index=False)
m.to_csv(out/'input/cells.tsv.gz',sep='\t',index=False)
pd.DataFrame({'feature':ids,'symbol':original,'canonical':genes,'mapped':valid}).to_csv(out/'tables/symbol_mapping.tsv.gz',sep='\t',index=False)
sm=m.groupby(['sample','group','day','subtype'],observed=True).size().reset_index(name='cells');sm.to_csv(out/'tables/input_design.tsv',sep='\t',index=False)
(out/'input_manifest.json').write_text(json.dumps(dict(source=str(source),metadata=str(p1/'objects/GSE183682_cell_scores.tsv.gz'),cells=len(m),genes=len(unique)),indent=2))
print('ANALYSIS_LOGIC_PASS EXTERNAL_COUNTS_FROZEN_METADATA',flush=True)

