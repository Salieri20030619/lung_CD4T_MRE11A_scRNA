from pathlib import Path
import sys
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse,io
out=Path(sys.argv[1]);src=Path(sys.argv[2]);a=ad.read_h5ad(src)
assert a.n_obs==10286 and a.var_names.is_unique and a.obs_names.is_unique
assert 'MRE11' in a.var_names
cl=a.obs.F08.astype(str);a=a[cl.isin(['0','1','12','2','3'])].copy()
x=sparse.csr_matrix(a.layers['counts']);x.eliminate_zeros()
assert np.isfinite(x.data).all() and (x.data>=0).all() and np.allclose(x.data,np.round(x.data))
m=a.obs.copy();m['cell']=a.obs_names.astype(str);m['F08']=m.F08.astype(str)
m['subtype']=m.F08.map({'0':'MemCD4','1':'NaiveCD4','12':'Treg','2':'Mix2','3':'Mix3'})
m['core']=m.F08.isin(['0','1','12']);m['MRE11_UMI']=x[:,list(a.var_names).index('MRE11')].toarray().ravel();m['library_UMI']=np.asarray(x.sum(axis=1)).ravel()
m['group']=m['sample'].map({'GSM5324866':'Exposed','GSM5324867':'Exposed','GSM5324868':'Silicosis','GSM5324869':'Silicosis','GSM5324870':'Silicosis'})
assert m.group.notna().all() and (m.library_UMI>0).all()
old=pd.read_csv(out/'input/expected_samples.tsv',sep='\t');old=old[old.dataset=='GSE174725']
for scope in ['Core','Expanded']:
 z=m[m.core] if scope=='Core' else m
 for ss,v in z.groupby('sample'):
  e=old[(old.scope==scope)&(old['sample']==ss)].iloc[0]
  assert len(v)==e.n_cells and v.MRE11_UMI.sum()==e.target_counts and np.isclose(v.library_UMI.sum(),e.total_counts)
assert m.core.sum()==2958 and len(m)==5481
m.to_csv(out/'input/cells.tsv.gz',sep='\t',index=False)
pd.DataFrame({'symbol':a.var_names}).to_csv(out/'input/genes.tsv',sep='\t',index=False)
io.mmwrite(out/'input/counts.mtx',x.T)
(out/'input/provenance.txt').write_text('Object='+str(src)+'\nCounts=layers[counts]; full measured feature background; pre-fixed F08 Core=0,1,12; Expanded adds 2,3.\n')
print('ANALYSIS_LOGIC_PASS INPUT_COUNTS_AND_SCOPE',flush=True)
