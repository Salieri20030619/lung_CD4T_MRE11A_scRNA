import sys
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse,io
out=Path(sys.argv[1]); src=Path(sys.argv[2]); a=ad.read_h5ad(src)
x=sparse.csr_matrix(a.layers['counts'])
assert a.n_obs==888 and x.shape==a.shape
assert np.isfinite(x.data).all() and (x.data>=0).all() and np.allclose(x.data,np.round(x.data))
ids=a.var_names.astype(str).str.replace(r'\.\d+$','',regex=True)
assert ids.is_unique and 'ENSMUSG00000031928' in ids
m=a.obs.copy();m['cell']=a.obs_names.astype(str);m['R08']=m.R08.astype(str);m['day']=pd.to_numeric(m.day)
assert m.cell.is_unique and m['sample'].nunique()==12 and set(m.day)=={0,7,28,56}
assert m.R08.isin(['0','1']).sum()==494
m['Mre_UMI']=x[:,list(ids).index('ENSMUSG00000031928')].toarray().ravel()
m['library_UMI']=np.asarray(x.sum(axis=1)).ravel()
assert (m.library_UMI>0).all()
m.to_csv(out/'input/cells.tsv',sep='\t',index=False)
pd.DataFrame({'gene_id':ids}).to_csv(out/'input/genes.tsv',sep='\t',index=False)
io.mmwrite(out/'input/counts.mtx',x.T)
print('EXPORT_PASS',a.shape,flush=True)
