import sys
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmwrite
src,out=map(Path,sys.argv[1:3])
a=sc.read_h5ad(src/'objects/CD4_core.h5ad')
assert a.n_obs==2311 and 'counts' in a.layers
p=out/'input';p.mkdir()
x=a.layers['counts'];assert np.all(x.data>=0) and np.allclose(x.data,np.round(x.data))
mmwrite(p/'counts.mtx',x.T)
pd.DataFrame({'gene_id':a.var_names,'symbol':a.var.symbol.values}).to_csv(p/'genes.tsv',sep='\t',index=False)
m=a.obs.copy();m.insert(0,'cell',a.obs_names)
m.to_csv(p/'cells.tsv',sep='\t',index=False)
print('EXPORT_PASS',flush=True)
