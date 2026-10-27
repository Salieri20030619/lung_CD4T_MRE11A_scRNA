import sys
from pathlib import Path
import numpy as np,pandas as pd,scanpy as sc
from scipy.io import mmwrite
from scipy.sparse import vstack,csr_matrix
root=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/02_analysis/GSE314943')
out=Path(sys.argv[1]);(out/'input').mkdir()
ma=sc.read_h5ad(root/'stage11_20260906_100314_pALwAB/objects/macrophages_core.h5ad')
t=sc.read_h5ad(root/'stage5_20260905_222457_a2dApT/objects/CD4_core.h5ad')
assert ma.n_obs==1796 and t.n_obs==2311
assert ma.var_names.equals(t.var_names) and not ma.obs_names.isin(t.obs_names).any()
assert ma.var.symbol.equals(t.var.symbol)
frames=[]
for a,lineage,field in [(ma,'Mac','MacType'),(t,'CD4','Subtype')]:
 frames.append(pd.DataFrame({'cell':a.obs_names,'sample':a.obs['sample'].values,'group':a.obs['group'].values,'lineage':lineage,'subtype':a.obs[field].astype(str).values}))
m=pd.concat(frames,ignore_index=True);assert m.cell.is_unique
x=vstack([ma.layers['counts'],t.layers['counts']],format='csr')
# Sum raw counts for duplicate gene symbols, then normalize in R.
symbols=ma.var['symbol'].astype('string').fillna('').to_numpy(dtype=str)
valid=np.array([bool(s.strip()) for s in symbols]);unique=pd.Index(sorted(set(symbols[valid])))
idx=unique.get_indexer(symbols[valid]);proj=csr_matrix((np.ones(valid.sum()),(np.flatnonzero(valid),idx)),shape=(len(symbols),len(unique)))
x=x@proj
mmwrite(out/'input/counts.mtx',x.T)
pd.DataFrame({'symbol':unique}).to_csv(out/'input/genes.tsv',sep='\t',index=False)
m.to_csv(out/'input/cells.tsv',sep='\t',index=False)
pd.DataFrame({'gene_id':ma.var_names,'symbol':symbols,'retained':valid}).to_csv(out/'input/gene_mapping.tsv',sep='\t',index=False)
print('EXPORT_PASS',flush=True)
