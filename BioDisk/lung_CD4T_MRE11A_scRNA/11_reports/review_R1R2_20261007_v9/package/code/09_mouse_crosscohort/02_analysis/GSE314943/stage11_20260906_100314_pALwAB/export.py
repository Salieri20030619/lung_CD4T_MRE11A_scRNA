import sys
from pathlib import Path
import pandas as pd,numpy as np,scanpy as sc
from scipy.io import mmwrite
root=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/02_analysis/GSE314943')
out=Path(sys.argv[1]);(out/'input').mkdir();(out/'objects').mkdir();(out/'tables').mkdir()
a=sc.read_h5ad(root/'stage9_20260906_080538_wN3sjO/objects/K15_trajectory.h5ad')
pt=pd.read_csv(root/'stage9_20260906_080538_wN3sjO/tables/pseudotime_all_runs.tsv.gz',sep='\t',index_col=0)
assert a.n_obs==1734 and a.obs_names.isin(pt.index).all()
m=a.obs.copy();m.insert(0,'cell',a.obs_names)
m['PT']=pt.loc[a.obs_names,'K15_R1'].values;m['AltPT']=pt.loc[a.obs_names,'K15_R2'].values
m.to_csv(out/'input/cells.tsv',sep='\t',index=False)
mmwrite(out/'input/counts.mtx',a.layers['counts'].T)
pd.DataFrame({'gene_id':a.var_names,'symbol':a.var.symbol.values}).to_csv(out/'input/genes.tsv',sep='\t',index=False)
mac=sc.read_h5ad(root/'stage10_20260906_094533_ZUUfJn/objects/myeloid_candidates.h5ad')
mp={'1':'AMlike','2':'Spp1Mac','5':'C1qMac','11':'Folr2Mac'}
mac.obs['CoreMac']=mac.obs.M06.astype(str).isin(mp)
mac.obs.to_csv(out/'tables/myeloid_membership.tsv.gz',sep='\t')
mac=mac[mac.obs.CoreMac].copy();mac.obs['MacType']=mac.obs.M06.astype(str).map(mp).astype('category')
assert mac.n_obs==1796
pd.crosstab(mac.obs.MacType,mac.obs['sample']).to_csv(out/'tables/macrophage_sample.tsv',sep='\t')
mac.write_h5ad(out/'objects/macrophages_core.h5ad',compression='gzip')
print('EXPORT_PASS',flush=True)
