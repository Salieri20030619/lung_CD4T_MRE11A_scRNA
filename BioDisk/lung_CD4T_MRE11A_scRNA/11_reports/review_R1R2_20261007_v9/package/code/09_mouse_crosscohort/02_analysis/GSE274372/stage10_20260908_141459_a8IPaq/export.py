import os
from pathlib import Path
import scanpy as sc,numpy as np,pandas as pd
root=Path(os.environ['ROOT']);run=Path(os.environ['RUN']);out=run/'tables';out.mkdir()
a=sc.read_h5ad(root/'02_analysis/GSE274372/stage9_20260908_140721_Yj56oq/objects/CD4_author_reclustered.h5ad')
assert a.n_obs==888
x=a.layers['counts'];assert np.all(x.data>=0) and np.allclose(x.data,np.round(x.data))
sy=a.var['symbol'].astype('string').fillna('');ix=np.flatnonzero(sy.eq('Mre11a'));assert len(ix)==1
cl=a.obs.R08.astype(str)
scopes={'TLike':cl.isin(['0','1']),'NoMix':~cl.isin(['4','5','6']),'Author':np.ones(a.n_obs,dtype=bool)}
rows=[]
for name,selection in scopes.items():
 cols={}
 for sample in sorted(a.obs['sample'].astype(str).unique()):
  mask=np.asarray(selection)&a.obs['sample'].astype(str).eq(sample).to_numpy()
  v=np.asarray(x[mask].sum(axis=0)).ravel().astype('int64');cols[sample]=v
  target=np.asarray(x[mask][:,ix].sum(axis=1)).ravel();n=int(mask.sum())
  rows.append([name,sample,int(a.obs.loc[a.obs['sample'].astype(str).eq(sample),'day'].iloc[0]),n,int((target>0).sum()),int(target.sum()),int(v.sum())])
 pd.DataFrame(cols,index=a.var_names).to_csv(out/f'{name}_counts.tsv.gz',sep='\t')
pd.DataFrame(rows,columns=['scope','sample','day','cells','Mre_positive','Mre_UMI','library_UMI']).to_csv(out/'sample_summary.tsv',sep='\t',index=False)
a.var[['symbol']].to_csv(out/'features.tsv.gz',sep='\t')
a.obs.to_csv(out/'cell_metadata.tsv.gz',sep='\t')
(out/'scope_notes.txt').write_text('TLike: R08 clusters0/1, activated-regulatory and naive T-like author CD4 candidates, not definitively pure CD4. NoMix: excludes clusters4/5/6; retains uncertain clusters2/3. Author: all888. No identities overwritten. Independent samples are analysis units. Minimum10 cells/sample is an operational inclusion rule, not a universal validity threshold. Per-day biological replicates remain small. No Mre11a-positive gating.\n')
print('EXPORT_PASS',flush=True)
