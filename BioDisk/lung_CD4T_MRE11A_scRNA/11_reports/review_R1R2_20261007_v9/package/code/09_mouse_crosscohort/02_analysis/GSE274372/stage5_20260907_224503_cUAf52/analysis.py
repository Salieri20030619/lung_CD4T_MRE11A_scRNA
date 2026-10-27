import os
from pathlib import Path
import numpy as np,pandas as pd,scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import adjusted_rand_score
root=Path(os.environ['ROOT']);run=Path(os.environ['RUN'])
for p in ['tables','figures','objects']:(run/p).mkdir()
a=sc.read_h5ad(root/'02_analysis/GSE274372/stage4_20260907_213653_9UOCET/objects/T_candidates.h5ad')
a.obs['stage4_T05']=a.obs.T05.astype(str)
a=a[a.obs.stage4_T05.isin(['1','2'])&~a.obs.dbl_020].copy()
assert a.n_obs==672
# Clear prior neighborhood and embedding, preserve original counts and metadata.
a.X=a.layers['counts'].copy();a.uns.clear()
for p in list(a.obsp):del a.obsp[p]
for p in list(a.obsm):del a.obsm[p]
sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a)
sc.pp.highly_variable_genes(a,n_top_genes=1500,flavor='seurat',batch_key='sample')
v=a[:,a.var.highly_variable].copy();sc.pp.scale(v,max_value=10)
sc.tl.pca(v,n_comps=30,svd_solver='arpack',random_state=42)
a.obsm['X_pca']=v.obsm['X_pca'].copy()
a.var[['symbol','highly_variable']].to_csv(run/'tables/HVG.tsv',sep='\t')
panel=['Cd3d','Cd3e','Cd3g','Bcl11b','Themis','Skap1','Cd4','Cd8a','Cd8b1','Zbtb7b','Runx3','Tcf7','Lef1','Ccr7','Sell','Il7r','Foxp3','Il2ra','Ctla4','Ikzf2','Icos','Pdcd1','Gata3','Rorc','Tbx21','Nkg7','Klrd1','Prf1','Gzmb','Mki67','Sftpc','Scgb1a1','Lyz2','Mre11a']
rows=[];top=[];keys=[]
for k in [10,20]:
 sc.pp.neighbors(a,n_neighbors=k,n_pcs=20,random_state=42)
 sc.tl.umap(a,random_state=42)
 for res in [0.4,0.8,1.2]:
  key=f'K{k}R{int(res*10):02d}';keys.append(key)
  sc.tl.leiden(a,resolution=res,key_added=key,random_state=42)
  sc.tl.rank_genes_groups(a,key,method='wilcoxon',use_raw=False,n_genes=50)
  d=sc.get.rank_genes_groups_df(a,group=None);d['symbol']=d.names.map(a.var.symbol);d['run']=key;top.append(d)
  pd.crosstab(a.obs[key],a.obs['sample']).to_csv(run/f'tables/{key}_sample.tsv',sep='\t')
  for cl in a.obs[key].cat.categories:
   mask=a.obs[key].eq(cl).to_numpy()
   for sym in panel:
    gi=np.flatnonzero(a.var.symbol.eq(sym).to_numpy())
    if not len(gi):continue
    vals=np.asarray(a.layers['counts'][mask][:,gi].sum(axis=1)).ravel()
    rows.append([key,cl,sym,int(mask.sum()),int((vals>0).sum()),100*(vals>0).mean()])
  sc.pl.umap(a,color=[key,'sample'],show=False,frameon=False)
  plt.savefig(run/f'figures/{key}.png',dpi=300,bbox_inches='tight');plt.close('all')
 a.obsm[f'X_umap_K{k}']=a.obsm['X_umap'].copy()
pd.concat(top).to_csv(run/'tables/top50.tsv',sep='\t',index=False)
pd.DataFrame(rows,columns=['run','cluster','symbol','cells','detected','pct']).to_csv(run/'tables/panel.tsv',sep='\t',index=False)
pd.DataFrame([[adjusted_rand_score(a.obs[i],a.obs[j]) for j in keys] for i in keys],index=keys,columns=keys).to_csv(run/'tables/ARI.tsv',sep='\t')
c=a.obs.copy()
for sym in panel:
 gi=np.flatnonzero(a.var.symbol.eq(sym).to_numpy())
 c[sym]=np.asarray(a.layers['counts'][:,gi].sum(axis=1)).ravel() if len(gi) else np.nan
c.to_csv(run/'tables/cell_evidence.tsv.gz',sep='\t')
a.uns['dataset']='GSE274372';a.write_h5ad(run/'objects/T_refined.h5ad',compression='gzip')
(run/'tables/methods.txt').write_text('Stage4 clusters1/2 only; exclude Scrublet score>=0.2 as sensitivity-defined working singlets, not confirmed doublets. No selection on Cd4/Cd8 expression. HVG1500 with sample key; PCA30, neighbors10/20 using20 PCs; Leiden0.4/0.8/1.2. No Harmony or day regression. KxxRyy labels: neighbors xx and resolution yy/10. Each dot is one candidate T cell; sample Dxx_Ry: day and replicate label. All solutions retained for marker and sample-support review, not selected by Mre11a significance. ARI: adjusted Rand index. Missing panel features are NA, not zero.\n')
print('GSE274372_STAGE5_ANALYSIS_PASS')
