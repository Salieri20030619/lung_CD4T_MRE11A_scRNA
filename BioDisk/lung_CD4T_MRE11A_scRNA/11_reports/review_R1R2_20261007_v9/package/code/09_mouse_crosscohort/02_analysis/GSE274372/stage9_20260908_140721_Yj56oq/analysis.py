import os
from pathlib import Path
import numpy as np,pandas as pd,scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(os.environ['ROOT']);run=Path(os.environ['RUN'])
for d in ['tables','figures','objects']:(run/d).mkdir()
a=sc.read_h5ad(root/'02_analysis/GSE274372/stage8_20260908_121103_A80K5F/objects/CD4_author_candidates.h5ad')
assert a.n_obs==888
# Recompute expression and geometry within all author candidates; preserve raw counts.
a.X=a.layers['counts'].copy()
sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a)
sc.pp.highly_variable_genes(a,n_top_genes=1500,flavor='seurat',batch_key='sample')
a.var.to_csv(run/'tables/features.tsv.gz',sep='\t')
b=a[:,a.var.highly_variable].copy()
sc.pp.scale(b,max_value=10);sc.tl.pca(b,n_comps=30,svd_solver='arpack',random_state=42)
a.obsm['X_pca']=b.obsm['X_pca'].copy()
sc.pp.neighbors(a,n_neighbors=15,n_pcs=20,random_state=42)
sc.tl.umap(a,random_state=42)
for res,key in [(0.4,'R04'),(0.8,'R08'),(1.2,'R12')]:
 sc.tl.leiden(a,resolution=res,key_added=key,random_state=42)
 sc.tl.rank_genes_groups(a,groupby=key,method='wilcoxon',use_raw=False,key_added='rank_'+key)
 tables=[]
 for c in a.obs[key].cat.categories:
  t=sc.get.rank_genes_groups_df(a,group=c,key='rank_'+key).head(50)
  t['cluster']=c;t['symbol']=t['names'].map(a.var['symbol']);tables.append(t)
 pd.concat(tables,ignore_index=True).to_csv(run/f'tables/{key}_top50.tsv',sep='\t',index=False)
 pd.crosstab(a.obs[key],a.obs['sample']).to_csv(run/f'tables/{key}_sample.tsv',sep='\t')
 pd.crosstab(a.obs[key],a.obs['author_identity']).to_csv(run/f'tables/{key}_author.tsv',sep='\t')
 a.obs.groupby(key,observed=True)[['total_counts','n_genes','pct_mt','dbl_score']].median().to_csv(run/f'tables/{key}_QC.tsv',sep='\t')
symbols=a.var['symbol'].astype('string').fillna('').to_numpy()
markers=['Cd4','Cd8a','Cd8b1','Zbtb7b','Runx3','Cd3d','Cd3e','Cd3g','Bcl11b','Themis','Skap1','Tcf7','Lef1','Ccr7','Sell','Il7r','Foxp3','Il2ra','Ctla4','Ikzf2','Icos','Nkg7','Ccl5','Prf1','Sftpc','Scgb1a1','Lyz2','Tyrobp','Lst1','Cd79a','Ms4a1','Pecam1','Col1a1','Mre11a']
e=a.obs.copy();rows=[]
for g in markers:
 ix=np.flatnonzero(symbols==g)
 values=np.asarray(a.layers['counts'][:,ix].sum(axis=1)).ravel() if len(ix) else np.full(a.n_obs,np.nan)
 e[g]=values
 for key in ['R04','R08','R12']:
  for c in a.obs[key].cat.categories:
   v=values[a.obs[key].eq(c).to_numpy()]
   rows.append([key,c,g,len(v),len(ix),100*np.mean(v>0) if len(ix) else np.nan,float(np.mean(v)) if len(ix) else np.nan])
pd.DataFrame(rows,columns=['resolution','cluster','gene','cells','feature_matches','pct','mean_UMI']).to_csv(run/'tables/marker_panel.tsv',sep='\t',index=False)
e.to_csv(run/'tables/cell_evidence.tsv.gz',sep='\t')
a.obs['Author']=a.obs.author_identity.astype(str).map({'Cd4+ T cells':'CD4','Regulatory T cells':'Treg'})
a.obs['DayLabel']=a.obs.day.map(lambda x:'D'+str(int(x))).astype('category')
for color in ['R04','R08','R12','Author','DayLabel','total_counts']:
 sc.pl.umap(a,color=color,show=False,title={'total_counts':'UMI'}.get(color,color),frameon=False)
 plt.savefig(run/f'figures/UMAP_{color}.png',dpi=300,bbox_inches='tight');plt.close('all')
a.write_h5ad(run/'objects/CD4_author_reclustered.h5ad',compression='gzip')
(run/'tables/methods.txt').write_text('888 author CD4/Treg candidates; same UMI>100 QC as prior stages. Raw counts preserved. Log1p CP10k, 1500 batch-aware HVGs, scaled PCA30, neighbors15/PC20, Leiden0.4/0.8/1.2, seed42. Author labels shown for comparison only. Cluster marker p-values describe clusters, not biological time effects. No automatic exclusions. R04/R08/R12=Leiden resolutions; Author=author annotation; CD4=author Cd4+ T cells; Treg=author Regulatory T cells; DayLabel=collection day; UMI=total UMI counts.\n')
print('GSE274372_STAGE9_ANALYSIS_PASS',flush=True)
