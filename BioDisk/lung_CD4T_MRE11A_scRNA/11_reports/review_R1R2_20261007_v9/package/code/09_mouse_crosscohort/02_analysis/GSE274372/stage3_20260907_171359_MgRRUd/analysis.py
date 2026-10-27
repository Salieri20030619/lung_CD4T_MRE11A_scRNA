import os
from pathlib import Path
import numpy as np,pandas as pd,scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(os.environ['ROOT']);run=Path(os.environ['RUN'])
for sub in ['tables','figures','objects']:(run/sub).mkdir()
a=sc.read_h5ad(root/'02_analysis/GSE274372/stage2_20260907_105825_0DB4hC/objects/whole_screened.h5ad')
labels={0:'AT2',1:'Myel',2:'Endoth',3:'AT1',4:'Club',5:'Fib',6:'AT2',7:'T',8:'B',9:'MonoDC',10:'Ciliated',11:'Mac',12:'Endoth',13:'Unresolved',14:'NKLike',15:'MesoLike',16:'LEC',17:'Neut',18:'MuralLike',19:'PlasmaLike'}
a.obs['Broad']=a.obs.C06.astype(str).map(lambda x:labels[int(x)]).astype('category')
pd.DataFrame({'C06':list(labels),'Broad':list(labels.values())}).to_csv(run/'tables/annotation.tsv',sep='\t',index=False)
pd.crosstab(a.obs.Broad,a.obs['sample']).to_csv(run/'tables/whole_sample.tsv',sep='\t')
sc.pl.umap(a,color='Broad',show=False,frameon=False,palette='tab20',title='Whole lung')
plt.savefig(run/'figures/whole_Broad.png',dpi=350,bbox_inches='tight');plt.close('all')
a.write_h5ad(run/'objects/whole_annotated.h5ad',compression='gzip')
b=a[a.obs.C06.astype(str).isin(['7','14'])].copy()
assert b.n_obs==2169
# Keep flagged cells to inspect their cluster distribution before exclusion.
b.X=b.layers['counts'].copy();b.uns.clear()
for key in list(b.obsp.keys()):del b.obsp[key]
for key in list(b.obsm.keys()):del b.obsm[key]
sc.pp.normalize_total(b,target_sum=1e4);sc.pp.log1p(b)
sc.pp.highly_variable_genes(b,n_top_genes=2000,flavor='seurat',batch_key='sample')
v=b[:,b.var.highly_variable].copy();sc.pp.scale(v,max_value=10)
sc.tl.pca(v,n_comps=30,svd_solver='arpack',random_state=42)
b.obsm['X_pca']=v.obsm['X_pca'].copy()
sc.pp.neighbors(b,n_neighbors=15,n_pcs=20,random_state=42);sc.tl.umap(b,random_state=42)
for res,key in [(0.3,'T03'),(0.5,'T05'),(0.8,'T08')]:sc.tl.leiden(b,resolution=res,key_added=key,random_state=42)
sc.tl.rank_genes_groups(b,'T05',method='wilcoxon',use_raw=False,n_genes=100)
top=sc.get.rank_genes_groups_df(b,group=None);top['symbol']=top.names.map(b.var.symbol)
top.to_csv(run/'tables/TNK_top100.tsv',sep='\t',index=False)
panel=['Cd3d','Cd3e','Cd3g','Trac','Trbc1','Trbc2','Cd247','Skap1','Themis','Bcl11b','Cd4','Cd8a','Cd8b1','Tcf7','Lef1','Ccr7','Sell','Il7r','Foxp3','Il2ra','Ctla4','Ikzf2','Tnfrsf4','Icos','Pdcd1','Gata3','Rorc','Il17a','Tbx21','Ifng','Nkg7','Klrd1','Klrk1','Eomes','Prf1','Gzma','Gzmb','Ccl5','Trdc','Trgc1','Trgc2','Zbtb16','Mki67','Top2a','Cd79a','Ms4a1','Csf1r','Lyz2','Sftpc','Scgb1a1','Pecam1','Mre11a']
rows=[]
for cl in b.obs.T05.cat.categories:
 mask=b.obs.T05.eq(cl).to_numpy();c=b.layers['counts'][mask]
 for sym in panel:
  gi=np.flatnonzero(b.var.symbol.eq(sym).to_numpy())
  if not len(gi):continue
  vals=np.asarray(c[:,gi].sum(axis=1)).ravel();expr=np.asarray(b.X[mask][:,gi].sum(axis=1)).ravel()
  rows.append([cl,sym,int(mask.sum()),int((vals>0).sum()),100*(vals>0).mean(),expr.mean()])
pd.DataFrame(rows,columns=['cluster','symbol','cells','detected','pct','mean_log']).to_csv(run/'tables/TNK_panel.tsv',sep='\t',index=False)
pd.crosstab(b.obs.T05,b.obs['sample']).to_csv(run/'tables/TNK_sample.tsv',sep='\t')
pd.crosstab(b.obs.T05,b.obs.C06).to_csv(run/'tables/TNK_origin.tsv',sep='\t')
b.obs.groupby('T05',observed=True)[['total_counts','n_genes','pct_mt','dbl_score','dbl_020']].mean().to_csv(run/'tables/TNK_QC.tsv',sep='\t')
b.obs.to_csv(run/'tables/TNK_cells.tsv.gz',sep='\t')
for col in ['T05','sample','Day','dbl_score']:
 sc.pl.umap(b,color=col,show=False,frameon=False,title=col)
 plt.savefig(run/f'figures/TNK_{col}.png',dpi=350,bbox_inches='tight');plt.close('all')
b.uns['dataset']='GSE274372';b.write_h5ad(run/'objects/TNK_reclustered.h5ad',compression='gzip')
(run/'tables/captions.txt').write_text('Broad: preliminary whole-lung lineage annotation. AT1/AT2: alveolar type I/II; Myel: myeloid; Endoth: vascular endothelial; Club: club epithelial; Fib: fibroblast; T/B: T/B lymphocytes; MonoDC: monocyte/dendritic unresolved; Ciliated: ciliated epithelial; Mac: macrophage; Unresolved: unresolved identity; NKLike: NK-like; MesoLike: mesothelial-like; LEC: lymphatic endothelial; Neut: neutrophil; MuralLike: mural-like; PlasmaLike: plasma-like. Like labels are provisional. T05: T/NK Leiden resolution0.5; Day: D00/D07/D28/D56; sample: day and replicate label. dbl_score: Scrublet score. Scores >=0.2 remain flagged, not removed. Low lineage-incongruent transcripts may reflect ambient RNA; evaluate coexpression before doublet decisions. No CD4 identity assigned in this stage.\n')
print(f'TNK_CELLS={b.n_obs};TNK_CLUSTERS={b.obs.T05.nunique()}')
print('GSE274372_STAGE3_ANALYSIS_PASS')
