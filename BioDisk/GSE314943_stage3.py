import os, sys, json, tarfile
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root = Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
out = Path(sys.argv[1])
for name in ['tables', 'figures', 'objects']:
    (out/name).mkdir(exist_ok=True)
sources = sorted((root/'02_analysis/GSE314943').glob('stage2_*/objects/GSE314943_whole_lung.h5ad'), key=lambda p:p.stat().st_mtime)
if not sources:
    raise SystemExit('STAGE2_OBJECT_NOT_FOUND')
source = sources[-1]
print('SOURCE=' + str(source), flush=True)
a = sc.read_h5ad(source)
assert a.n_obs == 37155 and a.n_vars == 29809, 'SOURCE_DIMENSION_MISMATCH'
assert set(a.obs.C06.astype(str)) == {str(i) for i in range(32)}
expected = {'Saline_1':5985,'Saline_2':6811,'Saline_3':4416,
            'Silica_1':6372,'Silica_2':6643,'Silica_3':6928}
assert a.obs['sample'].value_counts().to_dict() == expected
assert 'counts' in a.layers
sc.settings.n_jobs = 8
sc.settings.verbosity = 2
labels = ['B','Neut','Fib','B','T','T','AT2','NeutLike',
          'Mac','Endoth','DC','T','Mono','Mono','Mac','AT2Like',
          'Neut','NKLike','Plasma','Club','ILC2Like','SMC','Peri',
          'Ciliated','AT1','Endoth','Fib','LEC','Unresolved','BasoMast',
          'Neut','Platelet']
evidence = ['Cd79a/Ebf1','S100a8/S100a9/Csf3r','Col1a1/Pdgfra',
 'Cd79a/Ebf1','Cd3d/Cd3e/Trbc2','Cd3d/Trbc2/Lef1','Sftpc/Abca3',
 'S100a8/S100a9/Clec4d','Csf1r/Trem2/Cd68','Pecam1/Cdh5',
 'Flt3/Itgax/MHC-II','Cd3d/Cd3e/Tcrg-C1','Ly6c2/Ccr2/Csf1r',
 'Csf1r/Nr4a1/Ace','Mrc1/Chil3/Lpl','Sftpa1/Sftpd/Lamp3',
 'Csf3r/Hdc/S100a8','Nkg7/Klrd1/Ccl5','Jchain/Mzb1/Xbp1',
 'Scgb1a1/Cyp2f2','Gata3/Il1rl1/Il7r','Acta2/Tagln',
 'Pdgfrb/Notch3','Foxj1/Tppp3','Ager/Hopx','Cdh5/Kdr/Car4',
 'Col1a1/Col3a1','Flt4/Ccl21a/Mmrn1','Mki67/Top2a',
 'Cpa3/Gata2/Ms4a2','S100a8/S100a9/Ngp','Pf4/Ppbp']
mapping = dict(zip(map(str,range(32)), labels))
a.obs['Broad'] = pd.Categorical(a.obs.C06.astype(str).map(mapping))
a.obs['State'] = np.where(a.obs.C06.astype(str)=='28','Cycling','Unassigned')
notes = pd.DataFrame({'C06':list(map(str,range(32))), 'Broad':labels,
                      'markers':evidence,
                      'status':['Provisional' if i in [7,15,17,20,28,29] else 'Broad-supported' for i in range(32)]})
notes.to_csv(out/'tables/whole_annotation.tsv',sep='\t',index=False)
pd.crosstab(a.obs.Broad,a.obs['sample']).to_csv(out/'tables/whole_sample.tsv',sep='\t')

def plot_umap(obj, field, filename):
    fig, ax = plt.subplots(figsize=(10,7))
    sc.pl.umap(obj,color=field,ax=ax,show=False,size=6,frameon=False,title=field)
    fig.savefig(out/'figures'/f'{filename}.png',dpi=350,bbox_inches='tight')
    fig.savefig(out/'figures'/f'{filename}.pdf',bbox_inches='tight')
    plt.close(fig)

plot_umap(a,'Broad','whole_Broad')
a.write_h5ad(out/'objects/whole_annotated.h5ad',compression='gzip')
selected = ['4','5','11','17','20','28']
b = a[a.obs.C06.astype(str).isin(selected)].copy()
assert b.n_obs == 6303
del a
# Retain counts and metadata; rebuild all embeddings within this compartment.
b.X = b.layers['counts'].copy()
b.raw = None
b.uns.clear(); b.obsm.clear(); b.obsp.clear(); b.varm.clear()
sc.pp.normalize_total(b,target_sum=1e4)
sc.pp.log1p(b)
sc.pp.highly_variable_genes(b,flavor='seurat',n_top_genes=2500,batch_key='sample')
h = b[:,b.var.highly_variable].copy()
h.layers.clear()
sc.pp.scale(h,max_value=10)
sc.tl.pca(h,n_comps=40,svd_solver='arpack',random_state=42)
b.obsm['X_pca'] = h.obsm['X_pca'].copy()
del h
sc.pp.neighbors(b,n_neighbors=20,n_pcs=30,random_state=42)
for res,key in [(0.4,'T04'),(0.6,'T06'),(0.8,'T08')]:
    sc.tl.leiden(b,resolution=res,key_added=key,random_state=42,flavor='leidenalg')
sc.tl.umap(b,random_state=42)
sc.tl.rank_genes_groups(b,groupby='T06',method='wilcoxon',use_raw=False,n_genes=100,pts=True)
d = sc.get.rank_genes_groups_df(b,group=None)
d['symbol'] = d.names.map(b.var.symbol)
d.to_csv(out/'tables/TNK_top100.tsv',sep='\t',index=False)
pd.crosstab(b.obs.T06,b.obs['sample']).to_csv(out/'tables/TNK_sample.tsv',sep='\t')
pd.crosstab(b.obs.T06,b.obs.C06).to_csv(out/'tables/TNK_origin.tsv',sep='\t')
meta = b.obs.copy()
meta[['UMAP1','UMAP2']] = b.obsm['X_umap']
for gene in ['Cd3d','Cd3e','Cd4','Cd8a','Trdc','Nkg7','Epcam','Sftpc','Scgb1a1','Lyz2','Mre11a']:
    pos = np.flatnonzero(b.var.symbol.to_numpy()==gene)
    if len(pos):
        meta[gene+'_UMI'] = np.asarray(b.layers['counts'][:,pos].sum(axis=1)).ravel()
meta.to_csv(out/'tables/TNK_cells.tsv.gz',sep='\t')
panel = ['Cd3d','Cd3e','Trac','Trbc1','Trbc2','Cd4','Cd8a','Cd8b1',
 'Il7r','Ccr7','Tcf7','Lef1','Foxp3','Il2ra','Ctla4','Ikzf2','Tbx21',
 'Ifng','Gata3','Il4','Il13','Rorc','Il17a','Il17f','Bcl6','Cxcr5','Pdcd1',
 'Tox','Lag3','Tigit','Nkg7','Klrd1','Ncr1','Gzmb','Prf1','Trdc',
 'Trgc1','Trgc2','Tcrg-C1','Zbtb16','Il1rl1','Rora','Kit','Id2',
 'Mki67','Top2a','Ms4a1','Cd79a','Lyz2','Csf1r','Flt3','Epcam',
 'Sftpc','Sftpb','Scgb1a1','S100a8','S100a9','Mre11a']
rows, targets = [], []
symbols = b.var.symbol.to_numpy()
for cl in b.obs.T06.cat.categories:
    sub = b[b.obs.T06==cl]
    for gene in panel:
        pos = np.flatnonzero(symbols==gene)
        if not len(pos):
            continue
        counts = np.asarray(sub.layers['counts'][:,pos].sum(axis=1)).ravel()
        rows.append(dict(cluster=cl,gene=gene,n=sub.n_obs,
            pct=100*float((counts>0).mean()),pct_ge2=100*float((counts>=2).mean()),
            mean_umi=float(counts.mean()),mean_log=float(sub.X[:,pos].sum()/sub.n_obs)))
    for sample in sorted(expected):
        z = sub[sub.obs['sample']==sample]
        pos = np.flatnonzero(b.var_names=='ENSMUSG00000031928')
        v = np.asarray(z.layers['counts'][:,pos].sum(axis=1)).ravel()
        targets.append(dict(cluster=cl,sample=sample,cells=z.n_obs,
                            mre_pos=int((v>0).sum()),mre_umi=int(v.sum()),
                            mre_pct=100*float((v>0).mean()) if z.n_obs else np.nan))
ev = pd.DataFrame(rows)
ev.to_csv(out/'tables/TNK_panel.tsv',sep='\t',index=False)
pd.DataFrame(targets).to_csv(out/'tables/TNK_Mre11a.tsv',sep='\t',index=False)
for field in ['T06','sample','group','C06']:
    plot_umap(b,field,'TNK_'+field)
# Dot size = detection fraction; color = mean log expression within each gene.
genes = [g for g in panel if g in set(ev.gene)]
clusters = list(b.obs.T06.cat.categories)
mean = ev.pivot(index='cluster',columns='gene',values='mean_log').reindex(index=clusters,columns=genes)
pct = ev.pivot(index='cluster',columns='gene',values='pct').reindex(index=clusters,columns=genes)
scaled = mean.divide(mean.max(axis=0).replace(0,1),axis=1)
fig,ax = plt.subplots(figsize=(19,max(5,len(clusters)*0.35)))
xx,yy = np.meshgrid(np.arange(len(genes)),np.arange(len(clusters)))
points=ax.scatter(xx.ravel(),yy.ravel(),s=pct.to_numpy().ravel()*0.8,
                  c=scaled.to_numpy().ravel(),cmap='Reds',vmin=0,vmax=1)
ax.set_xticks(np.arange(len(genes)),genes,rotation=90)
ax.set_yticks(np.arange(len(clusters)),clusters)
ax.set_ylabel('T06'); ax.invert_yaxis()
fig.colorbar(points,ax=ax,label='RelExpr',shrink=0.5)
for p,label in [(25,'25%'),(50,'50%'),(100,'100%')]:
    ax.scatter([],[],s=p*0.8,c='grey',label=label)
ax.legend(title='Pct',loc='upper left',bbox_to_anchor=(1.02,1))
fig.tight_layout()
fig.savefig(out/'figures/TNK_panel.png',dpi=300,bbox_inches='tight')
fig.savefig(out/'figures/TNK_panel.pdf',bbox_inches='tight')
plt.close(fig)
b.write_h5ad(out/'objects/TNK_reclustered.h5ad',compression='gzip')
with open(out/'tables/captions.txt','w') as f:
    f.write('Whole-lung broad annotation from stage2 markers; provisional labels retained.\n')
    f.write('B=B cell; Neut=neutrophil; Fib=fibroblast; T=T cell; AT1/AT2=alveolar type 1/2 cell; Mac=macrophage; Endoth=endothelial; DC=dendritic; Mono=monocyte; NK=natural killer; ILC2=type 2 innate lymphoid; SMC=smooth muscle; Peri=pericyte; LEC=lymphatic endothelial; BasoMast=basophil/mast unresolved. Like=provisional identity.\n')
    f.write('TNK includes stage2 C4,C5,C11,C17,C20,C28; 6303 cells. C28 included to resolve cycling lineage.\n')
    f.write('T04/T06/T08=Leiden resolutions 0.4/0.6/0.8. TNK=T/NK/ILC/cycling candidate compartment.\n')
    f.write('Dot size=Pct (percentage with >0 UMI); color=RelExpr (mean log expression divided by gene maximum across clusters).\n')
    f.write('PCA/UMAP uncorrected; no ambient RNA correction performed. Background markers are audited at >0 and >=2 UMI.\n')
    f.write('Cluster marker P values are annotation evidence, not animal-level case-control inference.\n')
with open(out/'tables/run_info.json','w') as f:
    json.dump({'source':str(source),'output':str(out),'whole_cells':37155,'TNK_cells':b.n_obs},f,indent=2)
archive = root/'tmp/GSE314943_stage3_return.tar.gz'
archive.parent.mkdir(exist_ok=True)
with tarfile.open(str(archive)+'.part','w:gz') as t:
    t.add(out/'tables',arcname='tables'); t.add(out/'figures',arcname='figures')
os.replace(str(archive)+'.part',archive)
print('RETURN_PACKAGE='+str(archive),flush=True)
print('GSE314943_STAGE3_PASS',flush=True)
