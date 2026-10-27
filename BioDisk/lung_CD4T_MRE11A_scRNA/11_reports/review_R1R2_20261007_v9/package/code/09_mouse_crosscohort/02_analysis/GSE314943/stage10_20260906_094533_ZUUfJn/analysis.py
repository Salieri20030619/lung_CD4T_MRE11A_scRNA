import sys,os,json,tarfile
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import harmonypy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
source=root/'02_analysis/GSE314943/stage3_20260905_214228_VPtd0r/objects/whole_annotated.h5ad'
out=Path(sys.argv[1])
for d in ['tables','figures','objects']:(out/d).mkdir(exist_ok=True)
sc.settings.n_jobs=8;sc.settings.verbosity=2

def tab(t,n,index=False):t.to_csv(out/'tables'/n,sep='\t',index=index)
def figsave(fig,n,caption):
    for ext in ['png','pdf']:fig.savefig(out/'figures'/f'{n}.{ext}',dpi=350,bbox_inches='tight')
    plt.close(fig);(out/'figures'/f'{n}_caption.txt').write_text(caption+'\n')
print('SOURCE='+str(source),flush=True)
a=sc.read_h5ad(source)
expected={'Saline_1':5985,'Saline_2':6811,'Saline_3':4416,'Silica_1':6372,'Silica_2':6643,'Silica_3':6928}
assert a.n_obs==37155 and a.n_vars==29809 and 'counts' in a.layers
assert a.obs['sample'].value_counts().to_dict()==expected
selected=['8','10','12','13','14','28']
labels={'8':'Mac','10':'DC','12':'Mono','13':'Mono','14':'Mac','28':'Cycling'}
sel=a.obs.C06.astype(str).isin(selected)
tab(pd.crosstab(a.obs.C06.astype(str),a.obs['sample']),'whole_cluster_sample.tsv',True)
b=a[sel].copy();del a
b.obs['Origin']=b.obs.C06.astype(str)
b.obs['BroadPrior']=b.obs.Origin.map(labels).astype('category')
assert b.obs['sample'].nunique()==6 and b.n_obs>1000
print('MYELOID_CANDIDATES='+str(b.n_obs),flush=True)
tab(pd.crosstab(b.obs.Origin,b.obs['sample']),'candidate_sample.tsv',True)
b.raw=None;b.uns.clear();b.obsm.clear();b.obsp.clear();b.varm.clear()
b.X=b.layers['counts'].copy();sc.pp.normalize_total(b,target_sum=1e4);sc.pp.log1p(b)
sc.pp.highly_variable_genes(b,n_top_genes=3000,flavor='seurat',batch_key='sample')
sym=b.var.symbol.astype(str)
# Embedding exclusions are lineage-extraneous transcripts and technical families;
# retain genuine myeloid markers such as Lyz2 and S100a8/9 in this compartment.
background=['Sftpc','Sftpa1','Sftpb','Sftpd','Scgb1a1','Scgb3a2','Jchain']
exclude=sym.isin(background)|sym.str.match(r'(?i)^(mt-|Rpl|Rps|Hba|Hbb|Igh|Igl|Igk)')
b.var['embedding_excluded']=exclude
b.var['hvg_used']=b.var.highly_variable & ~exclude
assert b.var.hvg_used.sum()>500
tab(b.var,'embedding_features.tsv.gz',True)
h=b[:,b.var.hvg_used].copy();h.layers.clear();sc.pp.scale(h,max_value=10)
sc.tl.pca(h,n_comps=40,random_state=42,svd_solver='arpack');b.obsm['X_pca']=h.obsm['X_pca'].copy();del h
panel=['Ptprc','Lyz2','Csf1r','Adgre1','Fcgr1','Mertk','Cd68','Aif1','Tyrobp',
 'Siglecf','Marco','Pparg','Car4','Fabp4','Chil3','Mrc1','Folr2','Lyve1','Timd4',
 'C1qa','C1qb','C1qc','Apoe','Trem2','Spp1','Gpnmb','Lgals3','Lpl','Itgax',
 'Ly6c2','Ccr2','Plac8','Chil1','Nr4a1','Fcgr3','Ace','Csf3r','S100a8','S100a9','Retnlg',
 'Flt3','Zbtb46','Xcr1','Clec9a','Batf3','Clec10a','Cd209a','Sirpa','H2-Ab1','Cd74',
 'Siglech','Bst2','Tcf4','Mki67','Top2a','Cd3d','Cd3e','Cd79a','Ms4a1','Epcam','Sftpc','Scgb1a1',
 'Mre11a','Nlrp3','Pycard','Casp1','Casp4','Gsdmd','Il1b','Tnf','Il6','Tgfb1','Spp1','Ccl2','Ccl3','Ccl4','Cxcl9','Cxcl10']
panel=list(dict.fromkeys(panel));pos={g:np.flatnonzero(sym.to_numpy()==g) for g in panel}
missing=[g for g in panel if not len(pos[g])];tab(pd.DataFrame({'missing_symbol':missing}),'missing_markers.tsv')
mixing=[]
for mode,rep,key in [('Raw','X_pca','R06'),('Harmony','X_pca_harmony','M06')]:
    if mode=='Harmony':
        ho=harmonypy.run_harmony(b.obsm['X_pca'],b.obs,['sample'],random_state=42,max_iter_harmony=30)
        z=np.asarray(ho.Z_corr)
        if z.shape==b.obsm['X_pca'].T.shape:z=z.T
        assert z.shape==b.obsm['X_pca'].shape and np.isfinite(z).all()
        b.obsm[rep]=z;del ho
    sc.pp.neighbors(b,use_rep=rep,n_pcs=30,n_neighbors=20,random_state=42)
    sc.tl.leiden(b,resolution=.6,key_added=key,random_state=42,flavor='leidenalg')
    if mode=='Harmony':
        sc.tl.leiden(b,resolution=.4,key_added='M04',random_state=42,flavor='leidenalg')
        sc.tl.leiden(b,resolution=.8,key_added='M08',random_state=42,flavor='leidenalg')
    sc.tl.umap(b,random_state=42);b.obsm['X_umap_'+mode]=b.obsm['X_umap'].copy()
    for field in [key,'sample','group','BroadPrior']:
        fig,ax=plt.subplots(figsize=(10,7));sc.pl.umap(b,color=field,ax=ax,show=False,size=6,frameon=False,title=field)
        figsave(fig,mode+'_'+field,'Candidate myeloid compartment. Raw=uncorrected PCA; Harmony=sample-adjusted PCA using the same features and cells. Mac=macrophage; Mono=monocyte; DC=dendritic cell; Cycling=proliferating cells of unresolved lineage. Prior labels are inherited broad annotations, not final macrophage subtypes.')
    sc.tl.rank_genes_groups(b,groupby=key,method='wilcoxon',use_raw=False,n_genes=100,pts=True)
    top=sc.get.rank_genes_groups_df(b,group=None);top['symbol']=top.names.map(b.var.symbol)
    tab(top,mode+'_top100.tsv')
    cross=pd.crosstab(b.obs[key],b.obs['sample']);tab(cross,mode+'_sample.tsv',True)
    tab(pd.crosstab(b.obs[key],b.obs.Origin),mode+'_origin.tsv',True)
    evidence=[];by_sample=[]
    for cl in b.obs[key].cat.categories:
        sub=b[b.obs[key]==cl]
        f=cross.loc[cl].values.astype(float);f/=f.sum()
        mixing.append(dict(mode=mode,cluster=cl,cells=sub.n_obs,max_sample_frac=f.max(),samples_ge20=int((cross.loc[cl]>=20).sum())))
        for g in panel:
            if not len(pos[g]):continue
            v=np.asarray(sub.layers['counts'][:,pos[g]].sum(axis=1)).ravel()
            evidence.append(dict(cluster=cl,gene=g,cells=sub.n_obs,pct=100*(v>0).mean(),pct_ge2=100*(v>=2).mean(),mean_umi=v.mean()))
        if mode=='Harmony':
            for sample in sorted(expected):
                ss=sub[sub.obs['sample']==sample]
                for g in panel:
                    if not len(pos[g]):continue
                    v=np.asarray(ss.layers['counts'][:,pos[g]].sum(axis=1)).ravel()
                    by_sample.append(dict(cluster=cl,sample=sample,gene=g,cells=ss.n_obs,pct=100*(v>0).mean() if ss.n_obs else np.nan,mean_umi=v.mean() if ss.n_obs else np.nan))
    ev=pd.DataFrame(evidence);tab(ev,mode+'_panel.tsv')
    if mode=='Harmony':tab(pd.DataFrame(by_sample),'Harmony_panel_by_sample.tsv')
    # Gene labels below are short; no arbitrary subtype assignment from single markers.
    show=['Csf1r','Mertk','Adgre1','Siglecf','Pparg','Fabp4','Folr2','Lyve1','C1qa','Trem2','Spp1','Gpnmb','Ly6c2','Ccr2','Nr4a1','Flt3','Xcr1','Clec9a','Cd209a','Siglech','Csf3r','Cd3d','Cd79a','Sftpc','Mki67']
    show=[g for g in show if len(pos[g])];clusters=list(b.obs[key].cat.categories)
    z=ev[ev.gene.isin(show)].copy();z['color']=np.log1p(z.mean_umi)
    fig,ax=plt.subplots(figsize=(13,max(5,len(clusters)*.28)))
    im=ax.scatter(z.gene.map({g:i for i,g in enumerate(show)}),z.cluster.map({c:i for i,c in enumerate(clusters)}),s=z.pct*1.4,c=z.color,cmap='viridis',vmin=0)
    ax.set_xticks(range(len(show)),show,rotation=90);ax.set_yticks(range(len(clusters)),clusters);ax.set_ylabel('Cluster');fig.colorbar(im,ax=ax,label='log1p UMI')
    for pct in [10,50,100]:ax.scatter([],[],s=pct*1.4,c='grey',label=str(pct)+'%')
    ax.legend(title='Detected',bbox_to_anchor=(1.15,1),frameon=False)
    figsave(fig,mode+'_markers','Marker dot plot. Dot area=percentage of cells with at least one raw UMI; color=log1p mean raw UMI. These markers distinguish macrophage, monocyte, dendritic, granulocyte and contaminating lineages. Expression panels support subsequent annotation; raw UMI color can depend on depth.')
    print(mode+'_PASS',flush=True)
tab(pd.DataFrame(mixing),'sample_mixing.tsv')
tab(pd.crosstab(b.obs.R06,b.obs.M06),'Raw_Harmony.tsv',True)
for key in ['M04','M08']:tab(pd.crosstab(b.obs.M06,b.obs[key]),'M06_'+key+'.tsv',True)
metadata=b.obs.copy()
for mode in ['Raw','Harmony']:metadata[[mode+'_UMAP1',mode+'_UMAP2']]=b.obsm['X_umap_'+mode]
tab(metadata,'myeloid_cells.tsv.gz',True)
b.write_h5ad(out/'objects/myeloid_candidates.h5ad',compression='gzip')
(out/'tables/methods.txt').write_text('Source: '+str(source)+'\nOutput: '+str(out)+'\n'
 'Candidates: C06 clusters8,10,12,13,14,28 (macrophages, monocytes, dendritic cells and unresolved cycling cells). Do not equate all candidates with macrophages. Other whole-lung clusters remain in source.\n'
 'Same sample-aware HVG PCA features in Raw and Harmony arms. Remove epithelial background, immunoglobulin, ribosomal, mitochondrial and hemoglobin features only from embedding. Retain Lyz2/S100a8/S100a9 as legitimate myeloid markers. Original counts preserved; no decontamination.\n'
 'Harmony sample correction, 30 PCs, neighbors20; primary clustering resolution0.6 with 0.4 and0.8 comparisons. Marker tests are annotation aids, not condition DEG.\n'
 'Final macrophage inclusion requires concordant lineage markers and checking sample-level support. No cell communication significance is inferred at this stage.\n')
archive=root/'tmp/GSE314943_stage10_return.tar.gz'
with tarfile.open(str(archive)+'.part','w:gz') as tf:
    for d in ['tables','figures']:tf.add(out/d,arcname=d)
os.replace(str(archive)+'.part',archive)
print('GSE314943_STAGE10_PASS',flush=True)
