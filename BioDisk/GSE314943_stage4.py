import sys, os, tarfile, json
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import harmonypy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
source=root/'02_analysis/GSE314943/stage3_20260905_214228_VPtd0r/objects/TNK_reclustered.h5ad'
out=Path(sys.argv[1])
out.mkdir(parents=True,exist_ok=True)
for folder in ['tables','figures','objects']:
    (out/folder).mkdir(exist_ok=True)
print('SOURCE='+str(source),flush=True)
b=sc.read_h5ad(source)
assert b.n_obs==6303 and 'counts' in b.layers
assert b.obs['sample'].value_counts().to_dict()=={'Saline_1':1073,'Saline_2':1338,'Saline_3':665,'Silica_1':1223,'Silica_2':963,'Silica_3':1041}
sc.settings.n_jobs=8
sc.settings.verbosity=2
b.obs['OldT']=b.obs.T06.astype(str)
b.raw=None
b.uns.clear(); b.obsm.clear(); b.obsp.clear(); b.varm.clear()
b.X=b.layers['counts'].copy()
sc.pp.normalize_total(b,target_sum=1e4)
sc.pp.log1p(b)
sc.pp.highly_variable_genes(b,n_top_genes=3000,flavor='seurat',batch_key='sample')
# Feature sensitivity only: no ambient count correction and no gene deletion.
background=['Sftpc','Sftpa1','Sftpb','Sftpd','Scgb1a1','Scgb3a2',
            'S100a8','S100a9','Lyz2','Igkc','Ighm','Iglc1','Iglc2','Iglc3']
symbols=b.var.symbol.astype(str)
excluded=symbols.isin(background)|symbols.str.match(r'(?i)^(mt-|Rpl|Rps|Hba|Hbb)')
b.var['embedding_excluded']=excluded
b.var['hvg_used']=b.var.highly_variable & ~excluded
b.var.to_csv(out/'tables/embedding_genes.tsv.gz',sep='\t')
assert b.var.hvg_used.sum()>500
h=b[:,b.var.hvg_used].copy()
h.layers.clear()
sc.pp.scale(h,max_value=10)
sc.tl.pca(h,n_comps=40,svd_solver='arpack',random_state=42)
b.obsm['X_pca']=h.obsm['X_pca'].copy()
del h

def plot(field,stem):
    fig,ax=plt.subplots(figsize=(10,7))
    sc.pl.umap(b,color=field,ax=ax,show=False,size=8,frameon=False,title=field)
    for ext in ['png','pdf']:
        fig.savefig(out/'figures'/f'{stem}.{ext}',dpi=320,bbox_inches='tight')
    plt.close(fig)

panel=['Cd3d','Cd3e','Trac','Trbc1','Trbc2','Cd4','Cd8a','Cd8b1',
       'Il7r','Tcf7','Lef1','Ccr7','Foxp3','Il2ra','Ctla4','Ikzf2',
       'Tbx21','Ifng','Gata3','Il4','Il13','Rorc','Il17a','Bcl6','Cxcr5',
       'Pdcd1','Tox','Tigit','Nkg7','Ncr1','Klrd1','Prf1','Gzmb',
       'Trdc','Trgc1','Trgc2','Tcrg-C1','Zbtb16','Il1rl1','Rora','Kit',
       'Mki67','Top2a','Cd79a','Ms4a1','Csf1r','Flt3','Epcam','Sftpc',
       'Scgb1a1','Lyz2','Mre11a']
symbol_array=b.var.symbol.to_numpy()
diag=[]
for mode,rep,key in [('Raw','X_pca','R06'),('Harmony','X_pca_harmony','H06')]:
    if mode=='Harmony':
        fit=harmonypy.run_harmony(b.obsm['X_pca'],b.obs,['sample'],
              random_state=42,max_iter_harmony=30)
        corrected=np.asarray(fit.Z_corr)
        if corrected.shape==b.obsm['X_pca'].T.shape:
            corrected=corrected.T
        b.obsm['X_pca_harmony']=corrected
        assert b.obsm[rep].shape==b.obsm['X_pca'].shape
        assert np.isfinite(corrected).all()
        del fit
    sc.pp.neighbors(b,use_rep=rep,n_pcs=30,n_neighbors=20,random_state=42)
    sc.tl.leiden(b,resolution=0.6,key_added=key,random_state=42,flavor='leidenalg')
    sc.tl.umap(b,random_state=42)
    b.obsm['X_umap_'+mode]=b.obsm['X_umap'].copy()
    for field in [key,'sample','group','OldT']:
        plot(field,mode+'_'+field)
    sc.tl.rank_genes_groups(b,groupby=key,method='wilcoxon',use_raw=False,n_genes=100,pts=True)
    m=sc.get.rank_genes_groups_df(b,group=None)
    m['symbol']=m.names.map(b.var.symbol)
    m.to_csv(out/'tables'/f'{mode}_top100.tsv',sep='\t',index=False)
    cross=pd.crosstab(b.obs[key],b.obs['sample'])
    cross.to_csv(out/'tables'/f'{mode}_sample.tsv',sep='\t')
    pd.crosstab(b.obs[key],b.obs.OldT).to_csv(out/'tables'/f'{mode}_origin.tsv',sep='\t')
    ev,targets=[],[]
    for cl in b.obs[key].cat.categories:
        sub=b[b.obs[key]==cl]
        freq=cross.loc[cl].to_numpy(dtype=float)
        freq=freq/freq.sum()
        diag.append(dict(mode=mode,cluster=cl,cells=sub.n_obs,
              max_sample_frac=float(freq.max()),samples_ge20=int((cross.loc[cl]>=20).sum()),
              entropy=float(-(freq[freq>0]*np.log(freq[freq>0])).sum()/np.log(6))))
        for gene in panel:
            pos=np.flatnonzero(symbol_array==gene)
            if not len(pos):continue
            v=np.asarray(sub.layers['counts'][:,pos].sum(axis=1)).ravel()
            ev.append(dict(cluster=cl,gene=gene,n=sub.n_obs,pct=100*float((v>0).mean()),
                pct_ge2=100*float((v>=2).mean()),mean_umi=float(v.mean()),
                mean_log=float(sub.X[:,pos].sum()/sub.n_obs)))
        for sample in sorted(b.obs['sample'].unique()):
            z=sub[sub.obs['sample']==sample]
            pos=np.flatnonzero(b.var_names=='ENSMUSG00000031928')
            assert len(pos)==1
            v=np.asarray(z.layers['counts'][:,pos].sum(axis=1)).ravel()
            targets.append(dict(cluster=cl,sample=sample,cells=z.n_obs,
                 mre_pos=int((v>0).sum()),mre_umi=int(v.sum()),
                 mre_pct=100*float((v>0).mean()) if z.n_obs else np.nan))
    pd.DataFrame(ev).to_csv(out/'tables'/f'{mode}_panel.tsv',sep='\t',index=False)
    pd.DataFrame(targets).to_csv(out/'tables'/f'{mode}_Mre11a.tsv',sep='\t',index=False)
    print(mode+'_DONE',flush=True)
pd.DataFrame(diag).to_csv(out/'tables/sample_mixing.tsv',sep='\t',index=False)
pd.crosstab(b.obs.R06,b.obs.H06).to_csv(out/'tables/Raw_Harmony.tsv',sep='\t')
meta=b.obs.copy()
for mode in ['Raw','Harmony']:
    meta[[mode+'_UMAP1',mode+'_UMAP2']]=b.obsm['X_umap_'+mode]
meta.to_csv(out/'tables/cell_metadata.tsv.gz',sep='\t')
b.write_h5ad(out/'objects/TNK_harmony.h5ad',compression='gzip')
with open(out/'tables/methods.txt','w') as f:
    f.write('Source: '+str(source)+'\nOutput: '+str(out)+'\n')
    f.write('Same 6303 cells and same PCA features in both arms. Raw=feature-screened uncorrected PCA; Harmony=sample correction of the same PCA.\n')
    f.write('No group regression, no count modification, no ambient decontamination. Only embedding excludes specified background, ribosomal, mitochondrial and hemoglobin genes.\n')
    f.write('Raw and Harmony UMAP panels show R06/H06 (Leiden resolution 0.6), sample, group, OldT (stage3 cluster).\n')
    f.write('Mixing metrics are diagnostics; maximizing mixing is not the criterion for retaining biological identities.\n')
    f.write('CD4/CD8/gamma-delta/NK/ILC identities require concordant lineage markers after integration. Mre11a is not used as an inclusion criterion.\n')
archive=root/'tmp/GSE314943_stage4_return.tar.gz'
archive.parent.mkdir(parents=True,exist_ok=True)
with tarfile.open(str(archive)+'.part','w:gz') as t:
    t.add(out/'tables',arcname='tables');t.add(out/'figures',arcname='figures')
os.replace(str(archive)+'.part',archive)
print('RETURN_PACKAGE='+str(archive),flush=True)
print('GSE314943_STAGE4_PASS',flush=True)
