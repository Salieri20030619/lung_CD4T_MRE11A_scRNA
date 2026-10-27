import sys,os,tarfile,json,warnings
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse.csgraph import connected_components
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
src=root/'02_analysis/GSE314943/stage5_20260905_222457_a2dApT/objects/CD4_core.h5ad'
scorepath=root/'02_analysis/GSE314943/stage8_20260906_003153_6mFKmg/tables/cell_scores.tsv'
out=Path(sys.argv[1])
for d in ['tables','figures','objects']:(out/d).mkdir(exist_ok=True)
sc.settings.n_jobs=8
sc.settings.verbosity=2

def tab(df,name,index=False):df.to_csv(out/'tables'/name,sep='\t',index=index)
def figsave(fig,name,caption):
    for ext in ['png','pdf']:fig.savefig(out/'figures'/f'{name}.{ext}',dpi=350,bbox_inches='tight')
    plt.close(fig);(out/'figures'/f'{name}_caption.txt').write_text(caption+'\n')
def rho(x,y):
    x=np.asarray(x);y=np.asarray(y)
    if len(x)<3 or np.unique(x).size<2 or np.unique(y).size<2:return np.nan
    return float(spearmanr(x,y).statistic)
a=sc.read_h5ad(src)
assert a.n_obs==2311 and 'X_pca_harmony' in a.obsm
# Only conventional naive-like and activated CD4; Treg is retained in source.
a=a[a.obs.Subtype.isin(['Naive','ActCD4'])].copy()
assert a.n_obs==1734
scores=pd.read_csv(scorepath,sep='\t',index_col='cell')
assert a.obs_names.isin(scores.index).all() and scores.index.is_unique
scores=scores.loc[a.obs_names]
for col in scores:a.obs[col]=scores[col].values
# Fresh normalized expression from raw counts, avoiding stale log1p metadata.
a.X=a.layers['counts'].copy();a.uns.pop('log1p',None)
sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a)
samples=['Saline_1','Saline_2','Saline_3','Silica_1','Silica_2','Silica_3']
syms=a.var.symbol.astype(str).values
p=np.flatnonzero(np.isin(syms,['Tcf7','Lef1','Ccr7']))
assert len(set(syms[p]))==3
a.obs['RootScore']=np.asarray(a.X[:,p].mean(axis=1)).ravel()
for gene in ['Mre11a','Casp1','Casp4','Gsdmd','Gsdme','Tcf7','Ccr7','Pdcd1','Lag3']:
    pos=np.flatnonzero(syms==gene)
    if gene=='Mre11a':pos=np.flatnonzero(a.var_names=='ENSMUSG00000031928')
    if not len(pos):continue
    v=np.asarray(a.layers['counts'][:,pos].sum(axis=1)).ravel()
    a.obs[gene+'_det']=(v>0).astype(float)
    a.obs[gene+'_log']=np.asarray(a.X[:,pos].sum(axis=1)).ravel()
metrics=list(scores.columns)+['Mre11a_det','Mre11a_log','Casp1_det','Casp4_det','Gsdmd_det']
metrics=[m for m in metrics if m in a.obs]
# Predefined graph size comparison; no graph tuning on target trends.
assignments=[];runstatus=[];bins=[];relations=[];roots=[];pseudotimes={}
for k in [15,30]:
    b=a.copy()
    sc.pp.neighbors(b,n_neighbors=k,use_rep='X_pca_harmony',n_pcs=30,random_state=42)
    nc,lab=connected_components(b.obsp['connectivities'],directed=False)
    b.obs['Component']=lab
    for cell,component in zip(b.obs_names,lab):assignments.append(dict(graph=k,cell=cell,component=int(component)))
    eligible=[]
    for cl in range(nc):
        m=b.obs[lab==cl];ct=m.Subtype.value_counts();ns=m['sample'].value_counts().reindex(samples,fill_value=0)
        if ct.get('Naive',0)>=50 and ct.get('ActCD4',0)>=50 and (ns>=20).all():eligible.append((len(m),cl))
    if not eligible:
        runstatus.append(dict(graph=k,status='NO_SUPPORTED_COMPONENT',cells=0));continue
    _,chosen=max(eligible)
    # Keep its existing within-component graph, without reconnecting other components.
    b=b[b.obs.Component==chosen].copy()
    b.obs['Subtype']=b.obs.Subtype.astype('category').cat.remove_unused_categories()
    cand=np.flatnonzero(b.obs.Subtype.to_numpy()=='Naive')
    q=np.quantile(b.obs.RootScore.iloc[cand],.75)
    cand=cand[b.obs.RootScore.iloc[cand].to_numpy()>=q]
    pc=b.obsm['X_pca_harmony'][:,:30]
    center=pc[cand].mean(axis=0)
    ordered=cand[np.argsort(((pc[cand]-center)**2).sum(axis=1),kind='stable')]
    first=int(ordered[0]);different=[int(i) for i in ordered if b.obs['sample'].iloc[i]!=b.obs['sample'].iloc[first]]
    picks=[first]+different[:1]
    if len(picks)<2:
        runstatus.append(dict(graph=k,status='NO_SECOND_SAMPLE_ROOT',cells=b.n_obs));continue
    sc.tl.diffmap(b,n_comps=15,random_state=42)
    for ir,iroot in enumerate(picks,1):
        name=f'K{k}_R{ir}'
        b.uns['iroot']=iroot
        sc.tl.dpt(b,n_dcs=10,n_branchings=0)
        pt=b.obs.dpt_pseudotime.to_numpy().copy()
        assert np.isfinite(pt).all() and (pt>=0).all() and (pt<=1.00001).all()
        pseudotimes[name]=pd.Series(pt,index=b.obs_names)
        roots.append(dict(run=name,cell=b.obs_names[iroot],sample=b.obs['sample'].iloc[iroot],RootScore=float(b.obs.RootScore.iloc[iroot])))
        runstatus.append(dict(graph=k,run=name,status='SUPPORTED',cells=b.n_obs,fraction=b.n_obs/a.n_obs))
        m=b.obs.copy();m['PT']=pt;m['Bin']=np.minimum((pt*5).astype(int),4)
        tab(m,f'{name}_cells.tsv.gz',True)
        xy=b.obsm['X_umap_Harmony']
        fig,ax=plt.subplots(figsize=(7,5));im=ax.scatter(xy[:,0],xy[:,1],c=pt,cmap='viridis',s=7)
        ax.scatter(xy[iroot,0],xy[iroot,1],marker='*',s=120,c='red',edgecolors='black')
        fig.colorbar(im,ax=ax,label='DPT');ax.set(xlabel='UMAP1',ylabel='UMAP2');ax.set_xticks([]);ax.set_yticks([])
        figsave(fig,name+'_DPT','DPT=diffusion pseudotime. UMAP coordinates inherited from Stage4 Harmony. Star=root selected among high Tcf7/Lef1/Ccr7 naive-like CD4 cells, nearest their PCA centroid; second root uses another sample. This is an assumed cell-state ordering at one endpoint, not observed longitudinal time.')
        for sample in samples:
            sub=m[m['sample']==sample]
            for metric in metrics:
                relations.append(dict(run=name,sample=sample,metric=metric,cells=len(sub),rho=rho(sub.PT,sub[metric])))
            for bi in range(5):
                z=sub[sub.Bin==bi]
                for metric in metrics:
                    bins.append(dict(run=name,sample=sample,group='Saline' if sample.startswith('Saline') else 'Silica',bin=bi,
                        cells=len(z),metric=metric,PT=float(z.PT.mean()) if len(z) else np.nan,
                        value=float(z[metric].mean()) if len(z)>=10 else np.nan))
        tab(pd.crosstab(m['sample'],m.Subtype),f'{name}_sample_counts.tsv',True)
        tab(m.groupby('Subtype',observed=True).PT.agg(['count','median','mean']).reset_index(),f'{name}_state_order.tsv')
        if ir==1:b.write_h5ad(out/'objects'/f'K{k}_trajectory.h5ad',compression='gzip')
    print('GRAPH_DONE='+str(k),flush=True)
tab(pd.DataFrame(assignments),'components.tsv.gz')
tab(pd.DataFrame(runstatus),'analysis_status.tsv');tab(pd.DataFrame(roots),'roots.tsv')
if pseudotimes:
    pt=pd.DataFrame(pseudotimes);tab(pt,'pseudotime_all_runs.tsv.gz',True)
    rows=[]
    names=list(pt.columns)
    for i,n1 in enumerate(names):
        for n2 in names[i+1:]:
            z=pt[[n1,n2]].dropna();rows.append(dict(run1=n1,run2=n2,common_cells=len(z),rho=rho(z[n1],z[n2])))
    tab(pd.DataFrame(rows),'trajectory_stability.tsv')
    bt=pd.DataFrame(bins);rt=pd.DataFrame(relations)
    tab(bt,'sample_bin_trends.tsv');tab(rt,'within_sample_descriptive_rho.tsv')
    # Draw the predefined reference K15_R1 only; no switch based on favorable trends.
    ref='K15_R1'
    if ref in pseudotimes:
        z=bt[bt.run==ref]
        display=[m for m in ['Mre11a_det','Mre11a_log','OxPhos','MitoOrg','PyroResp','PyroExec'] if m in metrics]
        fig,axes=plt.subplots(2,3,figsize=(12,7))
        labels={'Mre11a_det':'MreDet','Mre11a_log':'MreLog'}
        colors=['#0072B2','#56B4E9','#009E73','#D55E00','#CC79A7','#E69F00']
        for ax,metric in zip(axes.ravel(),display):
            for sample,color in zip(samples,colors):
                v=z[(z['sample']==sample)&(z.metric==metric)].sort_values('bin')
                ax.plot(v.PT,v.value,'o-',color=color,label=sample.replace('Saline_','Sal').replace('Silica_','Si'),ms=4)
            ax.set(xlabel='DPT',ylabel=labels.get(metric,metric))
        for ax in axes.ravel()[len(display):]:ax.set_visible(False)
        axes.ravel()[0].legend(fontsize=8,frameon=False,ncol=2)
        fig.tight_layout()
        figsave(fig,'reference_trends','Reference=K15_R1, selected before observing trends. Five equal-width DPT bins shared by all samples; bins with fewer than 10 cells have missing values and line breaks. Each line is one sample; no cell-level P values. Sal1-3=Saline_1-3; Si1-3=Silica_1-3. MreDet=Mre11a detection fraction; MreLog=mean log1p(CP10k) Mre11a; OxPhos=oxidative phosphorylation UCell; MitoOrg=mitochondrion organization UCell; PyroResp=pyroptotic inflammatory response UCell; PyroExec=execution machinery transcript UCell. Other UCell definitions are inherited from Stage8.')
(out/'tables/methods.txt').write_text('Input: '+str(src)+'\nScores: '+str(scorepath)+'\n'
 'Trajectory candidates: 1734 Naive/ActCD4 cells; Treg remains in the source and is not forced onto this trajectory.\n'
 'Use Stage4 Harmony PCs 1-30; neighbors=15 and 30. Analyze largest connected component containing >=50 of each state and >=20 cells from every sample. Other components remain documented.\n'
 'Two naive roots: among top quartile Tcf7/Lef1/Ccr7 mean logCP10k, choose closest to their PCA centroid, then closest from a different sample. Neither group label, target expression nor functional scores selects roots.\n'
 'DPT with 10 diffusion components, no branch splitting. Naive-to-activated orientation is a prior to evaluate with state distributions, root and graph stability; it is not established lineage ancestry.\n'
 'Within-sample correlations are descriptive, without cell-level significance tests. Five fixed global DPT bins, minimum10 cells/bin. No pooled-cell regression or imputation.\n'
 'This dataset has a single endpoint. DPT is state ordering, not measured time, and does not establish Mre11a function.\n')
(out/'tables/versions.json').write_text(json.dumps({'scanpy':sc.__version__,'numpy':np.__version__,'pandas':pd.__version__},indent=2))
archive=root/'tmp/GSE314943_stage9_return.tar.gz'
with tarfile.open(str(archive)+'.part','w:gz') as tf:
    for d in ['tables','figures']:tf.add(out/d,arcname=d)
os.replace(str(archive)+'.part',archive)
print('TRAJECTORY_RUNS='+str(len(pseudotimes)),flush=True)
print('GSE314943_STAGE9_PASS',flush=True)
