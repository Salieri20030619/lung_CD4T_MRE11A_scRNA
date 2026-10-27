"""Freeze marker-supported core CD4 and export sample-level results."""
import sys, os, json, tarfile, itertools
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
source=root/'02_analysis/GSE314943/stage4_20260905_215942_3XPis4/objects/TNK_harmony.h5ad'
out=Path(sys.argv[1]); out.mkdir(parents=True,exist_ok=True)
for d in ['tables','figures','objects']:(out/d).mkdir(exist_ok=True)
def save(t,name,index=False):t.to_csv(out/'tables'/name,sep='\t',index=index)
def figure(fig,name,caption):
    for ext in ['png','pdf']:fig.savefig(out/'figures'/f'{name}.{ext}',dpi=350,bbox_inches='tight')
    plt.close(fig)
    (out/'figures'/f'{name}_caption.txt').write_text(caption+'\n')
print('SOURCE='+str(source),flush=True)
a=sc.read_h5ad(source)
assert a.n_obs==6303 and 'counts' in a.layers
mapping={'0':'Naive','2':'ActCD4','4':'Treg'}
a.obs['CoreCD4']=a.obs.H06.astype(str).isin(mapping)
save(a.obs,'TNK_membership.tsv.gz',True)
c=a[a.obs.CoreCD4].copy()
c.obs['Subtype']=pd.Categorical(c.obs.H06.astype(str).map(mapping),categories=list(mapping.values()))
assert c.n_obs==2311
samples=['Saline_1','Saline_2','Saline_3','Silica_1','Silica_2','Silica_3']
expected=[304,424,164,553,390,476]
assert c.obs['sample'].value_counts().reindex(samples).tolist()==expected
c.X=c.layers['counts'].copy();sc.pp.normalize_total(c,target_sum=1e4);sc.pp.log1p(c)
# Preserve the existing Harmony coordinates; do not fit another embedding.
c.obsm['X_umap']=c.obsm['X_umap_Harmony'].copy()
colors={'Naive':'#0072B2','ActCD4':'#D55E00','Treg':'#009E73'}
fig,ax=plt.subplots(figsize=(7,6))
for typ in colors:
    z=c.obs.Subtype==typ;xy=c.obsm['X_umap'][z]
    ax.scatter(xy[:,0],xy[:,1],s=5,c=colors[typ],label=typ,rasterized=True)
ax.legend(frameon=False,markerscale=3);ax.set(xlabel='UMAP1',ylabel='UMAP2');ax.set_xticks([]);ax.set_yticks([])
figure(fig,'CD4_UMAP','Core CD4: H0=Naive (naive-like CD4 T cells); H2=ActCD4 (activated CD4 T cells); H4=Treg (regulatory T cells). Coordinates are inherited from the Stage4 TNK Harmony embedding. Selection is based on cluster-level lineage markers, not individual Cd4 or Mre11a detection.')
n=pd.crosstab(c.obs['sample'],c.obs.Subtype).reindex(samples).fillna(0).astype(int)
frac=n.div(n.sum(axis=1),axis=0)*100
save(n,'subtype_counts.tsv',True);save(frac,'subtype_percent.tsv',True)
fig,ax=plt.subplots(figsize=(8,4));bottom=np.zeros(6)
for typ in colors:
    ax.bar(range(6),frac[typ],bottom=bottom,color=colors[typ],label=typ);bottom+=frac[typ].values
ax.set_xticks(range(6),['Sal1','Sal2','Sal3','Si1','Si2','Si3']);ax.set_ylabel('CD4 (%)');ax.legend(frameon=False,bbox_to_anchor=(1,1))
figure(fig,'CD4_comp','Subtype proportions within the selected core CD4 population. Naive=naive-like CD4; ActCD4=activated CD4; Treg=regulatory T cells. Sal1-3=Saline_1-3; Si1-3=Silica_1-3. Each bar is one sample. This denominator is core CD4, not all lung cells.')
# All pseudobulk values are sums of the unmodified counts layer, by sample.
genes=pd.DataFrame({'gene_id':c.var_names,'symbol':c.var.symbol.values});save(genes,'genes.tsv')
idx=np.flatnonzero(c.var_names=='ENSMUSG00000031928');assert len(idx)==1
j=int(idx[0]);print('TARGET_SYMBOL='+str(c.var.iloc[j]['symbol']),flush=True)
rows=[]
for typ in ['Core',*colors]:
    x=c if typ=='Core' else c[c.obs.Subtype==typ]
    matrix=[]
    for sample in samples:
        z=x[x.obs['sample']==sample]
        counts=np.asarray(z.layers['counts'].sum(axis=0)).ravel()
        assert np.isfinite(counts).all() and np.allclose(counts,np.round(counts))
        matrix.append(np.round(counts).astype(np.int64))
        v=np.asarray(z.layers['counts'][:,[j]].sum(axis=1)).ravel()
        rows.append(dict(set=typ,sample=sample,group='Saline' if sample.startswith('Saline') else 'Silica',
            cells=z.n_obs,library_umi=int(counts.sum()),Mre_pos=int((v>0).sum()),Mre_umi=int(v.sum()),
            Mre_pct=100*float((v>0).mean()) if z.n_obs else np.nan,
            Mre_CPM=float(counts[j]/counts.sum()*1e6) if counts.sum() else np.nan))
    save(pd.DataFrame(np.stack(matrix,axis=1),index=c.var_names,columns=samples),f'PB_{typ}_counts.tsv.gz',True)
meta=pd.DataFrame(rows);save(meta,'PB_sample_metadata.tsv')
core=meta[meta['set']=='Core'].set_index('sample').loc[samples]
assert core.Mre_pos.sum()==121
fig,axes=plt.subplots(1,2,figsize=(8,4))
for ax,metric,label in zip(axes,['Mre_pct','Mre_CPM'],['Mre11a (%)','Mre11a CPM']):
    for k,group in enumerate(['Saline','Silica']):
        y=core.loc[core.group==group,metric].values
        ax.scatter(k+np.array([-.08,0,.08]),y,s=50,color=['#0072B2','#D55E00'][k])
        ax.plot([k-.15,k+.15],[np.mean(y)]*2,c='black')
    ax.set_xticks([0,1],['Saline','Silica']);ax.set_ylabel(label)
figure(fig,'Mre11a_sample','Each dot is one sample; black lines show sample means. Left: fraction of core CD4 cells with at least one Mre11a UMI. Right: Mre11a counts per million total pseudobulk UMIs (CPM; library-size scaling only, not TMM). Mre11a is identified by ENSMUSG00000031928. No imputation is used.')
# Exact sample-label permutation: all 20 assignments for n=3 versus n=3.
# Descriptive exploratory tests; do not treat cells as biological replicates.
def permutation(v):
    v=np.asarray(v,dtype=float);obs=v[3:].mean()-v[:3].mean();null=[]
    for comb in itertools.combinations(range(6),3):
        mask=np.zeros(6,dtype=bool);mask[list(comb)]=True
        null.append(v[mask].mean()-v[~mask].mean())
    return obs,float(np.mean(np.abs(null)>=abs(obs)-1e-12))
results=[]
for metric,values in [('Mre_pct',core.Mre_pct),('Mre_CPM',core.Mre_CPM),*[(k+'_pct',frac[k]) for k in colors]]:
    effect,p=permutation(values)
    results.append(dict(metric=metric,Saline_mean=float(np.mean(values[:3])),Silica_mean=float(np.mean(values[3:])),
                        Silica_minus_Saline=effect,p_exact=p))
t=pd.DataFrame(results)
p=t.p_exact.to_numpy();order=np.argsort(p);q=np.empty(len(p));q[order]=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1].clip(0,1)
t['FDR_BH']=q;save(t,'sample_exploratory_tests.tsv')
save(c.obs,'CD4_metadata.tsv.gz',True)
c.write_h5ad(out/'objects/CD4_core.h5ad',compression='gzip')
(out/'tables/methods.txt').write_text('Source: '+str(source)+'\nOutput: '+str(out)+'\nCore clusters: H0,H2,H4. Cluster-level lineage selection; no per-cell Cd4 or Mre11a positivity filter.\nCore is a conservative subset, not a census of every CD4 state. H3,H10,H13 and other unresolved groups remain in the source TNK object.\nPseudobulk counts: unmodified RNA counts summed by sample; Core and each subtype exported separately. Metadata cells must be used to decide subtype DEG eligibility. Six sample identities preserved; biological-unit/pooling metadata must be checked before inferential DEG.\nExact two-sided sample-label permutation tests use 20 allocations; minimum attainable nonzero P=0.1. Five exploratory tests share a BH adjustment. Results are sample-level descriptive evidence.\nProportions use core CD4 as denominator. This cohort has one endpoint, so no observed time-course trend is estimated.\n')
archive=root/'tmp/GSE314943_stage5_return.tar.gz';archive.parent.mkdir(exist_ok=True)
with tarfile.open(str(archive)+'.part','w:gz') as tf:
    for d in ['tables','figures']:tf.add(out/d,arcname=d)
os.replace(str(archive)+'.part',archive)
print('RETURN_PACKAGE='+str(archive),flush=True)
print('GSE314943_STAGE5_PASS',flush=True)
