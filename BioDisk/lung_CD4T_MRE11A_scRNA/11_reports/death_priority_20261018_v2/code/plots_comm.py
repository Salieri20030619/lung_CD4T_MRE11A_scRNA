import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad,json,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/death_priority_20261018_v2');F=O/'figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
manifest=[]
def read(n):return pd.read_csv(O/'tables'/n,sep='\t')
def save(n,title,source):
 plt.savefig(F/(n+'.png'),bbox_inches='tight',facecolor='white');plt.close();manifest.append(dict(id=n,title=title,source=source,file=str(F/(n+'.png'))))
a=ad.read_h5ad(O/'objects/GSE183682_whole_clustered.h5ad');a.obs['Broad']=a.obs.cluster.astype(str).map(json.loads((O/'inputs/GSE183682_annotation.json').read_text()));z=a.obsm['X_umap']
fig,axs=plt.subplots(1,2,figsize=(12,5))
for ax,key in zip(axs,['Broad','sample']):
 values=list(a.obs[key].unique());cmap=plt.cm.tab20.colors
 for i,v in enumerate(values):
  take=a.obs[key]==v;ax.scatter(z[take,0],z[take,1],s=1,c=[cmap[i%20]],label=v,rasterized=True,alpha=.65)
 ax.set(xticks=[],yticks=[],xlabel='UMAP1',ylabel='UMAP2',title='GSE183682 | '+key);ax.legend(markerscale=5,fontsize=7,ncol=2,loc='best')
plt.tight_layout();save('N09','第二队列全肺图谱及样本覆盖','GSE183682_whole_clustered.h5ad; GSE183682_annotation.json')
m=read('GSE183682_cluster_markers.tsv');genes=['C1QA','ADGRE1','MERTK','MARCO','MSR1','PPARG','FABP4','CAR4','SIGLECF','LY6C2','CCR2','FLT3','CD209A','CD3E','CD79A']
clusters=[2,4,5,9];fig,ax=plt.subplots(figsize=(11,3.3))
for i,cl in enumerate(clusters):
 t=m[m.cluster==cl]
 for j,g in enumerate(genes):ax.scatter(j,i,s=15+180*t[g+'_pct'].mean(),c=[t[g].mean()],cmap='viridis',vmin=0,vmax=3,edgecolors='grey',lw=.25)
ax.set(xticks=range(len(genes)),xticklabels=genes,yticks=range(4),yticklabels=['2 Mac/AMlike','4 Mac/Mono-derived','5 Mono','9 DC/MonoAPC'],title='Marker mean across four biological samples; size = positive fraction');plt.setp(ax.get_xticklabels(),rotation=45,ha='right');ax.invert_yaxis();plt.tight_layout();save('N10','巨噬细胞范围由标志表达与样本覆盖确定','GSE183682_cluster_markers.tsv')
cover=read('sender_coverage.tsv');samples=cover[['dataset','sample','group','day']].drop_duplicates().sort_values(['dataset','sample']);samplelist=list(samples['sample']);fig,ax=plt.subplots(figsize=(11,3.5))
for j,source in enumerate(['AllMac','MacStrict']):
 d=cover[cover.source==source].set_index('sample').loc[samplelist];ax.bar(np.arange(10)+(j-.5)*.34,d.cells,width=.34,label=source,color=['#287c8e','#bc7457'][j])
labels=[f'{r.dataset.replace("GSE","")}\n{r.group}\nD{r.day}' for r in samples.itertuples()]
ax.set(xticks=np.arange(10),xticklabels=labels,ylabel='Macrophage-enriched cells',title='Ten biological samples; cell counts are not inferential n');ax.legend();plt.tight_layout();save('N11','完整5+5发送端覆盖与严格巨噬范围','sender_coverage.tsv')
lr=read('joint_LR_disease_effects.tsv');d=lr[lr.variant=='AllMac'].copy();top=d.sort_values(['P','effect'],ascending=[True,False]).program.drop_duplicates().head(18);states=['Naive','ActCD4','Treg']
mat=d.pivot(index='program',columns='scope',values='effect').reindex(index=top,columns=states);pv=d.pivot(index='program',columns='scope',values='P').reindex_like(mat)
fig,ax=plt.subplots(figsize=(9,7));v=mat.to_numpy();lim=np.nanmax(abs(v));im=ax.imshow(v,cmap='RdBu_r',vmin=-lim,vmax=lim,aspect='auto');ax.set(xticks=range(3),xticklabels=states,yticks=range(len(top)),yticklabels=top)
for i in range(len(top)):
 for j in range(3):
  if np.isfinite(pv.iloc[i,j]):ax.text(j,i,f'P={pv.iloc[i,j]:.3f}',ha='center',va='center',fontsize=8)
fig.colorbar(im,ax=ax,label='Silica - Saline, log1p LR proxy',shrink=.6);plt.tight_layout();save('N12','联合配体受体表达代理的疾病效应','joint_LR_disease_effects.tsv')
nn=read('joint_NicheNet_activity.tsv');fig,axs=plt.subplots(2,3,figsize=(12,7))
for i,direction in enumerate(['up','down']):
 for j,sc in enumerate(states):
  d=nn[(nn.source=='AllMac')&(nn.scope==sc)&(nn.response=='joint_disease_DE')&(nn.direction==direction)].sort_values('activity_rank').head(6).iloc[::-1];ax=axs[i,j];ax.barh(d.test_ligand,d.aupr_corrected,color='#287c8e' if i==0 else '#bf7656');ax.axvline(0,color='grey',lw=.6);ax.set(title=sc+' | '+direction,xlabel='Corrected AUPR (prediction)')
plt.tight_layout();save('N13','联合疾病响应的NicheNet配体排序','joint_NicheNet_activity.tsv')
whole=read('whole_lung_source_expression.tsv.gz');fig,axs=plt.subplots(1,2,figsize=(11,4))
for ax,ds in zip(axs,['GSE314943','GSE183682']):
 d=whole[(whole.dataset==ds)&(whole.gene=='IL2')].copy();d['logCP']=np.log1p(d.CP10k);s=d.groupby('Broad').logCP.mean().sort_values(ascending=False).head(10).iloc[::-1]
 ax.barh(s.index,s.values,color=['#bf7656' if str(v)=='Mac' else '#287c8e' for v in s.index]);ax.set(title=ds+' | IL2 sources',xlabel='Mean sample log1p(CP10k)')
plt.tight_layout();save('N14','IL2的全肺竞争来源','whole_lung_source_expression.tsv.gz')
loo=read('candidate_DE_leave_one_donor.tsv');de=read('candidate_existing_DE.tsv');pairs=[('HIF1A','Naive'),('IL2RA','Treg'),('ARL6IP1','Naive')]
fig,ax=plt.subplots(figsize=(10,3.5))
for i,(gene,scope) in enumerate(pairs):
 d=loo[(loo.gene==gene)&(loo.scope==scope)];r=de[(de.gene==gene)&(de.scope==scope)&(de.model=='joint_stratified')].iloc[0]
 ax.plot([d.logFC.min(),d.logFC.max()],[i,i],lw=3,c='#a7bfc6');ax.scatter(d.logFC,np.full(len(d),i)+np.linspace(-.09,.09,len(d)),s=20,c='#287c8e');ax.scatter(r.logFC,i,c='#c76141',marker='D',s=65,zorder=5)
 ax.text(max(d.logFC.max(),r.logFC)+.08,i,f'{(np.sign(d.logFC)==np.sign(r.logFC)).sum()}/10 same direction',va='center',fontsize=9)
ax.axvline(0,color='grey',lw=.7);ax.set(yticks=range(3),yticklabels=[g+' | '+sc for g,sc in pairs],xlabel='Disease log2FC; range across leave-one-donor fits (not CI)',xlim=(-.7,1.8));plt.tight_layout();save('N15','候选疾病差异的逐样本剔除稳定性','candidate_DE_leave_one_donor.tsv; candidate_existing_DE.tsv')
(F/'manifest_comm.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8');print('PLOTS_COMM',len(manifest))
