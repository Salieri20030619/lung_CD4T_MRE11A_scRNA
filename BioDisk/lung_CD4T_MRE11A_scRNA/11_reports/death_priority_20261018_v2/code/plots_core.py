import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np,pandas as pd,json
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/death_priority_20261018_v2';F=O/'figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
def read(n):return pd.read_csv(O/'tables'/n,sep='\t')
manifest=[]
def save(name,title,source):
 plt.savefig(F/(name+'.png'),bbox_inches='tight',facecolor='white');plt.close();manifest.append(dict(id=name,title=title,source=source,file=str(F/(name+'.png'))))
states=['Naive','ActCD4','Treg'];primary=[('HIF1A','Naive'),('IL2RA','Treg'),('ARL6IP1','Naive')]
colors={'Saline':'#287c8e','Silica':'#c75e42'};markers={'GSE314943':'o','GSE183682':'s'}
s=read('sample_scores.tsv');a=read('all_associations.tsv');e=read('program_effects.tsv');de=read('candidate_existing_DE.tsv')
fig,axs=plt.subplots(1,3,figsize=(11,3.8))
for ax,(g,sc) in zip(axs,primary):
 z=s[(s.scope==sc)&(s.variant=='Raw')&s.dataset.ne('GSE274372')]
 for ds in z.dataset.unique():
  for gr in ['Saline','Silica']:
   t=z[(z.dataset==ds)&(z.group==gr)];xx=np.full(len(t),['Saline','Silica'].index(gr))+np.linspace(-.15,.15,len(t));ax.scatter(xx,t['G_'+g],c=colors[gr],marker=markers[ds],s=50,label=ds if gr=='Saline' else None,edgecolors='white',linewidth=.5)
 dd=de[(de.model=='joint_stratified')&(de.scope==sc)&(de.gene==g)].iloc[0]
 ax.set(xticks=[0,1],xticklabels=['Saline','Silica'],title=f'{g} | {sc}\nlog2FC={dd.logFC:.2f}; P={dd.PValue:.3g}\nFDR={dd.FDR:.3g}',ylabel='log1p(pseudobulk CP10k)');ax.legend(fontsize=8)
plt.tight_layout();save('N01','主候选的联合疾病表达','candidate_existing_DE.tsv; sample_scores.tsv')
mods=['TAct','TReg','Inhib','Naive','Cycle','ApoExec','PyroExec','MitoROS','MitoISR','Mitophagy','Autophagy','OXPHOS','Glycolysis']
z=e[(e.cohort=='Joint')&(e.variant=='Raw')&e.scope.isin(states)];mat=z.pivot(index='scope',columns='program',values='effect').reindex(index=states,columns=mods);p=z.pivot(index='scope',columns='program',values='P').reindex_like(mat)
fig,ax=plt.subplots(figsize=(12,3.2));im=ax.imshow(mat,cmap='RdBu_r',vmin=-np.nanmax(abs(mat.to_numpy())),vmax=np.nanmax(abs(mat.to_numpy())),aspect='auto');ax.set(xticks=range(len(mods)),xticklabels=mods,yticks=range(3),yticklabels=states)
plt.setp(ax.get_xticklabels(),rotation=45,ha='right')
for i in range(3):
 for j in range(len(mods)):ax.text(j,i,f'{mat.iloc[i,j]:.2f}'+('*' if p.iloc[i,j]<.05 else ''),ha='center',va='center',fontsize=8)
fig.colorbar(im,ax=ax,label='Silica - Saline (score units)',shrink=.7);plt.tight_layout();save('N02','联合样本功能程序疾病效应','program_effects.tsv')
mat=[];ps=[]
for g,sc in primary:
 z=a[(a.cohort=='Joint')&(a.variant=='Raw')&(a.model=='JointStratum')&(a.gene==g)&(a.scope==sc)].set_index('program');mat.append(z.rho.reindex(mods).to_numpy());ps.append(z.P.reindex(mods).to_numpy())
mat=np.array(mat);ps=np.array(ps);fig,ax=plt.subplots(figsize=(12,3.3));im=ax.imshow(mat,vmin=-1,vmax=1,cmap='RdBu_r',aspect='auto');ax.set(xticks=range(len(mods)),xticklabels=mods,yticks=range(3),yticklabels=[g+' | '+sc for g,sc in primary]);plt.setp(ax.get_xticklabels(),rotation=45,ha='right')
for i in range(3):
 for j in range(len(mods)):ax.text(j,i,f'{mat[i,j]:.2f}'+('*' if ps[i,j]<.05 else ''),ha='center',va='center',fontsize=8)
fig.colorbar(im,ax=ax,label='Stratum-adjusted rank correlation');plt.tight_layout();save('N03','主候选在原提名亚群内的功能关联','all_associations.tsv')
chosen=[('HIF1A','Naive','PyroExec'),('ARL6IP1','Naive','MitoROS'),('IL2RA','Treg','OXPHOS'),('IL2RA','Treg','Glycolysis'),('IL2RA','Treg','ApoExec')]
variants=[('Raw','JointStratum','Primary'),('BG1','JointStratum','BG1'),('Raw','JointStratumCondition','Condition adjusted'),('BG1','JointStratumCondition','BG1 + adjusted')]
fig,axs=plt.subplots(1,5,figsize=(13,4),sharey=True)
for ax,(g,sc,mod) in zip(axs,chosen):
 for i,(var,model,label) in enumerate(variants):
  r=a[(a.cohort=='Joint')&(a.variant==var)&(a.model==model)&(a.gene==g)&(a.scope==sc)&(a.program==mod)].iloc[0]
  ax.plot([0,r.rho],[i,i],c='#cad4d8');ax.scatter(r.rho,i,s=50,c='#bf563c' if r.P<.05 else '#35778b');ax.text(-.96,i+.24,f'P={r.P:.3g}; q={r.FDR:.3g}',fontsize=7.5)
 ax.axvline(0,color='grey',lw=.7);ax.set(xlim=(-1,1),title=f'{g} | {sc}\n{mod}',xlabel='Partial rank rho',yticks=range(4),yticklabels=[v[2] for v in variants]);ax.invert_yaxis()
plt.tight_layout();save('N04','候选功能关联的条件调整与背景敏感性','all_associations.tsv')
tf=read('TF_disease_effects.tsv');focal=['HIF1A','FOXP1','ATF4','TRP53','STAT1','STAT3','NFKB1','NFKB2','RELA','RELB','JUN','JUNB','FOS','FOSB','NFATC1','NFATC2','FOXO1','FOXO3','TFEB']
z=tf[(tf.variant=='Raw')&tf.scope.isin(states)];focal=[x for x in focal if x in set(z.program)]
mat=z.pivot(index='scope',columns='program',values='effect').reindex(index=states,columns=focal);pv=z.pivot(index='scope',columns='program',values='P').reindex_like(mat)
fig,ax=plt.subplots(figsize=(12,3));lim=np.nanmax(abs(mat.to_numpy()));im=ax.imshow(mat,cmap='RdBu_r',vmin=-lim,vmax=lim,aspect='auto');ax.set(xticks=range(len(focal)),xticklabels=focal,yticks=range(3),yticklabels=states);plt.setp(ax.get_xticklabels(),rotation=45,ha='right')
for i in range(3):
 for j in range(len(focal)):
  if pv.iloc[i,j]<.05:ax.text(j,i,'*',ha='center',va='center',fontsize=14)
fig.colorbar(im,ax=ax,label='Silica - Saline (AUC units)');plt.tight_layout();save('N05','去重叠TF活性的联合疾病效应','TF_disease_effects.tsv')
gs=read('joint_GSEA_all.tsv');ss=read('Core_composition_adjusted_GSEA.tsv')
patterns=['apoptosis - mus','pyroptotic inflammatory response','canonical glycolysis','autophagy of mitochondrion','regulation of T cell apoptotic process','oxidative phosphorylation -']
selected=gs[gs.description.str.contains('|'.join(patterns),case=False,na=False)].copy();selected['label']=selected.description.str.replace(' - Mus musculus (house mouse)','',regex=False);labels=list(dict.fromkeys(selected.label))
fig,ax=plt.subplots(figsize=(10,5))
for r in selected.itertuples():
 xx=['Core']+states;ax.scatter(xx.index(r.scope),labels.index(r.label),s=30+50*min(-np.log10(max(r.pval,1e-30)),5),c=[r.NES],cmap='RdBu_r',vmin=-2.5,vmax=2.5,edgecolor='black' if r.padj<.05 else '#cccccc')
ax.set(xticks=range(4),xticklabels=['Core']+states,yticks=range(len(labels)),yticklabels=labels);ax.invert_yaxis();ax.grid(alpha=.15);plt.tight_layout();save('N06','联合GSEA的死亡与代谢相关通路','joint_GSEA_all.tsv')
z=gs[(gs.scope=='Core')&gs.description.str.contains('pyroptotic inflammatory response|Apoptosis - Mus|canonical glycolysis|lymphocyte apoptotic process',case=False,na=False)&~gs.description.str.contains('regulation')].merge(ss,on=['pathway','database'],suffixes=('_primary','_adjusted'))
fig,axs=plt.subplots(1,len(z),figsize=(12,3.5))
for ax,r in zip(np.atleast_1d(axs),z.itertuples()):
 ax.bar([0,1],[r.NES_primary,r.NES_adjusted],color=['#287c8e','#ca7155']);ax.axhline(0,color='grey',lw=.7)
 for j,tag in enumerate(['primary','adjusted']):ax.text(j,0,f'P={getattr(r,"pval_"+tag):.3g}\nq={getattr(r,"padj_"+tag):.3g}',ha='center',va='bottom',fontsize=8)
 ax.set(xticks=[0,1],xticklabels=['Primary','Adjusted'],title=r.description_primary.replace(' - Mus musculus (house mouse)','').replace(' ','\n'),ylabel='NES')
plt.tight_layout();save('N07','Core富集信号的细胞组成敏感性','joint_GSEA_all.tsv; Core_composition_adjusted_GSEA.tsv')
# Actual running ES curves, leading hits and source ranking.
for sc in states:
 rank=read('joint_rank_'+sc+'.tsv');zz=gs[(gs.scope==sc)&(gs.database=='KEGG')].sort_values(['pval','pathway']).head(5)
 # RDS memberships are exported below by an independent tiny R command; wait for file.
 pfile=O/'tables/pathway_symbol_members.tsv.gz'
 if not pfile.exists():continue
 members=pd.read_csv(pfile,sep='\t');stat=rank['rank'].to_numpy();names=rank.gene.to_numpy();n=len(stat)
 fig,axs=plt.subplots(3,1,figsize=(10,7),gridspec_kw={'height_ratios':[3,1,1.3]},sharex=True);cols=plt.cm.tab10.colors
 for i,r in enumerate(zz.itertuples()):
  hit=np.isin(names,members.loc[(members.database=='KEGG')&(members.pathway==r.pathway),'gene']);pos=np.flatnonzero(hit);inc=np.full(n,-1/(n-hit.sum()));inc[hit]=abs(stat[hit])/sum(abs(stat[hit]));es=np.cumsum(inc);observed=es.max() if es.max()>-es.min() else es.min();assert abs(observed-r.ES)<1e-6
  axs[0].plot(np.arange(n),es,c=cols[i],label=f'{r.pathway} | NES={r.NES:.2f}; P={r.pval:.3g}; q={r.padj:.3g}')
  axs[1].vlines(pos,i-.3,i+.3,color=cols[i],lw=.5)
 axs[0].legend(fontsize=8,loc='best');axs[0].set(ylabel='Running enrichment score',title=sc+' | KEGG top 5 by raw P');axs[0].axhline(0,color='grey',lw=.6);axs[1].set(ylabel='Hits',yticks=[]);axs[2].plot(np.arange(n),stat,c='#475c65');axs[2].set(xlabel='Full gene rank',ylabel='Signed sqrt(F)');plt.tight_layout();save('N08_'+sc,'联合'+sc+'的KEGG富集曲线','joint_GSEA_all.tsv; joint_rank_'+sc+'.tsv')
(F/'manifest_core.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8');print('PLOTS',len(manifest))
