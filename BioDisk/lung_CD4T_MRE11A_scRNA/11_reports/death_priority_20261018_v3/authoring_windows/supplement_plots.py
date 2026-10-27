from pathlib import Path
import pandas as pd,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).parent;T=P/'tables';F=P/'supplement_figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':210})
def read(n):return pd.read_csv(T/n,sep='\t')
def heat(mat,pv,title,label,name):
 fig,ax=plt.subplots(figsize=(7.2,8));v=mat.to_numpy();lim=np.nanmax(abs(v));im=ax.imshow(v,cmap='RdBu_r',vmin=-lim,vmax=lim,aspect='auto')
 ax.set(xticks=range(len(mat.columns)),xticklabels=mat.columns,yticks=range(len(mat)),yticklabels=mat.index,title=title)
 for i in range(len(mat)):
  for j in range(len(mat.columns)):
   if np.isfinite(v[i,j]):ax.text(j,i,f'{v[i,j]:.2f}'+('*' if pv.iloc[i,j]<.05 else ''),ha='center',va='center',fontsize=11)
 fig.colorbar(im,ax=ax,label=label,shrink=.6);fig.tight_layout();fig.savefig(F/(name+'.png'),bbox_inches='tight');plt.close(fig)
d=read('external_module_effects.tsv');z=d[d.contrast.str.startswith('GSE121883')].copy();z['context']=z.contrast.str.replace('GSE121883_','',regex=False)
mat=z.pivot(index='program',columns='context',values='effect');pv=z.pivot(index='program',columns='context',values='P').reindex_like(mat)
heat(mat,pv,'IL2RA KO - WT | independent mouse Tregs','Mean logCPM difference (* raw P < 0.05)','E01')
d=read('mouse_HIF1A_perturbation_GSEA.tsv');d=d[d.database=='Modules'];h=read('external_module_GSEA.tsv');h=h[h.contrast.str.startswith('GSE324038')]
z=pd.concat([d,h]);z['context']=z.contrast.map({'GSE253090_Hypoxia':'Mouse\nhypoxia','GSE253090_Normoxia':'Mouse\nnormoxia','GSE324038_HumanPanT_Hypoxia':'Human\nhypoxia'})
mat=z.pivot(index='pathway',columns='context',values='NES');pv=z.pivot(index='pathway',columns='context',values='pval').reindex_like(mat)
heat(mat,pv,'HIF1A KO - WT | targeted program enrichment','NES (* raw P < 0.05; separate families)','E02')
s=read('human_sample_scores.tsv');aa=read('human_candidate_associations.tsv');fig,axs=plt.subplots(3,1,figsize=(7,10.2))
for ax,(g,sub,mod) in zip(axs,[('IL2RA','Treg','Glycolysis'),('HIF1A','CD4-CCR7','PyroExec'),('ARL6IP1','CD4-CCR7','MitoROS')]):
 z=s[(s.series=='G227')&(s.subtype==sub)&(s.cells>=10)]
 for group,col in [('Control','#287c8e'),('CWP','#bf6445')]:
  q=z[z.group==group];ax.scatter(q['G_'+g],q[mod],s=60,c=col,label=group,edgecolors='white')
 r=aa[(aa.cohort=='G227')&(aa.gene==g)&(aa.program==mod)&(aa.min_cells==10)&(aa.model=='Source')].iloc[0]
 ax.set(xlabel=g+' mean log1p(CP10k)',ylabel=mod+' score',title=f'{sub}, n={len(z)} | rho={r.rho:.2f}, P={r.P:.3g}, FDR={r.FDR:.3g}');ax.legend(loc='best')
fig.tight_layout();fig.savefig(F/'E03.png',bbox_inches='tight');plt.close(fig)
if (T/'network_LOO_candidate_edges.tsv').exists():
 d=read('network_LOO_candidate_edges.tsv');base=read('network_baseline_seed17_edges.tsv.gz');d=d.merge(base,on=['arm','tf','target'],how='left');d['baseline_seed17_present']=d.baseline_seed17_present.fillna(False).astype(bool);d['self_signature_member']=d.tf.eq(d.target);d.to_csv(T/'network_LOO_candidate_edges_annotated.tsv',sep='\t',index=False,encoding='utf-8');d=d[(d.arm=='All')&d.tf.ne(d.target)&d.stable_2of3.astype(str).str.lower().eq('true')&d.target.isin(['HIF1A','IL2RA','ARL6IP1'])].sort_values(['target','folds_present','seed_support'],ascending=[True,False,False])
 fig,ax=plt.subplots(figsize=(7.2,max(7,.25*len(d))));labels=d.tf+' -> '+d.target
 ax.barh(np.arange(len(d)),d.folds_present,color=np.where(d.baseline_seed17_present,'#287c8e','#a8b4b8'));ax.set(yticks=np.arange(len(d)),yticklabels=labels,xticks=range(11),xlabel='Occurrences in 10 leave-one-mouse-out networks',title='Original All 2/3-seed motif edges\nTeal: present in full seed17; grey: absent in seed17');ax.invert_yaxis();fig.tight_layout();fig.savefig(F/'E04.png',bbox_inches='tight');plt.close(fig)
print('SUPPLEMENT_FIGURES_COMPLETE')
