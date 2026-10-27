from pathlib import Path
import pandas as pd,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/death_priority_20261018_v2');F=O/'figures'
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':200})
g=pd.read_csv(O/'tables/joint_GSEA_all.tsv',sep='\t');patterns=['apoptosis - mus','pyroptotic inflammatory response','canonical glycolysis','autophagy of mitochondrion','regulation of T cell apoptotic process','oxidative phosphorylation -']
z=g[g.description.str.contains('|'.join(patterns),case=False,na=False)].copy();z['label']=z.description.str.replace(' - Mus musculus (house mouse)','',regex=False);labels=list(dict.fromkeys(z.label));fig,ax=plt.subplots(figsize=(12,6))
for r in z.itertuples():ax.scatter(['Core','Naive','ActCD4','Treg'].index(r.scope),labels.index(r.label),s=30+50*min(-np.log10(max(r.pval,1e-30)),5),c=[r.NES],cmap='RdBu_r',vmin=-2.5,vmax=2.5,edgecolor='black' if r.padj<.05 else '#bbbbbb')
ax.set(xticks=range(4),xticklabels=['Core','Naive','ActCD4','Treg'],yticks=range(len(labels)),yticklabels=labels);ax.invert_yaxis();ax.grid(alpha=.15)
fig.colorbar(ScalarMappable(norm=Normalize(-2.5,2.5),cmap='RdBu_r'),ax=ax,label='NES',shrink=.7)
hs=[ax.scatter([],[],s=30+50*min(-np.log10(p),5),c='#bbbbbb',label=f'P={p:g}') for p in [.05,.01,.001]]
ax.legend(handles=hs,loc='upper left',bbox_to_anchor=(0,-.13),ncol=3,frameon=False,title='Raw P (black edge: FDR < 0.05)')
plt.tight_layout();plt.savefig(F/'N06.png',bbox_inches='tight');plt.close()
m=pd.read_csv(O/'tables/GSE183682_cluster_markers.tsv',sep='\t');genes=['C1QA','ADGRE1','MERTK','MARCO','MSR1','PPARG','FABP4','CAR4','SIGLECF','LY6C2','CCR2','FLT3','CD209A','CD3E','CD79A'];fig,ax=plt.subplots(figsize=(12,4.5))
for i,cl in enumerate([2,4,5,9]):
 t=m[m.cluster==cl]
 for j,gene in enumerate(genes):ax.scatter(j,i,s=15+180*t[gene+'_pct'].mean(),c=[t[gene].mean()],cmap='viridis',vmin=0,vmax=3,edgecolors='grey',lw=.25)
ax.set(xticks=range(len(genes)),xticklabels=genes,yticks=range(4),yticklabels=['2 Mac/AMlike','4 Mac/Mono-derived','5 Mono','9 DC/MonoAPC']);plt.setp(ax.get_xticklabels(),rotation=45,ha='right');ax.invert_yaxis()
fig.colorbar(ScalarMappable(norm=Normalize(0,3),cmap='viridis'),ax=ax,label='Mean log1p(CP10k)',shrink=.7)
hs=[ax.scatter([],[],s=15+180*p,c='#888888',label=f'{int(p*100)}%') for p in [.1,.5,1]]
ax.legend(handles=hs,loc='upper left',bbox_to_anchor=(0,-.35),ncol=3,frameon=False,title='Positive cells (sample mean)')
plt.tight_layout();plt.savefig(F/'N10.png',bbox_inches='tight');plt.close()
cover=pd.read_csv(O/'tables/sender_coverage.tsv',sep='\t');s=cover[['dataset','sample','group','day']].drop_duplicates().sort_values(['dataset','sample']);order=list(s['sample']);fig,ax=plt.subplots(figsize=(12,4))
for j,source in enumerate(['AllMac','MacStrict']):
 d=cover[cover.source==source].set_index('sample').loc[order];ax.bar(np.arange(10)+(j-.5)*.34,d.cells,width=.34,label=source,color=['#287c8e','#bc7457'][j])
labels=[f'{r.sample}\n{r.group}\nD{r.day}' for r in s.itertuples()];ax.set(xticks=range(10),xticklabels=labels,ylabel='Macrophage-enriched cells',title='Ten biological samples; cell counts are not inferential n');ax.tick_params(axis='x',labelsize=8);ax.legend();plt.tight_layout();plt.savefig(F/'N11.png',bbox_inches='tight');plt.close()
print('LEGENDS_UPDATED')
