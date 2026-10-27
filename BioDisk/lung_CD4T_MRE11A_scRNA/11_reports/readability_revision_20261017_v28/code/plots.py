from pathlib import Path
import pandas as pd,numpy as np,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
O=Path(__file__).parent;H=O.parent;B=H/'lung_mainline_revision_20261017_v27';F=O/'figures';T=O/'plot_data'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
a=pd.read_csv(B/'supplement/plot_data/F19.tsv',sep='\t');a.to_csv(T/'F19.tsv',sep='\t',index=False)
x=a.pivot(index='test_ligand',columns='scope',values='aupr_corrected').reindex(columns=['Naive','ActCD4','Treg']);r=a.pivot(index='test_ligand',columns='scope',values='activity_rank').reindex(index=x.index,columns=x.columns)
fig,ax=plt.subplots(figsize=(7,8));fig.subplots_adjust(left=.22,right=.79,bottom=.11,top=.92);cmap=plt.get_cmap('YlGnBu').copy();cmap.set_bad('#eeeeee');im=ax.imshow(x.to_numpy(),aspect='auto',cmap=cmap,vmin=0)
ax.set_xticks(range(3),x.columns);ax.set_yticks(range(len(x)),x.index);ax.tick_params(length=0,pad=8);ax.set_title('Ligand prediction | corrected AUPR',fontsize=13,pad=15)
for i in range(len(x)):
 for j in range(3):
  if pd.notna(x.iloc[i,j]):ax.text(j,i,f'#{int(r.iloc[i,j])}',ha='center',va='center',fontsize=11,color='white' if x.iloc[i,j]>.035 else '#152b3a')
ca=fig.add_axes([.85,.26,.025,.5]);fig.colorbar(im,cax=ca,label='Corrected AUPR');fig.text(.5,.025,'Selection, ranks and colors: corrected AUPR\nRanks are defined within each state; grey = no result',ha='center',fontsize=10)
for e in ['png','pdf']:fig.savefig(F/('F19.'+e),dpi=300,bbox_inches='tight',pad_inches=.15)
plt.close(fig)
a=pd.read_csv(B/'supplement/plot_data/F24.tsv',sep='\t');a.to_csv(T/'F24.tsv',sep='\t',index=False);assoc=pd.read_csv(H/'lung_append_analysis_20261017_v18/tables/all_associations.tsv',sep='\t');cols={'Saline':'#2e7cb4','Silica':'#dd6543'}
fig,axs=plt.subplots(3,1,figsize=(7.5,9.5));fig.subplots_adjust(left=.17,right=.96,bottom=.12,top=.79,hspace=.95)
fig.legend(handles=[Line2D([],[],marker='o',ls='',c='black',label='GSE314943'),Line2D([],[],marker='^',ls='',c='black',label='GSE183682'),Line2D([],[],marker='o',ls='',c=cols['Saline'],label='Saline'),Line2D([],[],marker='o',ls='',c=cols['Silica'],label='Silica')],loc='upper center',bbox_to_anchor=(.5,.995),ncol=2,frameon=False,fontsize=11)
fig.text(.5,.915,'Coordinates: original values (not residuals)',ha='center',va='top',fontsize=11)
fig.text(.5,.02,'M1: cohort/time stratum adjusted\nM2: M1 + disease status adjusted',ha='center',fontsize=11)
source=[]
for ax,(g,sc,mod) in zip(axs,[('IL1R2','Treg','Inhib'),('TTC39C','ActCD4','ApoExec'),('SKIL','Treg','TAct')]):
 z=a[a.scope==sc].drop_duplicates(['dataset','sample','scope'])
 for ds,mark in [('GSE314943','o'),('GSE183682','^')]:
  zz=z[z.dataset==ds];ax.scatter(zz['G_'+g],zz['M_'+mod],c=zz.group.map(cols),marker=mark,s=42)
 lines=[g+' / '+sc]
 for model,label in [('JointStratum','M1'),('JointStratumCondition','M2')]:
  rr=assoc[(assoc.cohort=='Joint')&(assoc.variant=='Raw')&(assoc.scope==sc)&(assoc.model==model)&(assoc.x=='G_'+g)&(assoc.y=='M_'+mod)];assert len(rr)==1;source.append(rr);v=rr.iloc[0];lines.append(f'{label}: rho={v.rho:.3f}; P={v.P:.3g}; FDR={v.FDR:.3g}')
 ax.set_title('\n'.join(lines),fontsize=12,pad=10);ax.set_xlabel(g+' log1p(CP10k)',fontsize=11);ax.set_ylabel(mod+' score',fontsize=11)
for e in ['png','pdf']:fig.savefig(F/('F24.'+e),dpi=300,bbox_inches='tight',pad_inches=.15)
pd.concat(source).to_csv(O/'evidence/two_model_plot_statistics.tsv',sep='\t',index=False);plt.close(fig)
print('PLOTS_DONE')