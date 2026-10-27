from pathlib import Path
import numpy as np,pandas as pd,json,shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.patches import Patch,Rectangle
O=Path(__file__).parent;B=O.parent/'lung_client_layout_20261017_v22';A=O.parent/'lung_client_rebuild_20261017_v20';F=O/'figures';plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10.5,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
a=pd.read_csv(A/'plot_data/F17.tsv',sep='\t');edges=pd.read_csv(O/'plot_data/trajectory_backbone_edges.tsv',sep='\t');nodes=pd.read_csv(O/'plot_data/trajectory_backbone_nodes.tsv',sep='\t');assert len(a)==1734
fig,axs=plt.subplots(2,1,figsize=(7,8.5));fig.subplots_adjust(left=.12,right=.76,bottom=.06,top=.91,hspace=.32)
im=axs[0].scatter(a.UMAP1,a.UMAP2,c=a.K15_R1,s=5,cmap='viridis',alpha=.65,rasterized=True);ca=fig.add_axes([.82,.57,.025,.3]);fig.colorbar(im,cax=ca,label='DPT K15_R1')
for st,co in [('Naive','#4c9e9c'),('ActCD4','#dc8451')]:
 z=a[a.subtype==st];axs[1].scatter(z.UMAP1,z.UMAP2,c=co,s=5,alpha=.6,label=st,rasterized=True)
axs[1].legend(loc='upper left',bbox_to_anchor=(1.02,1),frameon=False,fontsize=9)
for ax in axs:
 for r in edges.itertuples():
  ax.plot([r.x1,r.x2],[r.y1,r.y2],color='#263b50',lw=1.7,zorder=4,path_effects=[pe.Stroke(linewidth=3.3,foreground='white'),pe.Normal()]);u=.38;v=.7;ax.annotate('',xy=(r.x1+(r.x2-r.x1)*v,r.y1+(r.y2-r.y1)*v),xytext=(r.x1+(r.x2-r.x1)*u,r.y1+(r.y2-r.y1)*u),arrowprops=dict(arrowstyle='-|>',color='#263b50',lw=1.3,mutation_scale=13),zorder=6)
 ax.scatter(nodes.UMAP1,nodes.UMAP2,s=19,color='#263b50',edgecolors='white',linewidths=.6,zorder=7);ax.set(xlabel='UMAP 1',ylabel='UMAP 2')
axs[0].set_title('DPT with neighborhood backbone',fontsize=12);axs[1].set_title('Backbone across CD4 states',fontsize=12);fig.legend(handles=[Line2D([],[],c='#263b50',marker='o',ms=4,lw=1.5,label='Neighborhood backbone; arrows toward higher DPT')],loc='upper center',bbox_to_anchor=(.5,.99),frameon=False,fontsize=9)
for ext in ['png','pdf']:fig.savefig(F/('F17.'+ext),dpi=300,bbox_inches='tight',pad_inches=.15)
plt.close(fig);shutil.copy2(A/'plot_data/F17.tsv',O/'plot_data/F17.tsv')
a=pd.read_csv(A/'plot_data/F20.tsv',sep='\t');rows=a[['ligand','receptor']].drop_duplicates().sort_values(['ligand','receptor']).reset_index(drop=True);genes=['IL1R2','TTC39C','SKIL'];co={'NicheNet_predicted_target':'#397fa5','direct_receptor':'#d9994e'}
fig,ax=plt.subplots(figsize=(7,7));fig.subplots_adjust(left=.39,right=.97,bottom=.18,top=.90)
for i,r in rows.iterrows():
 for j,g in enumerate(genes):
  z=a[(a.ligand==r.ligand)&(a.receptor==r.receptor)&(a.gene==g)];color='#f0f1f2';lab=''
  if len(z):
   rr=z.iloc[0];kind=rr['route'];color=co[kind];lab=f"Rank {rr.target_rank:.0f}" if pd.notna(rr.target_rank) else 'Decoy'
  ax.add_patch(Rectangle((j-.48,i-.46),.96,.92,facecolor=color,edgecolor='white',linewidth=2));ax.text(j,i,lab,ha='center',va='center',fontsize=9,color='white' if color==co['NicheNet_predicted_target'] else '#3c3022')
ax.set_xlim(-.5,2.5);ax.set_ylim(len(rows)-.5,-.5);ax.set_xticks(range(3),['IL1R2\nTreg','TTC39C\nActCD4','SKIL\nTreg']);ax.set_yticks(range(len(rows)),[r.ligand+'\n('+r.receptor.replace(';',' + ')+')' for r in rows.itertuples()],fontsize=9);ax.tick_params(length=0,pad=8);ax.set_title('Ligand / receptor - candidate evidence',fontsize=12,pad=16)
for sp in ax.spines.values():sp.set_visible(False)
fig.legend(handles=[Patch(color=co['NicheNet_predicted_target'],label='NicheNet target prior (target rank)'),Patch(color=co['direct_receptor'],label='Direct decoy-receptor relationship'),Patch(color='#f0f1f2',label='Not supported in the displayed analysis')],loc='lower center',bbox_to_anchor=(.52,.005),ncol=1,frameon=False,fontsize=9)
for ext in ['png','pdf']:fig.savefig(F/('F20.'+ext),dpi=300,bbox_inches='tight',pad_inches=.15)
plt.close(fig);shutil.copy2(A/'plot_data/F20.tsv',O/'plot_data/F20.tsv');print('BACKBONE_AND_MATRIX_DONE',len(rows),'rows',len(a),'original_route_records')