from pathlib import Path
import pandas as pd,numpy as np,shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory
O=Path(__file__).parent;H=O.parent;F=O/'figures';T=O/'plot_data';plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
colors={'314':'#347dad','Joint':'#d9773c'}
def panel(ax,z,title,unit):
 labs=[]
 for i,sc in enumerate(['Naive','ActCD4','Treg']):
  for j,co in enumerate(['314','Joint']):
   yy=i*3+j;r=z[z.scope.eq(sc)&z.cohort.eq(co)].iloc[0];ax.scatter(r.effect,yy,c=colors[co],s=48,zorder=3);ax.text(1.04,yy,f'P={r.P:.3g}; FDR={r.FDR:.3g}',transform=blended_transform_factory(ax.transAxes,ax.transData),va='center',fontsize=11);labs.append((yy,sc+' / '+('3+3' if co=='314' else '5+5')))
 ax.set_yticks([y for y,l in labs],[l for y,l in labs],fontsize=10);ax.set_ylim(8,-1);ax.axvline(0,c='#999',lw=.8);ax.set_title(title,fontsize=13,pad=12);ax.set_xlabel(unit,fontsize=10);ax.margins(x=.2)
def save(fig,n):
 for e in ['png','pdf']:fig.savefig(F/(n+'.'+e),dpi=300,bbox_inches='tight',pad_inches=.18)
 plt.close(fig)
e=pd.read_csv(O/'tables/group_scores_both_models.tsv',sep='\t')
for name,mods in [('F15',['PyroExec','Inflam']),('F16',['ApoExec','MitoISR'])]:
 fig,axs=plt.subplots(2,1,figsize=(7.5,8.5));fig.subplots_adjust(left=.23,right=.58,bottom=.08,top=.85,hspace=.5)
 for ax,mod in zip(axs,mods):panel(ax,e[e.program.eq(mod)],mod,'Stratum-adjusted difference' if False else 'Disease effect (score units)')
 fig.legend(handles=[Line2D([],[],marker='o',ls='',color=colors['314'],label='Primary cohort: 3+3'),Line2D([],[],marker='o',ls='',color=colors['Joint'],label='Joint stratified: 5+5')],loc='upper center',bbox_to_anchor=(.53,.995),ncol=1,frameon=False,fontsize=11)
 save(fig,name);e[e.program.isin(mods)&e.scope.isin(['Naive','ActCD4','Treg'])].to_csv(T/(name+'.tsv'),sep='\t',index=False)
w=pd.read_csv(T/'F04.tsv',sep='\t');w=w.sort_values(['dataset','day','group','sample'],ascending=[False,True,True,True]);fig=plt.figure(figsize=(7.5,9));ax=fig.add_axes([.13,.59,.80,.32]);bottom=np.zeros(len(w));statecolors={'Naive':'#4c9e9c','ActCD4':'#dc8451','Treg':'#8b74aa'}
for sc,c in statecolors.items():ax.bar(np.arange(len(w)),w[sc+'_fraction'],bottom=bottom,color=c,label=sc);bottom+=w[sc+'_fraction'].to_numpy()
labs=[]
for r in w.itertuples():labs.append(('314 / '+r.sample.replace('Saline','C').replace('Silica','S')) if r.dataset=='GSE314943' else f'183 / D{r.day} / '+('C' if r.group=='Saline' else 'S'))
ax.set_xticks(range(len(w)),labs,rotation=45,ha='right',fontsize=9);ax.set_ylabel('Fraction of Core CD4');ax.set_title('CD4 composition in all 10 samples',fontsize=13);ax.legend(loc='upper center',bbox_to_anchor=(.5,1.23),ncol=3,frameon=False,fontsize=11)
ax=fig.add_axes([.23,.09,.35,.30]);pr=pd.read_csv(O/'tables/proportion_joint_supplement.tsv',sep='\t');panel(ax,pr,'Sample-level proportion effects','Difference in proportion');save(fig,'F04')
print('COMPARISON_FIGURES_DONE')