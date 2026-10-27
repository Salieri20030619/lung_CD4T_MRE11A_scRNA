from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd,numpy as np
W=Path('C:/Users/31937/lung_four_axes_20261017_v15');T=W/'tables';F=W/'figures';F.mkdir(exist_ok=True);plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
a=pd.read_csv(T/'05_intersection/stage_counts.tsv',sep='\t').set_index('scope').loc[['Naive','ActCD4','Treg']]
fig,axs=plt.subplots(1,2,figsize=(12,4.4));cats=['DE','Q_DE','damage_primary','mac_primary'];labels=['DE nominal','DE + quality','+ damage axis','+ macrophage axis'];col=['#8fa8b8','#2c7e9e','#d3994a','#4d947c'];xx=np.arange(3)
for i,(c,label,cl) in enumerate(zip(cats,labels,col)):
 y=a[c].to_numpy();axs[0].bar(xx+(i-1.5)*.2,y,width=.19,label=label,color=cl)
 for x,v in zip(xx+(i-1.5)*.2,y):axs[0].text(x,v+7,str(v),ha='center',fontsize=8)
axs[0].set_xticks(xx,a.index);axs[0].set_ylabel('Gene-state records');axs[0].set_title('Separate axes after DE + quality');axs[0].legend(fontsize=8,frameon=False);axs[0].set_ylim(0,840)
axs[1].bar(xx-.16,a.data_exploratory,width=.3,color='#c7d6dc',label='Pooled intersection (inclusive)');axs[1].bar(xx+.16,a.data_primary,width=.3,color='#197f85',label='Passes group-adjusted screen')
for x,p,n in zip(xx,a.data_primary,a.data_exploratory):axs[1].text(x-.16,n+.6,str(n),ha='center');axs[1].text(x+.16,p+.6,str(p),ha='center')
axs[1].set_xticks(xx,a.index);axs[1].set_ylim(0,39);axs[1].set_title('Same-state intersections before final source review');axs[1].set_ylabel('Gene-state records');axs[1].legend(fontsize=8,frameon=False);fig.tight_layout();fig.savefig(F/'F01_screen_axes.png');fig.savefig(F/'F01_screen_axes.pdf');plt.close(fig)
b=pd.read_csv(T/'05_intersection/DUSP5_sample_plot_values.tsv',sep='\t');raw=b[b.branch=='All'].copy();ids=raw['sample'].to_list();colors=['#3f91b6']*3+['#ce7357']*3
fig,axs=plt.subplots(1,3,figsize=(13,4.6))
for s,c in zip(ids,colors):
 z=b[b['sample']==s].set_index('branch').loc[['All','BG1','T_BG1']];axs[0].plot(range(3),z.DUSP5_CP10k,'o-',color=c,alpha=.8,lw=1,ms=4);axs[0].text(2.05,z.DUSP5_CP10k.iloc[-1],s.replace('Saline_','C').replace('Silica_','S'),fontsize=8,color=c)
axs[0].set_xticks(range(3),['All','BG1','T + BG1']);axs[0].set_xlim(-.2,2.6);axs[0].set_ylabel('Mean DUSP5 CP10k');axs[0].set_title('Treg expression across quality branches')
axs[1].scatter(raw.DUSP5_CP10k,raw.PyroExec_without_DUSP5,c=colors,s=42)
for r,c in zip(raw.itertuples(),colors):axs[1].annotate(r.sample.replace('Saline_','C').replace('Silica_','S'),(r.DUSP5_CP10k,r.PyroExec_without_DUSP5),xytext=(4,4),textcoords='offset points',fontsize=8,color=c)
axs[1].set_xlabel('Treg DUSP5 CP10k');axs[1].set_ylabel('PyroExec (candidate excluded)');axs[1].set_title('rho=-0.943\nPooled P=0.0167; BH=0.377\nAdjusted P=0.0278; BH=0.989')
src=pd.read_csv(T/'06_novelty/DUSP5_link_identity_counts.tsv',sep='\t');v=src[(src.lineage=='Mac')&(src.scope=='AllMac')&(src.gene=='CDH1')].set_index('sample').loc[ids]
axs[2].scatter(v.CP10k,raw.DUSP5_CP10k,c=colors,s=42)
for s,x,y,c in zip(ids,v.CP10k,raw.DUSP5_CP10k,colors):axs[2].annotate(s.replace('Saline_','C').replace('Silica_','S'),(x,y),xytext=(4,4),textcoords='offset points',fontsize=8,color=c)
axs[2].set_xlabel('AllMac CDH1 CP10k');axs[2].set_ylabel('Treg DUSP5 CP10k');axs[2].set_title('rho=1\nPooled P=0.00278; BH=0.668\nAdjusted P=0.0278; BH=1.0')
fig.tight_layout();fig.savefig(F/'F02_DUSP5_raw_evidence.png');fig.savefig(F/'F02_DUSP5_raw_evidence.pdf');plt.close(fig);print('FIGURES_CREATED')