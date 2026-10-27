from pathlib import Path
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path('/workspace/scratch/38d90e8ea00a');out=ROOT/'deliverables/scenic3_review';colors={'All':'#176B9B','BG1':'#C75327'}
e=pd.read_csv(ROOT/'scenic3_results/tables/mouse_TF_identity_adjustment.tsv',sep='\t');e=e[e.model=='DepthTreg']
fig,axes=plt.subplots(1,2,figsize=(9,4),sharey=True)
for ax,tf in zip(axes,['Ikzf2','Stat1']):
 z=e[e.tf==tf].sort_values(['filter','seed']).reset_index(drop=True)
 for j,r in z.iterrows():
  ax.hlines(j,r.loo_min,r.loo_max,color=colors[r['filter']],lw=2)
  ax.scatter(r.coef,j,c=colors[r['filter']],s=38,zorder=3)
 ax.axvline(0,color='#777777',lw=.8,ls='--');ax.set_yticks(range(len(z)),[f"{r['filter']}-S{r.seed}" for _,r in z.iterrows()]);ax.set(xlabel='AUC diff.',title=tf.upper());ax.invert_yaxis();ax.spines[['top','right']].set_visible(False)
fig.tight_layout()
for ext in ['png','pdf']:fig.savefig(out/f'TF_Identity.{ext}',dpi=300)
plt.close(fig)
m=pd.read_csv(ROOT/'scenic3_results/tables/human_cluster_marker_expression.tsv',sep='\t');genes=['TYROBP','C1QA','MARCO','FABP4','FCN1','S100A8','CD1C','CLEC9A','CD3D','EPCAM','SFTPC'];labs=['C0-AM','C1-Mono','C2-AM','C3-MoMac','C4-Mono','C5-Mono','C6-DC2','C7-MoMac','C8-EpiMix','C9-TMix','C10-Mixed','C11-AM','C12-DC1','C13-AM','C14-Mono']
pct=m.pivot(index='r4',columns='gene',values='pct').reindex(index=range(15),columns=genes)
val=np.log1p(m.pivot(index='r4',columns='gene',values='mean_CP10k').reindex(index=range(15),columns=genes));val=(val-val.min())/(val.max()-val.min())
fig,ax=plt.subplots(figsize=(9,7));xx,yy=np.meshgrid(range(len(genes)),range(15));s=ax.scatter(xx.ravel(),yy.ravel(),s=pct.to_numpy().ravel()*1.3,c=val.to_numpy().ravel(),cmap='viridis',vmin=0,vmax=1)
ax.set_xticks(range(len(genes)),genes,rotation=60,ha='right');ax.set_yticks(range(15),labs);ax.invert_yaxis();ax.set_title('Myel ID');ax.spines[['top','right']].set_visible(False);fig.colorbar(s,ax=ax,label='Expr scale',shrink=.6)
for v in [10,50,100]:ax.scatter([],[],s=v*1.3,c='#777777',label=str(v)+'%')
ax.legend(title='Detected',loc='upper left',bbox_to_anchor=(1.16,1),frameon=False);fig.tight_layout()
for ext in ['png','pdf']:fig.savefig(out/f'Myel_Identity.{ext}',dpi=300)
plt.close(fig)
(out/'figure_captions.txt').write_text('TF_Identity: mouse ActCD4 disease coefficient after adjustment for depth and TReg score excluding IKZF2. Dots: full-data coefficients, line segments: leave-one-subject-out coefficient ranges, NOT confidence intervals. n=3 Saline+3 Silica; residual df=2. S17/S29/S43: random seeds, not independent subjects. All: original cells; BG1: frozen background-filtered cells; AUC diff.: adjusted regulon AUCell difference.\nMyel_Identity: marker-reviewed human myeloid candidate clusters. Dot area proportional to detection percentage; color is per-gene min-max-scaled log1p(mean CP10k) across clusters. C: cluster at resolution .4. AM: alveolar/C1Q macrophage-like; Mono: monocyte-like; MoMac: monocyte/macrophage transitional; DC1/DC2: dendritic-like; EpiMix: epithelial/myeloid mixed signals; TMix: T/myeloid mixed signals; Mixed: endothelial/epithelial/myeloid mixed signals. These are working annotations; C2/C11 are dominated by individual subjects. Myel ID: myeloid identity; Expr scale: scaled expression; Detected: cells with positive count.\n',encoding='utf-8')
print('FIGURES_PASS')
