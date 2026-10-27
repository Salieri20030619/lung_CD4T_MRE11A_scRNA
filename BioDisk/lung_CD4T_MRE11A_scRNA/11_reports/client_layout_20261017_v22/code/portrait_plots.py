from pathlib import Path
import pandas as pd,numpy as np,json,textwrap,shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
O=Path(__file__).parent;B=O.parent/'lung_client_rebuild_20261017_v20';V=O.parent/'lung_client_editorial_20261017_v21';F=O/'figures';F.mkdir(exist_ok=True);T=O/'plot_data';T.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10.5,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
cols={'Saline':'#3278a8','Silica':'#d46746'};st=['Naive','ActCD4','Treg'];sc={'Naive':'#4c9e9c','ActCD4':'#dc8451','Treg':'#8b74aa'};checks=[]
def rd(id):return pd.read_csv(B/'plot_data'/(id+'.tsv'),sep='\t')
def save(fig,id):
 fig.canvas.draw();fb=fig.bbox;bad=[]
 for leg in fig.legends:
  bb=leg.get_window_extent(fig.canvas.get_renderer())
  if not fb.contains(bb.x0,bb.y0) or not fb.contains(bb.x1,bb.y1):bad.append('legend outside figure')
  for ax in fig.axes:
   if bb.overlaps(ax.get_window_extent()):bad.append('figure legend overlaps plot axes')
 checks.append({'figure':id,'legend_flags':bad});assert not bad,(id,bad)
 fig.savefig(F/(id+'.png'),dpi=300,bbox_inches='tight',pad_inches=.12);fig.savefig(F/(id+'.pdf'),bbox_inches='tight',pad_inches=.12);plt.close(fig)
 shutil.copy2(B/'plot_data'/(id+'.tsv'),T/(id+'.tsv'))
def leggroups(fig,y=.985):fig.legend(handles=[Line2D([],[],marker='o',ls='',color=cols[g],label=g) for g in cols],loc='upper center',bbox_to_anchor=(.5,y),ncol=2,frameon=False)
def dots(ax,z,y):
 for i,g in enumerate(cols):
  a=z[z.group==g];ax.scatter(i+np.linspace(-.08,.08,len(a)),a[y],c=cols[g],s=35);ax.hlines(a[y].mean(),i-.15,i+.15,colors='black',lw=1)
 ax.set_xticks([0,1],list(cols))
# Atlas legend gets its own reserved area beneath the plot.
a=rd('F01');types=a.Broad.value_counts();fig,ax=plt.subplots(figsize=(7,8));fig.subplots_adjust(left=.12,right=.97,bottom=.28,top=.94)
pal=list(plt.cm.tab20(np.linspace(0,1,20)))+list(plt.cm.Set2(np.linspace(0,1,3)))
for co,name in zip(pal,types.index):
 z=a[a.Broad==name];ax.scatter(z.UMAP1,z.UMAP2,s=2,c=[co],label=f'{name} ({len(z):,})',alpha=.7,rasterized=True)
ax.set(xlabel='UMAP 1',ylabel='UMAP 2',title='Whole-lung cell atlas | 37,155 cells');fig.legend(*ax.get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.5,.005),ncol=3,fontsize=9,markerscale=4,frameon=False,columnspacing=1.2,handletextpad=.4);save(fig,'F01')
a=rd('F02');genes=['CD3D','TRAC','MS4A1','CD79A','NKG7','LYZ2','CSF1R','S100A8','FCER1A','SFTPC','PECAM1','COL1A1'];fig,ax=plt.subplots(figsize=(7,8));fig.subplots_adjust(left=.19,right=.84,bottom=.18,top=.94);norm=Normalize(0,a.mean_logCP10k.max())
for r in a.itertuples():ax.scatter(genes.index(r.gene),list(types.index).index(r.state),s=3+110*r.fraction,c=[r.mean_logCP10k],cmap='viridis',norm=norm)
ax.set_xticks(range(len(genes)),genes,rotation=65,ha='right',fontsize=9);ax.set_yticks(range(len(types)),types.index,fontsize=9.5);ax.invert_yaxis();ax.set_title('Lineage marker expression',fontsize=12);cax=fig.add_axes([.88,.27,.025,.53]);fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap='viridis'),cax=cax,label='Mean log1p(CP10k)');fig.legend(handles=[plt.scatter([],[],s=3+110*x,c='#888888',label=str(int(x*100))+'%') for x in [.1,.5,1]],loc='lower center',bbox_to_anchor=(.5,.01),ncol=3,frameon=False,title='RNA-positive fraction');save(fig,'F02')
a=rd('F03');fig,axs=plt.subplots(2,1,figsize=(7,8.8));fig.subplots_adjust(left=.12,right=.73,bottom=.07,top=.94,hspace=.32)
for name in st:
 z=a[a.Subtype==name];axs[0].scatter(z.UMAP1,z.UMAP2,s=5,c=sc[name],label=f'{name} ({len(z):,})',alpha=.65,rasterized=True)
for g in cols:
 z=a[a.group==g];axs[1].scatter(z.UMAP1,z.UMAP2,s=5,c=cols[g],label=g,alpha=.6,rasterized=True)
for ax in axs:ax.set(xlabel='UMAP 1',ylabel='UMAP 2');ax.legend(loc='upper left',bbox_to_anchor=(1.02,1),frameon=False,fontsize=9,markerscale=2)
axs[0].set_title('Core CD4 states | 2,311 cells');axs[1].set_title('Disease-group distribution');save(fig,'F03')
a=rd('F04');samples=sorted(a['sample'].unique());wide=a.pivot(index='sample',columns='Subtype',values='fraction').reindex(samples);fig,axs=plt.subplots(2,1,figsize=(7,8));fig.subplots_adjust(left=.14,right=.96,bottom=.08,top=.90,hspace=.62);bt=np.zeros(6)
for name in st:axs[0].bar(range(6),wide[name],bottom=bt,color=sc[name],label=name);bt+=wide[name].to_numpy()
fig.legend(*axs[0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.5,.995),ncol=3,frameon=False);axs[0].set_xticks(range(6),samples,rotation=35,ha='right');axs[0].set_ylabel('Fraction within core CD4');axs[0].set_ylim(0,1)
stat=pd.read_csv(B/'tables/CD4_proportion_tests.tsv',sep='\t')
for i,name in enumerate(st):
 for off,g in [(-.13,'Saline'),(.13,'Silica')]:
  z=a[(a.Subtype==name)&(a.group==g)];axs[1].scatter(i+off+np.linspace(-.035,.035,len(z)),z.fraction,c=cols[g],s=40)
 p=stat[stat.state==name].P.iloc[0];axs[1].text(i,.94,f'P={p:.2g}',ha='center',fontsize=10)
axs[1].set_xticks(range(3),st);axs[1].set_ylabel('Independent-sample fraction');axs[1].set_ylim(0,1);fig.legend(handles=[Line2D([],[],marker='o',ls='',c=cols[g],label=g) for g in cols],loc='lower center',bbox_to_anchor=(.5,.005),ncol=2,frameon=False);save(fig,'F04')
a=rd('F05');fig,axs=plt.subplots(2,1,figsize=(7,8.5));fig.subplots_adjust(left=.14,right=.96,bottom=.08,top=.88,hspace=.4);focus={'Treg':['IL1R2','SKIL','FOXP1'],'ActCD4':['TTC39C','MDFIC','DDIT4','GSN']}
fig.legend(handles=[Patch(color='#bc4d42',label='FDR < 0.05'),Patch(color='#d6a342',label='P < 0.05 only'),Patch(color='#b9c1c7',label='Other')],loc='upper center',bbox_to_anchor=(.5,.995),ncol=3,frameon=False,fontsize=9)
for ax,name in zip(axs,['Treg','ActCD4']):
 z=a[a.scope==name];ax.scatter(z.logFC,z.neglogP,c=np.where(z.FDR<.05,'#bc4d42',np.where(z.PValue<.05,'#d6a342','#b9c1c7')),s=7,alpha=.65,rasterized=True);ax.axhline(-np.log10(.05),c='grey',ls='--',lw=.7);ax.axvline(0,c='grey',lw=.5)
 for i,g in enumerate(focus[name]):
  r=z[z.gene==g].iloc[0];offsets={'TTC39C':(55,-5),'MDFIC':(55,15),'DDIT4':(50,15),'GSN':(-50,16),'IL1R2':(8,8),'SKIL':(20,8),'FOXP1':(-55,10)};ax.annotate(g,(r.logFC,r.neglogP),xytext=offsets[g],textcoords='offset points',fontsize=9,arrowprops={'arrowstyle':'-','lw':.6,'color':'#444444'})
 ax.set(title=name+' | joint disease model (n=10)',xlabel='log2 fold change (Silica / Saline)',ylabel='-log10 raw P')
save(fig,'F05')
a=rd('F06');fig=plt.figure(figsize=(7,6.7));ax=fig.add_axes([.27,.12,.35,.75]);tx=fig.add_axes([.66,.12,.33,.75]);tx.axis('off');yy=np.arange(len(a));ax.set_ylim(len(a)-.5,-.5);tx.set_ylim(len(a)-.5,-.5);ax.set_yticks(yy,[g+' / '+s for g,s in zip(a.gene,a.scope)],fontsize=9);ax.axvline(0,c='grey',lw=.6)
for off,pre,co,lab in [(-.22,'314_disease','#3278a8','Original six'),(.22,'joint_stratified','#d46746','Joint ten')]:
 ax.scatter(a[pre+'_logFC'],yy+off,c=co,s=35,label=lab)
 for i,(_,r) in enumerate(a.iterrows()):tx.text(0,i+off,f'P={r[pre+"_PValue"]:.3g}\nFDR={r[pre+"_FDR"]:.3g}',fontsize=8,va='center',color=co,linespacing=1.05)
ax.set_xlabel('log2 fold change',fontsize=10);tx.set_title('Original P and FDR',fontsize=10);fig.legend(*ax.get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.55,.99),ncol=2,frameon=False);save(fig,'F06')
# GO and KEGG: raw P annotations; FDR remains in full statistical tables.
for id,db in [('F07','GO_BP'),('F08','KEGG')]:
 a=rd(id);fig,axs=plt.subplots(2,1,figsize=(7,9.5));fig.subplots_adjust(left=.49,right=.96,bottom=.06,top=.93,hspace=.38);fig.legend(handles=[Patch(color='#d46746',label='Up'),Patch(color='#3278a8',label='Down')],loc='upper center',bbox_to_anchor=(.65,.995),ncol=2,frameon=False)
 for ax,name in zip(axs,['Treg','ActCD4']):
  z=a[a.scope==name].reset_index(drop=True);labels=[textwrap.fill(s.replace(' - Mus musculus (house mouse)',''),31) for s in z.description];x=-np.log10(z.P.clip(lower=1e-300));ax.scatter(x,range(len(z)),s=20+z.k*3,c=z.direction.map({'up':'#d46746','down':'#3278a8'}));ax.set_yticks(range(len(z)),labels,fontsize=10);ax.invert_yaxis();ax.set_title(name+' | '+db,fontsize=12);ax.set_xlabel('-log10 raw P')
  for i,r in enumerate(z.itertuples()):ax.annotate(f'P={r.P:.2g}',(x.iloc[i],i),xytext=(12,0),textcoords='offset points',fontsize=9.5,va='center')
  ax.set_xlim(0,x.max()*1.55)
 save(fig,id)
# TF AUCell split into three panels, readable six-sample labels.
a=rd('F13');pv=a.pivot(index='tf',columns='column',values='AUC');z=pv.sub(pv.mean(axis=1),axis=0).div(pv.std(axis=1).replace(0,np.nan),axis=0).fillna(0);fig,axs=plt.subplots(3,1,figsize=(7,9));fig.subplots_adjust(left=.16,right=.82,bottom=.07,top=.95,hspace=.65)
for ax,name in zip(axs,st):
 sub=z[[name+' / '+s for s in samples]];im=ax.imshow(sub,cmap='RdBu_r',aspect='auto',vmin=-2.5,vmax=2.5);ax.set_yticks(range(len(sub)),sub.index,fontsize=9);ax.set_xticks(range(6),samples,rotation=30,ha='right',fontsize=9);ax.set_title(name,fontsize=12)
cax=fig.add_axes([.87,.32,.025,.37]);fig.colorbar(im,cax=cax,label='Row z score (display only)');save(fig,'F13')
a=rd('F14');tfs=sorted(a.tf.unique());targets=sorted(a.canonical.unique());fig,ax=plt.subplots(figsize=(7,8));fig.subplots_adjust(left=.02,right=.98,bottom=.04,top=.94);yp={g:i/max(1,len(tfs)-1) for i,g in enumerate(tfs)};yt={g:i/max(1,len(targets)-1) for i,g in enumerate(targets)}
for r in a.itertuples():ax.plot([0,1],[yp[r.tf],yt[r.canonical]],c='#96aebe',alpha=.55,lw=.8)
for g,y in yp.items():ax.scatter(0,y,s=65,c='#3278a8');ax.text(-.04,y,g,ha='right',va='center',fontsize=10)
for g,y in yt.items():ax.scatter(1,y,s=40,c='#d5a55c');ax.text(1.04,y,g,ha='left',va='center',fontsize=10)
ax.set_xlim(-.3,1.4);ax.set_ylim(-.06,1.06);ax.axis('off');ax.set_title('TF                          Predicted targets',fontsize=12);save(fig,'F14')
# Existing score plots are already portrait; keep plotting input and rerender font-safe panels.
E=pd.read_csv(B.parent/'lung_append_analysis_20261017_v18/tables/programme_group_effects.tsv',sep='\t')
for id,mods in [('F15',['PyroExec','Inflam']),('F16',['ApoExec','MitoISR'])]:
 a=rd(id);fig,axs=plt.subplots(3,2,figsize=(7,8.7));fig.subplots_adjust(left=.13,right=.97,bottom=.07,top=.95,hspace=.52,wspace=.40)
 for i,name in enumerate(st):
  for j,mod in enumerate(mods):
   ax=axs[i,j];zz=a[(a.scope==name)&(a.program==mod)];dots(ax,zz,'score');r=E[(E.cohort.astype(str)=='314')&(E.variant=='Raw')&(E.scope==name)&(E.program==mod)].iloc[0];ax.set_title(f'{name} | {mod}\nDelta={r.effect:+.3f}\nP={r.P:.3g}; FDR={r.FDR:.3g}',fontsize=9.5);ax.set_ylabel('Program score',fontsize=9.5)
 save(fig,id)
# DPT: stack maps and place categorical legend outside.
a=rd('F17');fig,axs=plt.subplots(2,1,figsize=(7,8.2));fig.subplots_adjust(left=.13,right=.78,bottom=.07,top=.94,hspace=.35);im=axs[0].scatter(a.UMAP1,a.UMAP2,c=a.K15_R1,s=5,cmap='viridis',rasterized=True);cax=fig.add_axes([.83,.56,.025,.32]);fig.colorbar(im,cax=cax,label='DPT K15_R1')
for name in st[:2]:
 z=a[a.subtype==name];axs[1].scatter(z.UMAP1,z.UMAP2,s=5,c=sc[name],label=name,alpha=.65,rasterized=True)
axs[1].legend(loc='upper left',bbox_to_anchor=(1.01,1),frameon=False,fontsize=9)
for ax in axs:ax.set(xlabel='UMAP 1',ylabel='UMAP 2')
axs[0].set_title('Naive / ActCD4 state ordering');axs[1].set_title('Same cells, CD4 state labels');save(fig,'F17')
a=rd('F18');tests=pd.read_csv(B/'sources/11_reports/customer_sequence_20261015_v2/03A_state_pseudotime/program_trajectory_tests.tsv',sep='\t');fig,axs=plt.subplots(2,2,figsize=(7,7.5));fig.subplots_adjust(left=.12,right=.97,bottom=.09,top=.88,hspace=.4,wspace=.36);leggroups(fig)
for i,name in enumerate(st[:2]):
 for j,mod in enumerate(['PyroExec','ApoExec']):
  ax=axs[i,j];z=a[(a.state==name)&(a.module==mod)]
  for sample,b in z.groupby('sample'):ax.plot(b.mean_pt,b.score,'o-',c=cols['Saline' if sample.startswith('Saline') else 'Silica'],alpha=.65,lw=1,ms=3)
  r=tests[(tests.state==name)&(tests.run=='K15_R1')&(tests.feature==mod)].iloc[0];ax.set_title(f'{name} | {mod}\nSlope={r.mean_slope:+.3f}\nP={r.P:.3g}; FDR={r.BH:.3g}',fontsize=9);ax.set_xlabel('Within-state DPT',fontsize=9);ax.set_ylabel('Matched-control score',fontsize=9)
save(fig,'F18')
a=rd('F19');pv=a.pivot_table(index='test_ligand',columns='scope',values='pearson').reindex(columns=st);fig,ax=plt.subplots(figsize=(7,7));fig.subplots_adjust(left=.24,right=.80,bottom=.1,top=.94);cm=plt.cm.YlGnBu.copy();cm.set_bad('#ededed');im=ax.imshow(pv,cmap=cm,aspect='auto');ax.set_xticks(range(3),st);ax.set_yticks(range(len(pv)),pv.index);ax.set_title('NicheNet ligand activity',fontsize=12);cax=fig.add_axes([.85,.20,.025,.58]);fig.colorbar(im,cax=cax,label='Prior-prediction Pearson');save(fig,'F19')
a=rd('F20');genes=['IL1R2','TTC39C','SKIL'];ligs=sorted(a.ligand.unique());fig,ax=plt.subplots(figsize=(7,6.7));fig.subplots_adjust(left=.02,right=.98,bottom=.07,top=.94);yp={g:i for i,g in enumerate(genes)};yl={g:i/max(1,len(ligs)-1)*2 for i,g in enumerate(ligs)}
for r in a.itertuples():ax.plot([0,1],[yl[r.ligand],yp[r.gene]],c='#96aebe',alpha=.6,lw=1)
for g,y in yl.items():ax.scatter(0,y,c='#3278a8',s=65);ax.text(-.04,y,g,ha='right',va='center',fontsize=10)
for g,y in yp.items():ax.scatter(1,y,c='#d46746',s=90);ax.text(1.04,y,g,ha='left',va='center',fontsize=10)
ax.set_xlim(-.3,1.4);ax.set_ylim(-.2,2.2);ax.axis('off');ax.set_title('Ligand                    Candidate',fontsize=12);save(fig,'F20')
a=rd('F21');fig,axs=plt.subplots(3,1,figsize=(7,8));fig.subplots_adjust(left=.16,right=.95,bottom=.07,top=.95,hspace=.48)
for ax,(g,name) in zip(axs,[('IL1R2','Treg'),('TTC39C','ActCD4'),('SKIL','Treg')]):
 z=a[(a.gene==g)&(a.scope==name)];dots(ax,z,'positive_fraction');ax.set_title(g+' / '+name,fontsize=12);ax.set_ylabel('RNA-positive fraction');ax.set_ylim(0,max(.35,z.positive_fraction.max()*1.2))
save(fig,'F21')
a=rd('F22');assoc=pd.read_csv(B.parent/'lung_append_analysis_20261017_v18/tables/all_associations.tsv',sep='\t');fig,axs=plt.subplots(3,1,figsize=(7,9));fig.subplots_adjust(left=.16,right=.95,bottom=.08,top=.87,hspace=.62)
fig.legend(handles=[Line2D([],[],marker='o',ls='',c='black',label='GSE314943'),Line2D([],[],marker='^',ls='',c='black',label='GSE183682'),Line2D([],[],marker='o',ls='',c=cols['Saline'],label='Saline'),Line2D([],[],marker='o',ls='',c=cols['Silica'],label='Silica')],loc='upper center',bbox_to_anchor=(.5,.995),ncol=2,frameon=False,fontsize=10)
for ax,(g,name,mod) in zip(axs,[('IL1R2','Treg','Inhib'),('TTC39C','ActCD4','ApoExec'),('SKIL','Treg','TAct')]):
 z=a[a.scope==name].drop_duplicates(['dataset','sample','scope'])
 for ds,marker in [('GSE314943','o'),('GSE183682','^')]:
  zz=z[z.dataset==ds];ax.scatter(zz['G_'+g],zz['M_'+mod],c=zz.group.map(cols),marker=marker,s=40)
 r=assoc[(assoc.cohort=='Joint')&(assoc.variant=='Raw')&(assoc.scope==name)&(assoc.model=='JointStratum')&(assoc.x=='G_'+g)&(assoc.y=='M_'+mod)].iloc[0];ax.set_title(f'{g} / {name}\nAdjusted rho={r.rho:.3f}; P={r.P:.3g}; FDR={r.FDR:.3g}',fontsize=10);ax.set_xlabel(g+' log1p(CP10k)');ax.set_ylabel(mod+' score')
save(fig,'F22')
(O/'python_layout_QA.json').write_text(json.dumps(checks,indent=2),encoding='utf-8');print('PORTRAIT_PLOTS',len(checks))