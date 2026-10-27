from pathlib import Path
import sys
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path(sys.argv[1]);m=pd.read_csv(out/'input/cells.tsv.gz',sep='\t');s=pd.read_csv(out/'tables/cell_scores.tsv',sep='\t');mods=list(s.columns[1:]);m=m.merge(s,on='cell',validate='one_to_one')
rows=[]
for scope in ['Core','Expanded']:
 z=m[m.core] if scope=='Core' else m
 for subtype in ['All']+sorted(z.subtype.unique().tolist()):
  zz=z if subtype=='All' else z[z.subtype==subtype]
  for ss,v in zz.groupby('sample'):
   for mod in mods:
    rows.append(dict(scope=scope,subtype=subtype,sample=ss,group=v.group.iloc[0],module=mod,cells=len(v),score_mean=v[mod].mean(),score_median=v[mod].median(),mean_UMI=v.library_UMI.mean(),MRE11_pct=100*(v.MRE11_UMI>0).mean(),MRE11_CPM=1e6*v.MRE11_UMI.sum()/v.library_UMI.sum()))
a=pd.DataFrame(rows);a.to_csv(out/'tables/sample_scores.tsv',sep='\t',index=False)
rows=[]
for (scope,mod),v in a[a.subtype=='All'].groupby(['scope','module']):
 for rule in [20,1]:
  z=v[v.cells>=rule];e=z[z.group=='Exposed'];c=z[z.group=='Silicosis'];effect=c.score_mean.mean()-e.score_mean.mean()
  rows.append(dict(scope=scope,module=mod,min_cells=rule,n_exposed=len(e),n_silicosis=len(c),exposed_mean=e.score_mean.mean(),silicosis_mean=c.score_mean.mean(),difference=effect,direction='Higher' if effect>0 else 'Lower' if effect<0 else 'Equal',p_status='descriptive_small_sample'))
r=pd.DataFrame(rows);r.to_csv(out/'tables/group_directions.tsv',sep='\t',index=False)
# Composition-standardized score, common weights from eligible Core samples.
base=m[m.core];sizes=base.groupby('sample').size();eligible=sizes[sizes>=20].index
ct=pd.crosstab(base['sample'],base.subtype);frac=ct.div(ct.sum(axis=1),axis=0);weights=frac.loc[eligible].mean();weights.rename('weight').to_csv(out/'tables/standardization_weights.tsv',sep='\t')
std=[]
for (ss,mod),v in a[(a.scope=='Core')&(a.subtype!='All')].groupby(['sample','module']):
 vv=v.set_index('subtype').reindex(weights.index);complete=vv.score_mean.notna().all()
 std.append(dict(sample=ss,group=v.group.iloc[0],module=mod,cells=sizes[ss],eligible20=ss in eligible,score_std=(vv.score_mean*weights).sum() if complete else np.nan,status='complete' if complete else 'missing_subtype'))
std=pd.DataFrame(std);std.to_csv(out/'tables/core_standardized_scores.tsv',sep='\t',index=False)
samples=['GSM5324866','GSM5324867','GSM5324868','GSM5324869','GSM5324870'];short=['Exp1','Exp2','Sil1','Sil2','Sil3']
fig,axs=plt.subplots(2,len(mods),figsize=(3*len(mods),6),squeeze=False)
for j,scope in enumerate(['Core','Expanded']):
 for i,mod in enumerate(mods):
  ax=axs[j,i];v=a[(a.scope==scope)&(a.subtype=='All')&(a.module==mod)].set_index('sample').reindex(samples)
  for k,(_,rr) in enumerate(v.iterrows()):
   color='#31688e' if k<2 else '#c55635';ax.scatter(k,rr.score_mean,facecolors=color if rr.cells>=20 else 'none',edgecolors=color,s=42);ax.annotate('n='+str(int(rr.cells)),(k,rr.score_mean),xytext=(0,7),textcoords='offset points',ha='center',fontsize=7)
  ax.set(title=scope+' / '+mod,xticks=range(5),xticklabels=short,ylabel='UCell');ax.margins(y=.3);ax.spines[['top','right']].set_visible(False)
fig.tight_layout()
for ext in ['png','pdf']:fig.savefig(out/'figures'/('human_scores.'+ext),dpi=350,bbox_inches='tight')
plt.close(fig)
fig,ax=plt.subplots(figsize=(5,3.5));q=r[r.min_cells==20].pivot(index='module',columns='scope',values='difference').reindex(mods);vmax=max(np.abs(q.to_numpy()).max(),1e-8);im=ax.imshow(q.to_numpy(),cmap='RdBu_r',vmin=-vmax,vmax=vmax,aspect='auto');ax.set(xticks=range(2),xticklabels=q.columns,yticks=range(len(q)),yticklabels=q.index)
for i in range(len(q)):
 for j in range(2):ax.text(j,i,f'{q.iloc[i,j]:+.4f}',ha='center',va='center',color='black' if abs(q.iloc[i,j])<vmax*.6 else 'white')
fig.colorbar(im,ax=ax,label='Sil - Exp');fig.tight_layout()
for ext in ['png','pdf']:fig.savefig(out/'figures'/('scope_effects.'+ext),dpi=350,bbox_inches='tight')
plt.close(fig)
(out/'captions.txt').write_text('Human scores: x-axis independent samples Exp1/Exp2=exposed comparators GSM5324866/67, Sil1/Sil2/Sil3=silicosis GSM5324868/69/70. Y-axis mean per-cell UCell score within sample/scope. Filled points have >=20 cells; hollow points have <20 cells; n denotes cell number. Core=F08 clusters 0/1/12, Expanded adds mixed CD4/CD8 clusters 2/3. OxPhos=oxidative phosphorylation; OxStress=response to oxidative stress; Inflam=inflammasome-mediated signaling pathway; PyroResp=pyroptotic inflammatory response. Scope effects: rows modules, columns cell scope, color and printed value are equally weighted sample mean Silicosis minus Exposed for >=20-cell samples. Positive indicates higher score in silicosis. Descriptive differences, not significance. Standardized scores use common Core subtype weights from eligible samples; absent subtypes yield missing values. MemCD4=memory-like CD4, NaiveCD4=naive-like CD4, Treg=regulatory T, Mix2/Mix3=mixed CD4/CD8 clusters.\n')
(out/'methods.txt').write_text('Four GO programs fixed before human scoring, same GO definitions as mouse cohorts; human GOALL symbols derived from org.Hs.eg.db, not capitalization or claimed identical cross-species gene membership. Actual annotated/present/detected genes and annotation versions exported. MRE11 and MRE11A excluded. UCell 2.14.0 maxRank=1500, full-feature raw counts, signature restricted to present features exactly as main-cohort approach. Minimum 3 present features; coverage and detection must be interpreted. Full measured feature background retained. Core >=20 cells primary descriptive direction; Expanded and all-sample summaries sensitivity. No new differential-expression, group-significance or MRE11-correlation tests for this five-subject supplement. All subject points retained; groups weighted equally. Composition-standardization uses fixed Core subtype proportions. Scores describe transcriptional programs; no cross-species absolute-score comparison.\n')
print('GSE174725_SUPPC2_ANALYSIS_PASS',flush=True)
