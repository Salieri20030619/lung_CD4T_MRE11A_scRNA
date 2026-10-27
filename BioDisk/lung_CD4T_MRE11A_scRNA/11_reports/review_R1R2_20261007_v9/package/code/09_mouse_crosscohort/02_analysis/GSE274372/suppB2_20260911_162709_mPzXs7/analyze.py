import sys,itertools,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path(sys.argv[1]); rng=np.random.default_rng(42)
def wr(d,n):d.to_csv(out/'tables'/n,sep='\t',index=False)
def bh(p):
 p=np.asarray(p,float);q=np.full(len(p),np.nan);k=np.flatnonzero(np.isfinite(p));o=k[np.argsort(p[k])];q[o]=np.minimum(1,np.minimum.accumulate((p[o]*len(o)/np.arange(1,len(o)+1))[::-1])[::-1]);return q
m=pd.read_csv(out/'input/cells.tsv',sep='\t',dtype={'R08':str})
s=pd.read_csv(out/'tables/cell_scores.tsv',sep='\t');mods=list(s.columns[1:]);m=m.merge(s,on='cell',validate='one_to_one')
rows=[]
for scope in ['TLike','Author']:
 z=m[m.R08.isin(['0','1'])] if scope=='TLike' else m
 for sample,v in z.groupby('sample',sort=True):
  for mod in mods:rows.append(dict(scope=scope,sample=sample,day=int(v.day.iloc[0]),cells=len(v),module=mod,score=v[mod].mean(),Mre_CPM=v.Mre_UMI.sum()/v.library_UMI.sum()*1e6,Mre_pct=(v.Mre_UMI>0).mean()*100,mean_UMI=v.library_UMI.mean()))
a=pd.DataFrame(rows);wr(a,'sample_scores.tsv')
# Independent samples: permutation omnibus; exact 3-vs-3 or 2-vs-3 baseline contrasts.
tr=[]; cr=[]
def between(y,g):return sum(np.sum(g==k)*(y[g==k].mean()-y.mean())**2 for k in np.unique(g))
for (scope,mod),v0 in a.groupby(['scope','module']):
 for rule in ['All','Min10']:
  v=v0 if rule=='All' else v0[v0.cells>=10]
  y=v.score.to_numpy();g=v.day.to_numpy();n=len(v)
  ob=between(y,g);null=np.array([between(y,rng.permutation(g)) for _ in range(9999)])
  p=(1+np.sum(null>=ob-1e-15))/10000 if np.ptp(y)>0 else np.nan
  tr.append(dict(scope=scope,module=mod,rule=rule,contrast='TimeAll',n=n,effect=np.nan,p=p,method='sample_permutation_9999'))
  for day in [7,28,56]:
   y0=y[g==0];yd=y[g==day];yy=np.r_[y0,yd];obs=yd.mean()-y0.mean();vals=[]
   for idx in itertools.combinations(range(len(yy)),len(yd)):
    take=np.zeros(len(yy),bool);take[list(idx)]=True;vals.append(yy[take].mean()-yy[~take].mean())
   p=np.mean(np.abs(vals)>=abs(obs)-1e-15) if np.ptp(yy)>0 else np.nan
   tr.append(dict(scope=scope,module=mod,rule=rule,contrast=f'D{day}-D0',n=len(yy),effect=obs,p=p,method='exact_sample_permutation'))
  # Rank associations; time-adjusted rank regression with categorical day.
  Z=np.column_stack([np.ones(n)]+[(g==d).astype(float) for d in [7,28,56]])
  for metric in ['Mre_CPM','Mre_pct']:
   xx=stats.rankdata(v[metric]);yy=stats.rankdata(y)
   for adjust in ['None','Day']:
    if adjust=='Day':
     xr=xx-Z@np.linalg.lstsq(Z,xx,rcond=None)[0];yr=yy-Z@np.linalg.lstsq(Z,yy,rcond=None)[0];df=n-np.linalg.matrix_rank(Z)-1
    else:xr=xx-xx.mean();yr=yy-yy.mean();df=n-2
    valid=np.linalg.norm(xr)>1e-10 and np.linalg.norm(yr)>1e-10 and df>0
    r=float(np.clip(np.dot(xr,yr)/np.linalg.norm(xr)/np.linalg.norm(yr),-1,1)) if valid else np.nan
    p=2*stats.t.sf(abs(r)*np.sqrt(df/max(1e-15,1-r*r)),df) if valid else np.nan
    cr.append(dict(scope=scope,module=mod,rule=rule,metric=metric,adjust=adjust,n=n,df=df,rho=r,p=p,method='rank_correlation_t_approx',status='exploratory' if valid else 'insufficient_variation'))
t=pd.DataFrame(tr);c=pd.DataFrame(cr)
# Each scope and inclusion rule is a separate prespecified family.
for d in [t,c]:
 d['FDR']=np.nan
 for _,ix in d.groupby(['scope','rule']).groups.items():d.loc[ix,'FDR']=bh(d.loc[ix,'p'])
wr(t,'time_tests.tsv');wr(c,'Mre_associations.tsv')
# Whole-lung annotations reused as supplied; denominator explicit.
w=pd.read_csv(out/'input/whole_sample.tsv',sep='\t').set_index('Broad')
assert w.to_numpy().sum()==29831
immune=['Myel','T','B','MonoDC','Mac','NKLike','Neut','PlasmaLike']
rr=[]
for sample in w.columns:
 for label in immune:
  rr.append(dict(sample=sample,day=int(sample[1:3]),label=label,cells=w.loc[label,sample],all_QC=w[sample].sum(),immune_QC=w.loc[immune,sample].sum(),pct_all=100*w.loc[label,sample]/w[sample].sum(),pct_immune=100*w.loc[label,sample]/w.loc[immune,sample].sum()))
wlong=pd.DataFrame(rr);wr(wlong,'whole_immune_composition.tsv')
def plotgrid(data,category,value,filename,ylab):
 cats=list(data[category].unique());fig,axes=plt.subplots((len(cats)+3)//4,4,figsize=(12,5.8),squeeze=False)
 for ax,cat in zip(axes.ravel(),cats):
  z=data[data[category]==cat]
  for i,day in enumerate([0,7,28,56]):
   y=z[z.day==day][value].to_numpy();ax.scatter(i+np.linspace(-.08,.08,len(y)),y,s=19,color='#2d718e');ax.plot([i-.18,i+.18],[y.mean()]*2,color='#c55635',lw=2)
  ax.set(title=cat,xticks=range(4),xticklabels=['D0','D7','D28','D56'],ylabel=ylab);ax.spines[['top','right']].set_visible(False)
 for ax in axes.ravel()[len(cats):]:ax.set_visible(False)
 fig.tight_layout()
 for ext in ['png','pdf']:fig.savefig(out/'figures'/f'{filename}.{ext}',dpi=350,bbox_inches='tight')
 plt.close(fig)
for scope in ['TLike','Author']:plotgrid(a[a.scope==scope],'module','score',scope+'_scores','UCell')
plotgrid(wlong,'label','pct_all','immune_all','% all QC')
plotgrid(wlong,'label','pct_immune','immune_only','% immune')
(out/'tables/methods.txt').write_text('GSE274372 supplementation B2. TLike=R08 0/1 (494); Author=888; identity limitations retained. All samples descriptive and exploratory; Min10 sensitivity excludes sample-scope groups with <10 cells. Actual GSE314943 stage8 Ensembl gene members reused, Mre11a removed, UCell maxRank=1500 on full-feature raw counts. Package version recorded; cross-platform absolute scores not compared. Time omnibus and baseline contrasts use independent sample permutations; FDR over all modules/contrasts per scope/rule. Mre correlations use sample ranks and approximate t p values, with categorical-day partial rank sensitivity; small-sample exploratory, no causal inference. Correlation FDR includes both Mre metrics and both adjustments per scope/rule. No contemporaneous later-day controls: time/age/exposure effects cannot be separated. Inflam/Pyro scores describe transcripts, not inflammasome activity or death. Module coverage and component detection must inform interpretation. Whole immune broad labels reused from stage3; denominators all QC cells or selected immune labels, not absolute tissue abundance. Dots independent samples; orange bars means. OxPhos=oxidative phosphorylation; MitoOrg=mitochondrion organization; OxStress=oxidative stress response; Inflam=inflammasome signaling; PyroResp=pyroptotic inflammatory response; PyroExec=execution machinery; InflamRsp=inflammatory response; DDR=DNA damage response. Myel=myeloid; MonoDC=monocyte/DC; Mac=macrophage; NKLike=NK-like; Neut=neutrophil; PlasmaLike=plasma-like.\n')
print('ANALYSIS_PASS',flush=True)
