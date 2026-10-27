import numpy as np
import pandas as pd
from scipy import stats

def canon(g):
 a={'MRE11A':'MRE11','ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','PARK2':'PRKN'}
 s=str(g).upper();return a.get(s,s)

def design(d,model,rank_depth=False):
 cols=[np.ones(len(d)),d.case.to_numpy(float)]
 if model in ['Source','SourceDep','SourceComp']:
  for v in sorted(d.series.unique())[1:]:cols.append((d.series==v).to_numpy(float))
 if model in ['Depth','SourceDep','SourceComp','Comp']:
  x=d.library.to_numpy(float);cols.append(stats.rankdata(x) if rank_depth else np.log1p(x))
 if model in ['SourceComp','Comp']:
  comps=[c for c in d.columns if c.startswith('prop_')]
  for c in comps[:-1]:cols.append(d[c].to_numpy(float))
 return np.column_stack(cols)

def effect(y,d,model):
 y=np.asarray(y,float);assert len(y)==len(d)
 nc=int((d.case==0).sum());nt=int((d.case==1).sum())
 r=dict(n_control=nc,n_case=nt,delta=np.nan,coef=np.nan,p=np.nan,df=np.nan,reason='')
 if min(nc,nt)<2:r['reason']='FEWER_THAN_TWO_PER_GROUP';return r
 r['delta']=float(y[d.case==1].mean()-y[d.case==0].mean())
 X=design(d,model);df=len(y)-np.linalg.matrix_rank(X);r['df']=df
 if np.linalg.matrix_rank(X)<X.shape[1]:r['reason']='RANK_DEFICIENT';return r
 if df<=0:r['reason']='NO_RESIDUAL_DF';return r
 if not np.isfinite(y).all():r['reason']='NONFINITE';return r
 b=np.linalg.lstsq(X,y,rcond=None)[0];r['coef']=b[1]
 s2=float(np.sum((y-X@b)**2)/df);se=np.sqrt(max(0,np.linalg.inv(X.T@X)[1,1]*s2))
 if se>1e-15:r['p']=2*stats.t.sf(abs(b[1]/se),df)
 else:r['reason']='ZERO_RESIDUAL_VARIANCE'
 return r

def association(x,y,d,model):
 x=np.asarray(x,float);y=np.asarray(y,float);assert len(x)==len(y)==len(d)
 r=dict(n=len(x),rho=np.nan,p=np.nan,df=np.nan,reason='')
 if not np.isfinite(x).all() or not np.isfinite(y).all():r['reason']='NONFINITE';return r
 if len(np.unique(x))<2 or len(np.unique(y))<2:r['reason']='CONSTANT';return r
 if model=='Spearman':
  z=stats.spearmanr(x,y);r['rho']=z.statistic
  if len(x)>=10:r['p']=z.pvalue
  r['df']=len(x)-2
 else:
  X=design(d,model,rank_depth=True);df=len(x)-np.linalg.matrix_rank(X)-1;r['df']=df
  if np.linalg.matrix_rank(X)<X.shape[1] or df<=0:r['reason']='UNIDENTIFIABLE';return r
  rx=stats.rankdata(x);ry=stats.rankdata(y);rx=rx-X@np.linalg.lstsq(X,rx,rcond=None)[0];ry=ry-X@np.linalg.lstsq(X,ry,rcond=None)[0]
  if np.linalg.norm(rx)<1e-10 or np.linalg.norm(ry)<1e-10:r['reason']='CONSTANT_RESIDUAL';return r
  rho=float(np.clip(np.corrcoef(rx,ry)[0,1],-1,1));r['rho']=rho
  if len(x)>=10:r['p']=2*stats.t.sf(abs(rho)*np.sqrt(df/max(1-rho*rho,1e-15)),df)
 assert not np.isfinite(r['rho']) or abs(r['rho'])<=1
 return r

def bh(p):
 p=np.asarray(p,float);o=np.argsort(p);q=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];z=np.empty(len(p));z[o]=np.minimum(q,1);return z

def adjust(table,keys,kind):
 table['FDR']=np.nan
 for _,ii in table.groupby(keys,observed=True,dropna=False).groups.items():
  rows=table.loc[ii];ok=rows.p.notna()
  if kind=='effect':ok &= (rows.n_control>=10)&(rows.n_case>=10)
  else:ok &= rows.n>=10
  idx=rows.index[ok]
  if len(idx):table.loc[idx,'FDR']=bh(table.loc[idx,'p'])
 return table

