import numpy as np,pandas as pd,itertools
from scipy import stats
PERMS={}
COVCOL={'G':'group','S':'series','R':'receiver_library','D':'sender_library','E':'sender_epi','T':'TregNoTF'}
def bh(p):
 p=np.asarray(p,float);o=np.argsort(p);q=np.empty(len(p));q[o]=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];return np.minimum(q,1)
def design(z,covs):
 cols=[np.ones(len(z))];names=['Intercept']
 for c in covs:
  if c=='G':cols.append(z.group.isin(['CWP','Silica']).to_numpy(float));names.append('Case')
  elif c=='S':
   for v in sorted(z.series.unique())[1:]:cols.append(z.series.eq(v).to_numpy(float));names.append('Source_'+v)
  else:
   q=z[c].to_numpy(float) if c.startswith('prop_') else stats.rankdata(z[COVCOL[c]].to_numpy(float))
   sd=np.std(q);cols.append((q-q.mean())/sd if sd>0 else np.zeros(len(q)));names.append(c)
 return np.column_stack(cols),names
def fit(z,xcol,ycol,covs=(),compute_p=True):
 n=len(z);out=dict(n=n,rho=np.nan,p=np.nan,p_t_approx=np.nan,df=np.nan,status='FIT',p_method='partial_rank_t_exploratory')
 if n<3:out['status']='TOO_FEW_SAMPLES';return out
 if not np.isfinite(z[[xcol,ycol]].to_numpy(float)).all():out['status']='NONFINITE';return out
 X,names=design(z,covs);rank=np.linalg.matrix_rank(X);df=n-rank-1
 out.update(df=df,design_rank=rank,design_columns='|'.join(names),condition_number=float(np.linalg.cond(X)))
 if z[xcol].nunique()<2 or z[ycol].nunique()<2:out['status']='CONSTANT_VARIABLE';return out
 if rank<X.shape[1] or df<=0:out['status']='UNIDENTIFIABLE';return out
 x=stats.rankdata(z[xcol]);y=stats.rankdata(z[ycol]);rx=x-X@np.linalg.lstsq(X,x,rcond=None)[0];ry=y-X@np.linalg.lstsq(X,y,rcond=None)[0]
 if np.linalg.norm(rx)<1e-10 or np.linalg.norm(ry)<1e-10:out['status']='CONSTANT_RESIDUAL';return out
 rho=float(np.clip(np.corrcoef(rx,ry)[0,1],-1,1));out['rho']=rho
 if not compute_p:return out
 p=float(2*stats.t.sf(abs(rho)*np.sqrt(df/max(1-rho*rho,1e-15)),df));out.update(p=p,p_t_approx=p)
 if not covs:
  out['p_method']='rank_correlation_t_exploratory'
  if n<=9:
   if n not in PERMS:PERMS[n]=np.array(list(itertools.permutations(range(n))),dtype=np.int16)
   rx=x-x.mean();ry=y-y.mean();pr=ry[PERMS[n]]@rx/(np.linalg.norm(rx)*np.linalg.norm(ry))
   out['p']=float(np.mean(np.abs(pr)>=abs(rho)-1e-12));out['p_method']='exact_two_sided_rank_permutation'
 elif list(covs)==['G'] and n<=9:
  groups=[np.flatnonzero(z.group.to_numpy()==g) for g in sorted(z.group.unique())]
  permutations=[list(itertools.permutations(ii)) for ii in groups];pr=[]
  for combo in itertools.product(*permutations):
   idx=np.arange(n)
   for ii,jj in zip(groups,combo):idx[ii]=jj
   pr.append(float(np.dot(rx,ry[idx])/(np.linalg.norm(rx)*np.linalg.norm(ry))))
  out['p']=float(np.mean(np.abs(pr)>=abs(rho)-1e-12));out['p_method']='exact_within_group_rank_permutation';out['permutation_count']=len(pr)
 out['near_perfect_low_df']=bool(covs and df<=3 and abs(rho)>.999999)
 out['small_n_or_low_df']=n<10 or df<=3
 return out
def subsets(c):return [v for k in range(len(c)+1) for v in itertools.combinations(c,k)]
