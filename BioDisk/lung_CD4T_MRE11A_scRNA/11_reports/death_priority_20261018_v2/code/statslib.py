from pathlib import Path
import pandas as pd,numpy as np,itertools,math
from scipy import stats
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/death_priority_20261018_v2')
pcache={};results=[];effectrows=[];testids=set()
def design(d,model):
 cols=[np.ones(len(d))]
 if model.startswith('Joint'):
  lab=d.dataset.astype(str)+'_'+d.day.astype(str)
  cols.extend(pd.get_dummies(lab,drop_first=True).to_numpy(dtype=float).T)
 if model=='TimeDayAdjusted':cols.extend(pd.get_dummies(d.day.astype(str),drop_first=True).to_numpy(dtype=float).T)
 if 'Condition' in model:cols.append(d.group.eq('Silica').to_numpy(dtype=float))
 return np.column_stack(cols)
def blocks(d,model):
 if model.startswith('Joint'):v=(d.dataset.astype(str)+'_'+d.day.astype(str)).tolist()
 elif model=='TimeDayAdjusted':v=d.day.astype(str).tolist()
 else:v=['all']*len(d)
 if 'Condition' in model:v=[a+'_'+b for a,b in zip(v,d.group.astype(str))]
 return tuple(v)
def permutations(labels):
 if labels in pcache:return pcache[labels]
 n=len(labels);groups=[np.flatnonzero(np.array(labels)==s) for s in sorted(set(labels))];count=math.prod(math.factorial(len(g)) for g in groups)
 if count<=10000:
  out=np.tile(np.arange(n),(count,1));choices=[list(itertools.permutations(g)) for g in groups]
  for i,part in enumerate(itertools.product(*choices)):
   for g,values in zip(groups,part):out[i,g]=values
  exact=True
 else:
  rng=np.random.default_rng(20261017);out=np.tile(np.arange(n),(10000,1))
  for i in range(len(out)):
   for g in groups:out[i,g]=rng.permutation(g)
  exact=False
 pcache[labels]=(out,exact,count);return pcache[labels]
def coefficient(x,y,d,model):
 C=design(d,model);rank=np.linalg.matrix_rank(C)
 if rank!=C.shape[1] or len(d)-rank-1<1:return np.nan,None,None,None,None
 M=np.eye(len(d))-C@np.linalg.pinv(C);rx=M@stats.rankdata(x);ry=M@stats.rankdata(y);den=np.linalg.norm(rx)*np.linalg.norm(ry)
 if den<1e-12:return np.nan,M,rx,ry,rank
 return float(np.clip(np.dot(rx,ry)/den,-1,1)),M,rx,ry,rank

def assoc(d,xcol,ycol,cohort,variant,kind,scope,model,source='none',extra=None):
 key=(cohort,variant,kind,scope,model,source,xcol,ycol)
 if key in testids:return
 testids.add(key);base=dict(cohort=cohort,variant=variant,kind=kind,scope=scope,model=model,sender=source,x=xcol,y=ycol,family='|'.join([cohort,variant,model,kind,source]));base.update(extra or {})
 if xcol not in d or ycol not in d:results.append({**base,'status':'NOT_AVAILABLE','n':0});return
 mask=np.isfinite(d[xcol].to_numpy(float))&np.isfinite(d[ycol].to_numpy(float));z=d.loc[mask].reset_index(drop=True);x=z[xcol].to_numpy(float);y=z[ycol].to_numpy(float);base.update(n=len(z),sample_list=';'.join(z['sample'].astype(str)),dropped_samples=';'.join(d.loc[~mask,'sample'].astype(str)))
 if len(z)<4:results.append({**base,'status':'TOO_FEW_INDEPENDENT_SAMPLES'});return
 r,M,rx,ry,rank=coefficient(x,y,z,model)
 if not np.isfinite(r):results.append({**base,'status':'CONSTANT_OR_NO_RESIDUAL_DF'});return
 ids,exact,total=permutations(blocks(z,model));permy=stats.rankdata(y)[ids]@M;den=np.linalg.norm(permy,axis=1)*np.linalg.norm(rx);valid=den>1e-12;null=(permy[valid]@rx)/den[valid];extreme=int((np.abs(null)>=abs(r)-1e-12).sum());p=extreme/len(null) if exact else (extreme+1)/(len(null)+1);df=len(z)-rank-1;pt=2*stats.t.sf(abs(r)*np.sqrt(df/max(1-r*r,1e-15)),df);loo=[]
 for i in range(len(z)):
  take=np.arange(len(z))!=i;v=coefficient(x[take],y[take],z.loc[take].reset_index(drop=True),model)[0]
  if np.isfinite(v):loo.append(v)
 base.update(status='OK',rho=r,P=p,P_model_t=pt,residual_df=df,permutations_evaluated=len(null),permutations_total=total,permutation_exact_enumeration=exact,loo_valid=len(loo),loo_same_sign=float(np.mean(np.sign(loo)==np.sign(r))) if loo else np.nan,loo_rho_min=min(loo) if loo else np.nan,loo_rho_max=max(loo) if loo else np.nan);results.append(base)

def group_effect(d,col,cohort,variant,scope):
 r=dict(cohort=cohort,variant=variant,scope=scope,program=col[2:],family=cohort+'|'+variant+'|programme_disease_effect');v=d[col].to_numpy(float)
 if not np.isfinite(v).all() or len(np.unique(v))<2:effectrows.append({**r,'status':'CONSTANT_OR_MISSING'});return
 if cohort.startswith('274'):
  C=np.column_stack([np.ones(len(d)),pd.get_dummies(d.day.astype(str),drop_first=True).to_numpy(float)]);b=np.linalg.lstsq(C,v,rcond=None)[0];rss=np.sum((v-C@b)**2);rss0=np.sum((v-v.mean())**2);df1=np.linalg.matrix_rank(C)-1;df2=len(d)-np.linalg.matrix_rank(C);F=max(0,((rss0-rss)/df1)/(rss/df2)) if rss>1e-15 else np.nan;p=stats.f.sf(F,df1,df2) if np.isfinite(F) else np.nan
  effectrows.append({**r,'status':'OK','n':len(d),'P':p,'method':'time_omnibus_F_model','effect':np.nan,'F':F,'n_control':int((d.day==0).sum()),'n_case':int((d.day!=0).sum()),'sample_list':';'.join(d['sample'])});return
 cond=d.group.eq('Silica').to_numpy(float);C=design(d,'JointStratum' if cohort=='Joint' else 'Pooled');M=np.eye(len(d))-C@np.linalg.pinv(C);resid=M@v;xc=M@cond;beta=float(xc@resid/(xc@xc));labels=blocks(d,'JointStratum' if cohort=='Joint' else 'Pooled');groups=[np.flatnonzero(np.array(labels)==g) for g in sorted(set(labels))];choices=[list(itertools.combinations(g,int(cond[g].sum()))) for g in groups];null=[]
 for picked in itertools.product(*choices):
  xx=np.zeros(len(d))
  for ids in picked:xx[list(ids)]=1
  zz=M@xx
  if zz@zz>1e-12:null.append(float(zz@resid/(zz@zz)))
 p=float(np.mean(np.abs(null)>=abs(beta)-1e-12));effectrows.append({**r,'status':'OK','n':len(d),'n_control':int((cond==0).sum()),'n_case':int((cond==1).sum()),'effect':beta,'P':p,'method':'stratified_assignment_permutation','allocations':len(null),'sample_list':';'.join(d['sample'])})

def finalize(rows,name):
 d=pd.DataFrame(rows);d['FDR']=np.nan;d['family_planned_n']=0;d['family_evaluable_n']=0
 for fam,ids in d.groupby('family').groups.items():
  ok=d.loc[ids,'P'].notna();ii=np.array(list(ids))[ok.to_numpy()];d.loc[ids,'family_planned_n']=len(ids);d.loc[ids,'family_evaluable_n']=len(ii)
  if len(ii):d.loc[ii,'FDR']=stats.false_discovery_control(d.loc[ii,'P'].to_numpy())
 d.to_csv(O/'tables'/name,sep='\t',index=False);return d
