from pathlib import Path
import os,json,hashlib,fcntl,datetime,itertools
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np,pandas as pd
from scipy import stats
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/front_half_20261004_v1');R=O/'results'
lock=open(O/'analysis.lock','a+');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
assert (O/'provenance/extraction_complete.json').exists() and not (O/'provenance/models_complete.json').exists()
def rd(n):return pd.read_csv(R/n,sep='\t')
def wr(d,n):pd.DataFrame(d).to_csv(R/n,sep='\t',index=False,na_rep='NA')
rec=rd('receiver_audit_inputs.tsv');sen=rd('sender_state_and_ligands.tsv');paired=rd('paired_inputs.tsv');comp=rd('sender_composition.tsv')
panels=['State_'+p for p in ['Treg','Activation','Naive','Inhibitory','Effector']]
models=[];loo=[];members=[];matrix=[];rosters=[]
oldroot=O.parent/'missing_links_20261003_v1/results'
oldmodels=pd.read_csv(oldroot/'all_new_models_final.tsv',sep='\t',low_memory=False)
oldmodels=oldmodels[oldmodels.species.eq('Human')&oldmodels.kind.eq('receiver')]
oldmodels['covariates']=oldmodels.covariates.fillna('None')
oldindex={}
for _,v in oldmodels.iterrows():
 key=(v.scope,v.x,v.y,v.covariates,v.sample_ids)
 oldindex.setdefault(key,v)
score_alias={'State_Treg':'TregNoTF_NoRegUnion','State_Activation':'Activation_NoRegUnion'}
score_audit=rd('duplicate_score_audit.tsv');assert score_audit.identical.all() and score_audit.max_abs_error.eq(0).all()
def needed(covs):return ['group' if c=='G' else 'series' if c=='S' else c for c in covs]
def design(z,covs):
 a=[np.ones(len(z))];names=['Intercept']
 for c in covs:
  if c=='G':a.append(z.group.eq('CWP').to_numpy(float));names.append('Group_CWP')
  elif c=='S':
   for s in sorted(z.series.unique())[1:]:a.append(z.series.eq(s).to_numpy(float));names.append('Source_'+str(s))
  else:
   v=stats.rankdata(z[c]);v=(v-v.mean())/v.std() if v.std()>0 else np.zeros(len(v));a.append(v);names.append(c)
 return np.column_stack(a),names
def fit(z,x,y,covs,pvalue=True):
 n=len(z);out=dict(n=n,rho=np.nan,P=np.nan,status='FIT',df=np.nan,design_rank=np.nan,test_hash='',condition_number=np.nan,design_columns='',P_method='partial_rank_t_exploratory' if covs else 'rank_t_exploratory')
 if n<10:out['status']='FEWER_THAN_10_SUBJECTS';return out
 D,names=design(z,covs);rank=np.linalg.matrix_rank(D);df=n-rank-1;out.update(df=int(df),design_rank=int(rank),design_columns=';'.join(names),condition_number=float(np.linalg.cond(D)))
 rx=stats.rankdata(z[x]);ry=stats.rankdata(z[y]);out['test_hash']=hashlib.sha256(json.dumps([z['sample'].tolist(),rx.tolist(),ry.tolist(),np.round(D,12).tolist()]).encode()).hexdigest()
 if np.ptp(rx)==0 or np.ptp(ry)==0:out['status']='CONSTANT_ENDPOINT';return out
 if rank<D.shape[1]:out['status']='RANK_DEFICIENT';return out
 if df<3:out['status']='RESIDUAL_DF_BELOW_3';return out
 rx=rx-D@np.linalg.lstsq(D,rx,rcond=None)[0];ry=ry-D@np.linalg.lstsq(D,ry,rcond=None)[0]
 if min(np.linalg.norm(rx),np.linalg.norm(ry))<1e-10:out['status']='CONSTANT_RESIDUAL';return out
 rho=float(np.clip(rx@ry/np.linalg.norm(rx)/np.linalg.norm(ry),-1,1));out['rho']=rho
 if pvalue:out['P']=float(2*stats.t.sf(abs(rho)*np.sqrt(df/max(1-rho*rho,1e-15)),df))
 return out
def one(z,x,y,covs,meta,policy,family,eligible):
 cols=list(dict.fromkeys([x,y]+needed(covs)));zz=z.dropna(subset=cols).sort_values('sample').copy()
 key=(meta['scope'],x,score_alias.get(y,''),'+'.join(covs) or 'None',';'.join(zz['sample']))
 prior=oldindex.get(key) if meta['kind']=='TF' else None
 reused=prior is not None and prior.status=='FIT' and len(zz)>=10 and prior.df>=3
 if reused:
  D,names=design(zz,covs);rx=stats.rankdata(zz[x]);ry=stats.rankdata(zz[y])
  signature=hashlib.sha256(json.dumps([zz['sample'].tolist(),rx.tolist(),ry.tolist(),np.round(D,12).tolist()]).encode()).hexdigest()
  out=dict(n=len(zz),rho=prior.rho,P=prior.P,status='FIT',df=prior.df,design_rank=prior['rank'],test_hash=signature,condition_number=prior.condition_number,design_columns=';'.join(names),P_method=prior.P_method)
 else:out=fit(zz,x,y,covs)
 mid='FH'+str(len(models)+1).zfill(5);row=dict(model_id=mid,**meta,x=x,y=y,covariates='+'.join(covs) or 'None',sample_policy=policy,FDR_family=family,eligible_n=len(eligible),sample_ids=';'.join(zz['sample']),excluded_ids=';'.join(sorted(set(eligible['sample'])-set(zz['sample']))),n_case=int(zz.group.eq('CWP').sum()),n_control=int(zz.group.eq('Control').sum()),**out)
 lv=[]
 if out['status']=='FIT' and not reused:
  for s in zz['sample']:
   v=fit(zz[zz['sample'].ne(s)],x,y,covs,False);loo.append(dict(model_id=mid,omitted_sample=s,rho=v['rho'],status=v['status'],n=v['n'],df=v['df']))
   if np.isfinite(v['rho']):lv.append(v['rho'])
 row.update(LOO_valid=len(lv),LOO_min=min(lv) if lv else np.nan,LOO_max=max(lv) if lv else np.nan,LOO_same_direction=bool(lv and all(np.sign(v)==np.sign(out['rho']) for v in lv)),LOO_all_subjects_valid=bool(lv and len(lv)==len(zz)))
 row.update(statistic_origin='reused_missing_links' if reused else 'new_finite_analysis',legacy_model_id=prior.model_id if reused else '',legacy_FDR=prior.FDR if reused else np.nan)
 if reused:
  row.update(LOO_valid=prior.LOO_valid,LOO_min=prior.LOO_rho_min,LOO_max=prior.LOO_rho_max,LOO_same_direction=prior.LOO_same_direction,LOO_all_subjects_valid=prior.LOO_valid==len(zz))
  loo.append(dict(model_id=mid,omitted_sample='legacy_summary_only',rho=np.nan,status='REUSED_RANGE_SEE_MODEL_TABLE',n=np.nan,df=np.nan))
 models.append(row)
 for _,v in zz.iterrows():
  r=dict(model_id=mid,sample=v['sample'],group=v.group,series=v.series,x_value=v[x],y_value=v[y])
  for c in cols:
   if c not in [x,y]:r[c]=v[c]
  members.append(r)
 for _,v in eligible.iterrows():rosters.append(dict(model_id=mid,sample=v['sample'],included=v['sample'] in set(zz['sample']),reason='' if v['sample'] in set(zz['sample']) else policy+'_gate_or_missing'))
def sequence(z,x,y,kind,scope,sender,extra=False):
 z=z.sort_values('sample').copy();assert z['sample'].is_unique
 steps=[[],['G'],['G','S'],['G','S','receiver_library']]
 if kind in ['ligand','sender_state']:steps.append(steps[-1]+['sender_library'])
 if kind=='sender_state':steps.append(steps[-1]+['MacAM_fraction_MacWide'])
 if kind=='TF' and scope=='Core':steps.append(steps[-1]+['Treg_fraction'])
 meta=dict(kind=kind,scope=scope,sender=sender)
 matrix.append(dict(**meta,x=x,y=y,steps=' | '.join('+'.join(c) or 'None' for c in steps)))
 allcols=list(dict.fromkeys([x,y]+sum([needed(c) for c in steps],[])));fixed=z.dropna(subset=allcols)
 family='sensitivity' if sender=='MacWide' else 'primary'
 one(z,x,y,[],meta,'available_baseline',family,z)
 for cov in steps:one(fixed,x,y,cov,meta,'fixed_complete',family,z)
 for policy,take in [('cells20',z.receiver_cells.ge(20)&(z.sender_cells.ge(20) if 'sender_cells' in z else True)),('NoDom',~z['sample'].isin(['CTL068','CTL070']))]:
  one(fixed[fixed['sample'].isin(z.loc[take,'sample'])],x,y,steps[-1],meta,policy,'sensitivity',z)
 if kind=='ligand':one(fixed,x,y,steps[-1]+['TNK_IFNG'],meta,'TNK_competition','sensitivity',z)
 if kind=='TF' and x.endswith('_RNA'):
  matched=fixed[fixed[x.replace('_RNA','_AUC_17')].notna()]
  for cov in steps:one(matched,x,y,cov,meta,'RNA_matched_AUC_roster','sensitivity',z)
for (scope,sender),z in paired.groupby(['scope','sender']):
 for ligand,y in itertools.product(['IL15','TNF','IL18','CXCL16'],panels):sequence(z,'L_'+ligand,y,'ligand',scope,sender)
core=rec[rec.scope.eq('Core')]
cp=core.merge(comp[['sample','MacAM_fraction_candidate','MacWide_fraction_candidate']],on='sample',validate='one_to_one')
for x,y in itertools.product(['MacAM_fraction_candidate','MacWide_fraction_candidate'],panels):sequence(cp,x,y,'composition','Core','composition')
for sender,states in [('MacAM',['AM','C1qMac']),('MonoMac',['SPP1Mac','FOLR2Mac','Mono'])]:
 z=paired[paired.scope.eq('Core')&paired.sender.eq(sender)]
 for x,y in itertools.product(states,['State_Activation','State_Effector','State_Naive']):sequence(z,x,y,'sender_state','Core',sender)
for scope,z in rec.groupby('scope'):
 for tf,mod,y in itertools.product(['IKZF2','STAT1'],['RNA','AUC_17'],panels):sequence(z,tf+'_'+mod,y,'TF',scope,'receiver')
def bh(p):
 p=np.asarray(p,float);ix=np.argsort(p);q=np.empty(len(p));q[ix]=np.minimum.accumulate((p[ix]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];return np.minimum(q,1)
m=pd.DataFrame(models);m['FDR']=np.nan
for fam,idx in m.groupby('FDR_family').groups.items():
 u=m.loc[idx];u=u[u.P.notna()].drop_duplicates('test_hash');q=dict(zip(u.test_hash,bh(u.P)));m.loc[idx,'FDR']=m.loc[idx,'test_hash'].map(q)
m['nominal_P_lt_05']=m.P.lt(.05);m['FDR_lt_05']=m.FDR.lt(.05)
wr(m,'new_link_models.tsv');wr(loo,'new_link_LOO.tsv.gz');wr(members,'new_link_model_inputs.tsv.gz');wr(rosters,'new_link_rosters.tsv.gz');wr(matrix,'finite_test_matrix.tsv')
print('LINK_MODELS',len(m),m.status.value_counts().to_dict(),flush=True)
# Group context OLS, kept distinct from original count model inference.
groupmodels=[];gloo=[];gin=[]
def gfit(z,y,covs):
 out=dict(n=len(z),n_case=int(z.group.eq('CWP').sum()),n_control=int(z.group.eq('Control').sum()),effect=np.nan,P=np.nan,df=np.nan,status='FIT')
 if len(z)<10 or min(out['n_case'],out['n_control'])<2:out['status']='INSUFFICIENT_SUBJECTS_OR_GROUP';return out
 D,names=design(z,['G']+covs);rank=np.linalg.matrix_rank(D);df=len(z)-rank;out.update(df=int(df),design_rank=int(rank),design_columns=';'.join(names))
 if rank<D.shape[1] or df<3:out['status']='UNIDENTIFIABLE';return out
 yy=z[y].to_numpy();b=np.linalg.lstsq(D,yy,rcond=None)[0];v=np.sum((yy-D@b)**2)/df*np.linalg.inv(D.T@D)[1,1];out['effect']=float(b[1])
 if v<=1e-20:out['status']='ZERO_RESIDUAL_VARIANCE';return out
 out['P']=float(2*stats.t.sf(abs(b[1]/np.sqrt(v)),df));return out
def groupseq(z,y,kind,scope,depth=None):
 covs=[[],['S']]+([['S',depth]] if depth else []);z=z.dropna(subset=[y,'group','series']+([depth] if depth else [])).sort_values('sample')
 for cov in covs:
  o=gfit(z,y,cov);mid='GC'+str(len(groupmodels)+1).zfill(4);vv=[]
  for ss in z['sample']:
   v=gfit(z[z['sample'].ne(ss)],y,cov);gloo.append(dict(model_id=mid,omitted_sample=ss,**v))
   if np.isfinite(v['effect']):vv.append(v['effect'])
  groupmodels.append(dict(model_id=mid,kind=kind,scope=scope,y=y,covariates='Group'+('+'+'+'.join(cov) if cov else ''),sample_ids=';'.join(z['sample']),effect_unit='CWP minus Control in original fraction/gene-z score/meanCP10k units',LOO_valid=len(vv),LOO_min=min(vv) if vv else np.nan,LOO_max=max(vv) if vv else np.nan,LOO_same_direction=bool(vv and all(np.sign(v)==np.sign(o['effect']) for v in vv)),**o))
  for _,v in z.iterrows():gin.append(dict(model_id=mid,sample=v['sample'],group=v.group,series=v.series,y_value=v[y],depth_value=v[depth] if depth else np.nan))
for scope,z in rec.groupby('scope'):
 for y in panels:groupseq(z,y,'receiver_state',scope,'receiver_library')
for y in [c for c in core if c.startswith('Fraction_')]:groupseq(core,y,'receiver_composition','Core')
for y in ['MacAM_fraction_candidate','MacWide_fraction_candidate','MacAM_fraction_MacWide']:groupseq(comp,y,'sender_composition','candidate_myeloid')
for scope,z in sen.groupby('sender'):
 for y in ['AM','C1qMac','SPP1Mac','FOLR2Mac','Mono','L_IL15','L_TNF','L_IL18','L_CXCL16']:groupseq(z,y,'sender_state_or_ligand',scope,'sender_library')
gm=pd.DataFrame(groupmodels);gm['FDR']=np.nan;idx=gm.P.notna();gm.loc[idx,'FDR']=bh(gm.loc[idx,'P']);gm['nominal_P_lt_05']=gm.P.lt(.05);gm['FDR_lt_05']=gm.FDR.lt(.05)
wr(gm,'group_context_models.tsv');wr(gloo,'group_context_LOO.tsv');wr(gin,'group_context_inputs.tsv')
summary=m.groupby(['kind','FDR_family','status']).agg(rows=('model_id','size'),unique_tests=('test_hash','nunique'),nominal=('nominal_P_lt_05','sum'),FDR05=('FDR_lt_05','sum')).reset_index();wr(summary,'model_status_summary.tsv')
assert m.model_id.is_unique and m.n.max()<=19
assert not m[(m.status=='FIT')&((m.n<10)|(m.df<3))].shape[0]
assert len(m)==1675 and len(matrix)==185
(O/'provenance/models_complete.json').write_text(json.dumps(dict(time=datetime.datetime.now().isoformat(),new_link_specs=len(m),new_hypotheses=len(matrix),group_context_specs=len(gm),status=m.status.value_counts().to_dict(),reused_statistic_rows=int(m.statistic_origin.eq('reused_missing_links').sum()),no_old_models_rerun=True),indent=2))
print('MODELS_COMPLETE',flush=True)
