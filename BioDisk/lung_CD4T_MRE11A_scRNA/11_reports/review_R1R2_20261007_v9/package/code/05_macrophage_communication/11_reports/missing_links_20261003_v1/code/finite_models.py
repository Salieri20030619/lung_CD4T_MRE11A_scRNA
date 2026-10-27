"""Run the frozen, finite new hypothesis matrix after coverage review."""
from pathlib import Path
import os, json, itertools, hashlib, fcntl, datetime
import numpy as np
import pandas as pd
from scipy import stats
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/missing_links_20261003_v1')
R=O/'results'
lock=open(O/'analysis.lock','a+');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
assert (O/'provenance/coverage_review.json').exists(),'Coverage review is required before associations'
assert not (O/'provenance/complete.json').exists(),'Refusing completed run overwrite'
rec=pd.read_csv(R/'receiver_paired_modalities.tsv',sep='\t')
pair=pd.read_csv(R/'sender_receiver_paired_inputs.tsv',sep='\t')
sen=pd.read_csv(R/'sender_all_samples.tsv',sep='\t')
comp=sen[sen.sender.eq('TNK_other')][['species','sample','L_IFNG']].rename(columns={'L_IFNG':'TNK_IFNG'})
pair=pair.merge(comp,on=['species','sample'],how='left',validate='many_to_one')
perms={};models=[];members=[];rosters=[];matrix=[]
def bh(p):
 p=np.asarray(p,float);ii=np.argsort(p);q=np.empty(len(p));q[ii]=np.minimum.accumulate((p[ii]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];return np.minimum(q,1)
def design(z,covs):
 cols=[np.ones(len(z))];names=['Intercept']
 for c in covs:
  if c=='G':cols.append(z.group.isin(['CWP','Silica']).to_numpy(float));names.append('Group')
  elif c=='S':
   for val in sorted(z.series.unique())[1:]:cols.append(z.series.eq(val).to_numpy(float));names.append('Source_'+val)
  else:
   v=stats.rankdata(z[c]);sd=np.std(v);cols.append((v-v.mean())/sd if sd>0 else np.zeros(len(v)));names.append(c)
 return np.column_stack(cols),names
def fit(z,xcol,ycol,covs,pvalue=True):
 n=len(z);out=dict(n=n,rho=np.nan,P=np.nan,P_t_approx=np.nan,df=np.nan,rank=np.nan,status='FIT',P_method='not_computed',test_hash='')
 if n<3:out['status']='TOO_FEW_PAIRED_SAMPLES';return out
 D,names=design(z,covs);rank=np.linalg.matrix_rank(D);df=n-rank-1
 out.update(rank=int(rank),df=int(df),design_columns=';'.join(names),condition_number=float(np.linalg.cond(D)))
 rx=stats.rankdata(z[xcol]);ry=stats.rankdata(z[ycol])
 sig=dict(species=z.species.iloc[0],sample=z['sample'].tolist(),x=rx.tolist(),y=ry.tolist(),D=np.round(D,12).tolist())
 out['test_hash']=hashlib.sha256(json.dumps(sig,sort_keys=True).encode()).hexdigest()
 if z[xcol].nunique()<2 or z[ycol].nunique()<2:out['status']='CONSTANT_VARIABLE';return out
 if rank<D.shape[1]:out['status']='RANK_DEFICIENT';return out
 if df<1 or (covs and df<2):out['status']='INSUFFICIENT_RESIDUAL_DF';return out
 rx=rx-D@np.linalg.lstsq(D,rx,rcond=None)[0];ry=ry-D@np.linalg.lstsq(D,ry,rcond=None)[0]
 if min(np.linalg.norm(rx),np.linalg.norm(ry))<1e-10:out['status']='CONSTANT_RESIDUAL';return out
 rho=float(np.clip(np.dot(rx,ry)/(np.linalg.norm(rx)*np.linalg.norm(ry)),-1,1));out['rho']=rho
 if not pvalue:return out
 p=float(2*stats.t.sf(abs(rho)*np.sqrt(df/max(1-rho*rho,1e-15)),df));out.update(P=p,P_t_approx=p,P_method='partial_rank_t_exploratory' if covs else 'rank_t_exploratory')
 if n<=9 and not covs:
  if n not in perms:perms[n]=np.asarray(list(itertools.permutations(range(n))),dtype=np.int16)
  pr=ry[perms[n]]@rx/(np.linalg.norm(rx)*np.linalg.norm(ry));out.update(P=float(np.mean(np.abs(pr)>=abs(rho)-1e-12)),P_method='exact_two_sided_rank_permutation',permutations=len(pr))
 elif n<=9 and covs==['G']:
  groups=[np.flatnonzero(z.group.eq(g)) for g in sorted(z.group.unique())];pr=[]
  for combo in itertools.product(*[list(itertools.permutations(i)) for i in groups]):
   ix=np.arange(n)
   for i,j in zip(groups,combo):ix[i]=j
   pr.append(np.dot(rx,ry[ix])/(np.linalg.norm(rx)*np.linalg.norm(ry)))
  out.update(P=float(np.mean(np.abs(pr)>=abs(rho)-1e-12)),P_method='exact_within_group_rank_permutation',permutations=len(pr))
 out['small_n_or_low_df']=bool(n<10 or df<=3)
 return out
def covcols(covs):return ['group' if c=='G' else 'series' if c=='S' else c for c in covs]
def sequence(z,x,y,kind,axis,source='receiver'):
 sp=z.species.iloc[0];scope=z.scope.iloc[0]
 steps=[[],['G']]
 if sp=='Human':steps += [['G','S'],['G','S','receiver_library']]
 else:steps += [['G','receiver_library']]
 base=['G','S','receiver_library'] if sp=='Human' else ['G']
 if kind=='sender':
  if sp=='Human':steps.append(base+['sender_library'])
  if axis in ['IL27','IFNB1'] and source=='AM':steps.append(base+['TNK_IFNG'])
  if y.startswith('STAT1_'):steps.append(base+['IFN_state_full'])
 if kind=='receiver' and x.startswith('IKZF2_'):
  for c in ['TregNoTF_full','Activation_full','Cycle_full']:
   if not y.startswith(c.split('_')[0]):steps.append(base+[c])
  if scope=='Core':steps.append(base+['Treg_fraction'])
 # Exact step list frozen; no combinatorial covariate search.
 steps=list(dict.fromkeys(tuple(v) for v in steps));steps=[list(v) for v in steps]
 matrix.append(dict(species=sp,scope=scope,kind=kind,axis=axis,sender=source,x=x,y=y,steps=' | '.join('+'.join(v) or 'None' for v in steps)))
 for threshold in [1,20]:
  eligibility=z.receiver_cells.ge(threshold)
  if kind=='sender':eligibility &= z.sender_cells.ge(threshold)
  for row in z.itertuples():
   include=row.receiver_cells>=threshold and (kind!='sender' or row.sender_cells>=threshold)
   rosters.append(dict(species=sp,scope=scope,kind=kind,axis=axis,sender=source,x=x,y=y,threshold=threshold,sample=row.sample,included=include,receiver_cells=row.receiver_cells,sender_cells=getattr(row,'sender_cells',np.nan),reason='' if include else 'prespecified_cell_count_threshold'))
  q=z[eligibility].copy().sort_values('sample')
  allcols=list(dict.fromkeys([x,y]+[c for v in steps for c in covcols(v)]))
  fixed=q.dropna(subset=allcols)
  for policy in ['available','fixed_complete']:
   for covs in steps:
    cols=list(dict.fromkeys([x,y]+covcols(covs)))
    zz=(q if policy=='available' else fixed).dropna(subset=cols).copy()
    out=fit(zz,x,y,covs)
    base_sens=('NoRegUnion' in y or x.endswith(('_29','_43')) or y.endswith(('_29','_43')) or threshold!=1 or policy!='available')
    family='sensitivity' if base_sens else 'primary_new'
    row=dict(model_id='NEW'+str(len(models)+1).zfill(5),species=sp,scope=scope,kind=kind,axis=axis,sender=source,x=x,y=y,covariates='+'.join(covs) or 'None',threshold=threshold,sample_policy=policy,eligible_n=len(q),missing_n=len(q)-len(zz),sample_ids=';'.join(zz['sample']),excluded_ids=';'.join(sorted(set(q['sample'])-set(zz['sample']))),FDR_family=family,**out)
    effects=[]
    if out['status']=='FIT':
     effects=[fit(zz[zz['sample']!=sample],x,y,covs,False)['rho'] for sample in zz['sample']];effects=[v for v in effects if np.isfinite(v)]
    row.update(LOO_valid=len(effects),LOO_rho_min=min(effects) if effects else np.nan,LOO_rho_max=max(effects) if effects else np.nan,LOO_same_direction=bool(effects and all(np.sign(v)==np.sign(out['rho']) for v in effects)))
    models.append(row)
    for _,v in zz.iterrows():
     member=dict(model_id=row['model_id'],species=sp,sample=v['sample'],group=v.group,series=v.series,x_value=v[x],y_value=v[y],receiver_cells=v.receiver_cells)
     for c in covcols(covs):member[c]=v[c]
     members.append(member)

axes={'IL15':('L_IL15_complex',['STAT5_proxy_full','STAT5_proxy_NoRegUnion']),
 'TNF':('L_TNF',['NFkB_proxy_full','NFkB_proxy_NoRegUnion','RELB_components_full','RELB_components_NoRegUnion']),
 'IL18':('L_IL18',['NFkB_proxy_full','NFkB_proxy_NoRegUnion','AP1_components_full','AP1_components_NoRegUnion','IFNG_RNA']),
 'CXCL16':('L_CXCL16',['R_CXCL16']),
 'IL27':('L_IL27_complex',['STAT1_RNA','STAT1_AUC_17','STAT1_AUC_29','STAT1_AUC_43','IFN_state_full','IFN_state_NoRegUnion']),
 'IFNB1':('L_IFNB1',['STAT1_RNA','STAT1_AUC_17','STAT1_AUC_29','STAT1_AUC_43','IFN_state_full','IFN_state_NoRegUnion']),
 'IFNG':('L_IFNG',['STAT1_RNA','STAT1_AUC_17','STAT1_AUC_29','STAT1_AUC_43','IFN_state_full','IFN_state_NoRegUnion'])}
for (sp,scope,source),z in pair.groupby(['species','scope','sender']):
 print('SENDER_MATRIX',sp,scope,source,flush=True)
 for axis,(x,ys) in axes.items():
  for y in ys:sequence(z,x,y,'sender',axis,source)
for (sp,scope),z in rec.groupby(['species','scope']):
 print('STATE_MATRIX',sp,scope,flush=True)
 for tf in ['STAT1','IKZF2']:
  for x in [tf+'_RNA']+[tf+'_AUC_'+str(seed) for seed in [17,29,43]]:
   names=['IFN_state'] if tf=='STAT1' else ['TregNoTF','Activation','Cycle','OXPHOS_nuclear']
   for name in names:
    for variant in ['full','NoRegUnion']:sequence(z,x,name+'_'+variant,'receiver',tf)
m=pd.DataFrame(models);m['FDR']=np.nan
for family,inds in m.groupby('FDR_family').groups.items():
 u=m.loc[inds];u=u[u.P.notna()].drop_duplicates('test_hash')
 if len(u):
  mapping=dict(zip(u.test_hash,bh(u.P)));m.loc[inds,'FDR']=m.loc[inds,'test_hash'].map(mapping)
m['nominal_positive']=m.P.lt(.05);m['FDR_positive']=m.FDR.lt(.05)
m.to_csv(R/'all_new_models.tsv',sep='\t',index=False,na_rep='NA')
pd.DataFrame(members).to_csv(R/'model_sample_inputs.tsv.gz',sep='\t',index=False,na_rep='NA',compression='gzip')
pd.DataFrame(rosters).drop_duplicates().to_csv(R/'cell_filter_rosters.tsv',sep='\t',index=False,na_rep='NA')
pd.DataFrame(matrix).to_csv(R/'finite_test_matrix.tsv',sep='\t',index=False)
summary=[]
for keys,q in m.groupby(['FDR_family','species','kind','axis','status'],dropna=False):
 summary.append(dict(zip(['FDR_family','species','kind','axis','status'],keys),rows=len(q),unique_tests=q.test_hash.replace('',np.nan).nunique(),nominal=int(q.P.lt(.05).sum()),FDR_positive=int(q.FDR.lt(.05).sum())))
pd.DataFrame(summary).to_csv(R/'model_status_summary.tsv',sep='\t',index=False)
assert m.model_id.is_unique and not m[(m.species=='Mouse')&m.n.gt(6)].shape[0]
assert not m[(m.kind=='sender')&m.axis.isin(['IL15','TNF'])&m.y.str.contains('AUC')].shape[0]
assert not m[(m.status=='FIT')&m.covariates.ne('None')&m.df.lt(2)].shape[0]
(O/'provenance/models_complete.json').write_text(json.dumps(dict(time=datetime.datetime.now().isoformat(),rows=len(m),status=m.status.value_counts().to_dict(),unique_primary=m[m.FDR_family.eq('primary_new')&m.P.notna()].test_hash.nunique(),unique_sensitivity=m[m.FDR_family.eq('sensitivity')&m.P.notna()].test_hash.nunique(),mouse_max_n=int(m[m.species.eq('Mouse')].n.max()),old_IL15_TNF_AUC_rerun=False),indent=2))
print((O/'provenance/models_complete.json').read_text(),flush=True)
