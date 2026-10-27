from pathlib import Path
import numpy as np,pandas as pd,json,hashlib,itertools
from scipy import stats

ROOT=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/stepwise_diagnostic_20261002_v1')
I=ROOT/'inputs';O=ROOT/'results'
d=pd.read_csv(I/'receiver_all_positive_cell_samples.tsv',sep='\t');d=d[d.species=='Human']
s=pd.read_csv(I/'sender_IL18_MacAM_all_samples.tsv',sep='\t')
COVORDER=['G','S','R','D','E']
COVCOL={'G':'group','S':'series','R':'receiver_library','D':'sender_library','E':'sender_epi'}
specs={}
def reg(covs,label):
 key=tuple(c for c in COVORDER if c in covs);specs.setdefault(key,[]).append(label)
reg([],'forward_0')
for k in range(1,6):reg(COVORDER[:k],f'forward_{k}')
for c in COVORDER:reg([c],'single_add_'+c);reg([a for a in COVORDER if a!=c],'full_drop_'+c)
for k in range(1,6):reg(COVORDER[::-1][:k],f'reverse_{k}')
PERMS={}
def bh(p):
 p=np.asarray(p,float);o=np.argsort(p);q=np.empty(len(p));q[o]=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];return np.minimum(q,1)
def fit(z,covs):
 n=len(z);cols=[np.ones(n)];names=['Intercept']
 for c in covs:
  if c=='G':cols.append(z.group.eq('CWP').to_numpy(float));names.append('CWP')
  elif c=='S':
   for level in sorted(z.series.unique())[1:]:cols.append(z.series.eq(level).to_numpy(float));names.append('Source_'+level)
  else:
   q=stats.rankdata(z[COVCOL[c]].to_numpy(float));sd=np.std(q)
   cols.append((q-q.mean())/sd if sd>0 else np.zeros(n));names.append(c)
 X=np.column_stack(cols);rank=np.linalg.matrix_rank(X);df=n-rank-1
 r={'n':n,'rho':np.nan,'p':np.nan,'p_t_approx':np.nan,'df':df,'design_rank':rank,'design_columns':'|'.join(names),'status':'FIT','p_method':'partial_rank_t_exploratory','condition_number':float(np.linalg.cond(X))}
 x=stats.rankdata(z.sender_mean_CP10k);y=stats.rankdata(z.mean_auc)
 if rank<X.shape[1] or df<=0:r['status']='UNIDENTIFIABLE';return r
 rx=x-X@np.linalg.lstsq(X,x,rcond=None)[0];ry=y-X@np.linalg.lstsq(X,y,rcond=None)[0]
 if np.linalg.norm(rx)<1e-10 or np.linalg.norm(ry)<1e-10:r['status']='CONSTANT_RESIDUAL';return r
 rho=float(np.clip(np.corrcoef(rx,ry)[0,1],-1,1));approx=float(2*stats.t.sf(abs(rho)*np.sqrt(df/max(1-rho*rho,1e-15)),df))
 r.update(rho=rho,p=approx,p_t_approx=approx,x_residual_fraction=float(np.sum(rx*rx)/np.sum((x-x.mean())**2)),y_residual_fraction=float(np.sum(ry*ry)/np.sum((y-y.mean())**2)))
 full=np.column_stack([X,x]);hat=full@np.linalg.pinv(full);r['max_leverage']=float(np.diag(hat).max())
 if not covs and n<=9:
  if n not in PERMS:PERMS[n]=np.asarray(list(itertools.permutations(range(n))),dtype=np.int16)
  x=x-x.mean();y=y-y.mean();pr=y[PERMS[n]]@x/(np.linalg.norm(x)*np.linalg.norm(y))
  r['p']=float(np.mean(np.abs(pr)>=abs(rho)-1e-12));r['p_method']='exact_two_sided_rank_permutation'
 return r

results=[];rosters=[];paired=[]
def record(z,covs,key,available_n,base_ids,policy,roles):
 z=z.sort_values('sample').copy();fitr=fit(z,covs)
 ident={'scope':key['scope'],'filter':key['filter'],'seed':int(key['seed']),'tf':key['tf'],'covariates':'+'.join(covs) or 'None','sample_ids':';'.join(z['sample'])}
 vals=z[['sample','sender_mean_CP10k','mean_auc','group','series','receiver_library','sender_library','sender_epi']].to_csv(index=False)
 unique=hashlib.sha256((json.dumps(ident,sort_keys=True)+vals).encode()).hexdigest()[:20]
 row=dict(**key,**fitr,covariates=ident['covariates'],roles=';'.join(roles),sample_policy=policy,sample_ids=ident['sample_ids'],excluded_from_threshold_roster=';'.join(sorted(set(base_ids)-set(z['sample']))),available_case_n=available_n,complete_case_loss=available_n-len(z) if policy=='fixed_complete_cases' else 0,n_control=int(z.group.eq('Control').sum()),n_case=int(z.group.eq('CWP').sum()),group_distribution=json.dumps(z.group.value_counts().to_dict(),sort_keys=True),source_distribution=json.dumps(z.series.value_counts().to_dict(),sort_keys=True),receiver_cells_total=int(z.receiver_cells.sum()),sender_cells_total=int(z.sender_cells.sum()),receiver_cells_min=int(z.receiver_cells.min()) if len(z) else np.nan,sender_cells_min=int(z.sender_cells.min()) if len(z) else np.nan,tf_detection_median_pct=float(z.tf_pct.median()),sender_detection_median_pct=float(z.sender_pct.median()),small_sample_or_low_df=len(z)<10 or fitr['df']<=3,unique_test_id=unique)
 row['nominal_positive']=bool(np.isfinite(row['p']) and row['p']<.05);row['model_id']='S'+str(len(results)+1).zfill(5)
 keep=['sample','group','series','receiver_cells','sender_cells','tf_mean_CP10k','tf_pct','sender_mean_CP10k','sender_pct','mean_auc','receiver_library','sender_library','sender_epi']
 rr=z[keep].copy();rr.insert(0,'model_id',row['model_id']);rosters.append(rr);results.append(row)

for (arm,scope,seed,tf),z in d.groupby(['filter','scope','seed','tf'],sort=True):
 for variant,sarm in [('Matched','All'),('NoDom','All'),('NoParent16','NoParent16')]:
  zz=z.merge(s[s.sender_filter==sarm],on='sample',suffixes=('','_sender'),validate='one_to_one')
  assert (zz.group==zz.group_sender).all() and (zz.series==zz.series_sender).all()
  if variant=='NoDom':zz=zz[~zz['sample'].isin(['CTL068','CTL070'])]
  zz=zz.replace([np.inf,-np.inf],np.nan)
  paired.append(zz.assign(variant=variant))
  for threshold in ([1,5,10,20] if variant=='Matched' else [1,20]):
   b=zz[(zz.receiver_cells>=threshold)&(zz.sender_cells>=threshold)].dropna(subset=['sender_mean_CP10k','mean_auc'])
   cc=b.dropna(subset=list(COVCOL.values()))
   for covs,roles in specs.items():
    avail=b.dropna(subset=[COVCOL[c] for c in covs])
    key=dict(species='Human',scope=scope,filter=arm,seed=int(seed),tf=tf,sender='MacAM',ligand='IL18',variant=variant,threshold=threshold,analysis='SENDER_TF_STEPS',omitted_sample='')
    record(cc,covs,key,len(avail),b['sample'],'fixed_complete_cases',roles)
    record(avail,covs,key,len(avail),b['sample'],'available_cases',roles)
  if variant=='Matched' and scope=='CD4-IL7R':
   b=zz[(zz.receiver_cells>=1)&(zz.sender_cells>=1)].dropna(subset=['sender_mean_CP10k','mean_auc']+list(COVCOL.values()))
   for covs in [(),('G','S','R'),('G','S','R','D','E')]:
    for who in b['sample']:
     key=dict(species='Human',scope=scope,filter=arm,seed=int(seed),tf=tf,sender='MacAM',ligand='IL18',variant='Matched',threshold=1,analysis='SENDER_TF_ONE_SAMPLE_REMOVAL',omitted_sample=who)
     record(b[b['sample']!=who],covs,key,len(b)-1,b['sample'],'fixed_complete_cases',['all_subjects_systematic_LOO'])
 print('SENDER_DIAG_DONE',arm,scope,seed,tf,flush=True)
r=pd.DataFrame(results);u=r.drop_duplicates(['analysis','unique_test_id']).copy();u['FDR_exploratory']=np.nan
for _,idx in u.groupby('analysis').groups.items():
 ix=u.loc[idx].p.dropna().index;u.loc[ix,'FDR_exploratory']=bh(u.loc[ix].p)
r=r.merge(u[['analysis','unique_test_id','FDR_exploratory']],on=['analysis','unique_test_id'],validate='many_to_one')
r['FDR_family']='New_stepwise/'+r.analysis
r.to_csv(O/'sender_stepwise_models.tsv',sep='\t',index=False)
pd.concat(rosters,ignore_index=True).to_csv(O/'sender_stepwise_model_samples.tsv.gz',sep='\t',index=False)
pd.concat(paired,ignore_index=True).to_csv(O/'sender_all_paired_samples.tsv.gz',sep='\t',index=False)
focus=r[(r['filter']=='All')&(r.scope=='CD4-IL7R')&(r.variant=='Matched')&(r.sample_policy=='fixed_complete_cases')&(r.analysis=='SENDER_TF_STEPS')]
focus.to_csv(O/'sender_IL7R_All_diagnostic_focus.tsv',sep='\t',index=False)
print('SENDER_STEPS_PASS',len(r),'rows;',len(u),'unique tests')
print('FILTER_ONLY_UNADJUSTED')
print(focus[focus.covariates=='None'][['tf','seed','threshold','n','n_control','n_case','rho','p','df','sample_ids']].to_string(index=False))
print('NESTED_ADJUSTMENT_FIXED_SAMPLE')
print(focus[focus.threshold.isin([1,20])&focus.roles.str.contains('forward')][['tf','seed','threshold','covariates','n','rho','p','df','complete_case_loss']].to_string(index=False))
notes={'covariates':COVCOL,'forward_order':COVORDER,'reverse_order':COVORDER[::-1],'all_single_add_and_full_drop':True,'cell_thresholds':[1,5,10,20],'thresholds_chosen_before_model_results':True,'FDR':'New exploratory BH across all unique planned tests in SENDER_TF_STEPS; identical data/design repetitions counted once. Systematic single-subject omission family separate. Original FDR unchanged.','partial_p':'partial-rank t approximation at positive residual df, small n/low df flagged; no n<10 withholding','baseline_small_n':'unadjusted n<=9 uses exact rank permutation','mean_auc_frozen':True,'max_complete_case_loss':int(r.complete_case_loss.max()),'model_rows':len(r),'unique_tests':len(u),'causal_inference':False}
(O/'sender_stepwise_methods.json').write_text(json.dumps(notes,indent=2))
