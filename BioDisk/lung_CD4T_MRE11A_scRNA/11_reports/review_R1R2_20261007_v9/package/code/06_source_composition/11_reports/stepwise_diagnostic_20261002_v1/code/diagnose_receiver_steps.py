from pathlib import Path
import numpy as np,pandas as pd,json,hashlib,itertools
from scipy import stats
ROOT=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/stepwise_diagnostic_20261002_v1')
d=pd.read_csv(ROOT/'inputs/receiver_all_positive_cell_samples.tsv',sep='\t');O=ROOT/'results'
MODULES=['OXPHOS','MitoISR','MitoROS','Mitophagy','PyroExec','ApoExec','TAct','TReg','Inhib','Inflam','Naive','Cycle']
COVCOL={'G':'group','S':'series','R':'receiver_library','T':'TregNoTF'}
PERMS={};rows=[];members=[]
def bh(p):
 p=np.asarray(p,float);o=np.argsort(p);q=np.empty(len(p));q[o]=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];return np.minimum(q,1)
def cov_design(z,covs,ranked):
 cols=[np.ones(len(z))];names=['Intercept']
 for c in covs:
  if c=='G':cols.append(z.group.isin(['CWP','Silica']).to_numpy(float));names.append('Case')
  elif c=='S':
   for v in sorted(z.series.unique())[1:]:cols.append(z.series.eq(v).to_numpy(float));names.append('Source_'+v)
  else:
   q=z[COVCOL[c]].to_numpy(float)
   q=stats.rankdata(q) if ranked else np.log1p(q) if c=='R' else q
   sd=np.std(q);cols.append((q-q.mean())/sd if sd>0 else np.zeros(len(q)));names.append(c)
 return np.column_stack(cols),names
def fit(z,evidence,ycol,covs):
 association=evidence=='TF_MODULE';y=z[ycol].to_numpy(float)
 if evidence=='RNA_GROUP':y=np.log1p(y)
 X,names=cov_design(z,covs,association)
 if not association:X=np.column_stack([X,z.group.isin(['CWP','Silica']).to_numpy(float)]);names+=['Case_effect']
 rank=np.linalg.matrix_rank(X);df=len(z)-rank-(1 if association else 0)
 r=dict(estimate=np.nan,p=np.nan,p_t_approx=np.nan,df=df,design_rank=rank,design_columns='|'.join(names),status='FIT',condition_number=float(np.linalg.cond(X)),p_method='partial_rank_t_exploratory' if association else 'OLS_t_exploratory',effect_type='partial_Spearman_rho' if association else 'case_minus_control_OLS')
 if rank<X.shape[1] or df<=0:r['status']='UNIDENTIFIABLE';return r
 if association:
  x=stats.rankdata(z.mean_auc);y=stats.rankdata(y);rx=x-X@np.linalg.lstsq(X,x,rcond=None)[0];ry=y-X@np.linalg.lstsq(X,y,rcond=None)[0]
  if np.linalg.norm(rx)<1e-10 or np.linalg.norm(ry)<1e-10:r['status']='CONSTANT_RESIDUAL';return r
  rho=float(np.clip(np.corrcoef(rx,ry)[0,1],-1,1));p=float(2*stats.t.sf(abs(rho)*np.sqrt(df/max(1-rho*rho,1e-15)),df))
  r.update(estimate=rho,p=p,p_t_approx=p,x_residual_fraction=float(np.sum(rx*rx)/np.sum((x-x.mean())**2)),y_residual_fraction=float(np.sum(ry*ry)/np.sum((y-y.mean())**2)))
  if not covs and len(z)<=9:
   n=len(z)
   if n not in PERMS:PERMS[n]=np.asarray(list(itertools.permutations(range(n))),dtype=np.int16)
   x-=x.mean();y-=y.mean();pr=y[PERMS[n]]@x/(np.linalg.norm(x)*np.linalg.norm(y))
   r['p']=float(np.mean(np.abs(pr)>=abs(rho)-1e-12));r['p_method']='exact_two_sided_rank_permutation'
 else:
  if np.ptp(y)<1e-14:r['status']='CONSTANT_OUTCOME';return r
  b=np.linalg.lstsq(X,y,rcond=None)[0];res=y-X@b;se=float(np.sqrt(np.sum(res*res)/df*np.linalg.inv(X.T@X)[-1,-1]))
  r.update(estimate=float(b[-1]),se=se)
  if se<=1e-14:r['status']='DEGENERATE_RESIDUAL_VARIANCE';return r
  r.update(p=float(2*stats.t.sf(abs(b[-1]/se),df)),CI_low=float(b[-1]-stats.t.ppf(.975,df)*se),CI_high=float(b[-1]+stats.t.ppf(.975,df)*se))
 return r
def run(z,evidence,feature,tf,seed,ycol):
 sp=z.species.iloc[0];scope=z.scope.iloc[0];arm=z['filter'].iloc[0]
 if evidence=='TF_MODULE':candidates=['G','S','R'] if sp=='Human' else ['G','R']
 elif sp=='Human':candidates=['S','R']
 elif scope=='ActCD4' and evidence in ['RNA_GROUP','AUC_GROUP']:candidates=['R','T']
 else:candidates=['R']
 specs=[q for k in range(len(candidates)+1) for q in itertools.combinations(candidates,k)]
 need=[ycol]+(['mean_auc'] if evidence=='TF_MODULE' else [])
 z=z.replace([np.inf,-np.inf],np.nan)
 for threshold in [1,5,10,20]:
  b=z[z.receiver_cells>=threshold].dropna(subset=need)
  cc=b.dropna(subset=[COVCOL[c] for c in candidates])
  for covs in specs:
   avail=b.dropna(subset=[COVCOL[c] for c in covs])
   for policy,q in [('fixed_complete_cases',cc),('available_cases',avail)]:
    q=q.sort_values('sample').copy();r=fit(q,evidence,ycol,covs)
    ids=';'.join(q['sample']);covname='+'.join(covs) or 'None'
    ident=dict(species=sp,scope=scope,filter=arm,tf=tf,seed=str(seed),evidence=evidence,feature=feature,covariates=covname,sample_ids=ids)
    vals=q[['sample',ycol]+(['mean_auc'] if evidence=='TF_MODULE' else [])+[COVCOL[c] for c in covs]].to_csv(index=False)
    uid=hashlib.sha256((json.dumps(ident,sort_keys=True)+vals).encode()).hexdigest()[:20]
    row=dict(**ident,**r,threshold=threshold,forward_order='+'.join(candidates),all_subsets_for_order_check=True,sample_policy=policy,n=len(q),n_control=int(q.group.isin(['Control','Saline']).sum()),n_case=int(q.group.isin(['CWP','Silica']).sum()),available_case_n=len(avail),complete_case_loss=len(avail)-len(q) if policy=='fixed_complete_cases' else 0,excluded_from_threshold_roster=';'.join(sorted(set(b['sample'])-set(q['sample']))),group_distribution=json.dumps(q.group.value_counts().to_dict(),sort_keys=True),source_distribution=json.dumps(q.series.value_counts().to_dict(),sort_keys=True),receiver_cells_total=int(q.receiver_cells.sum()),receiver_cells_min=int(q.receiver_cells.min()) if len(q) else np.nan,receiver_cells_max=int(q.receiver_cells.max()) if len(q) else np.nan,tf_detection_median_pct=float(q.tf_pct.median()) if tf else np.nan,module_detection_median=float(q['detect_'+feature].median()) if 'detect_'+feature in q else np.nan,small_sample_or_low_df=len(q)<10 or r['df']<=3,unique_test_id=uid)
    row['nominal_positive']=bool(np.isfinite(row['p']) and row['p']<.05);row['model_id']='R'+str(len(rows)+1).zfill(6)
    keep=['sample','group','series','receiver_cells','receiver_library','tf_mean_CP10k','tf_pct','mean_auc']+([ycol] if ycol not in ['mean_auc','tf_mean_CP10k'] else [])+(['TregNoTF'] if 'T' in candidates else [])
    mm=q[keep].copy();mm.insert(0,'model_id',row['model_id']);mm['modeled_feature_value']=q[ycol].to_numpy();members.append(mm);rows.append(row)

for (sp,arm,scope,tf),z in d.groupby(['species','filter','scope','tf'],sort=True):run(z[z.seed==17],'RNA_GROUP',tf,tf,'NA','tf_mean_CP10k')
print('RECEIVER_RNA_STEPS_DONE',len(rows),flush=True)
for (sp,arm,scope,seed,tf),z in d.groupby(['species','filter','scope','seed','tf'],sort=True):
 run(z,'AUC_GROUP',tf,tf,seed,'mean_auc')
 for mod in MODULES:run(z,'TF_MODULE',mod,tf,seed,'score_'+mod)
 print('RECEIVER_TF_STEPS_DONE',sp,arm,scope,seed,tf,len(rows),flush=True)
for (sp,arm,scope),z in d[(d.seed==17)&(d.tf=='IKZF2')].groupby(['species','filter','scope'],sort=True):
 for mod in MODULES:run(z,'MODULE_GROUP',mod,'','NA','score_'+mod)
r=pd.DataFrame(rows);u=r.drop_duplicates('unique_test_id').copy();u['FDR_exploratory']=np.nan
u['FDR_family']='New_stepwise/'+u.species+'/'+u.evidence
for _,idx in u.groupby('FDR_family').groups.items():
 valid=u.loc[idx].p.dropna().index;u.loc[valid,'FDR_exploratory']=bh(u.loc[valid].p)
r=r.merge(u[['unique_test_id','FDR_exploratory','FDR_family']],on='unique_test_id',validate='many_to_one')
r.to_csv(O/'receiver_stepwise_models.tsv',sep='\t',index=False)
pd.concat(members,ignore_index=True).to_csv(O/'receiver_stepwise_model_samples.tsv.gz',sep='\t',index=False)
notes={'status':'PASS','rows':len(r),'unique_tests':len(u),'thresholds':[1,5,10,20],'all_covariate_subsets':'Up to 3 predeclared applicable covariate terms, fully retained to inspect order dependence','Human_group_effects':['Source','receiver log library'],'Human_TF_module':['Group','Source','receiver rank library'],'Mouse_group_TF_ActCD4':['receiver log library','frozen TregNoTF identity score'],'Mouse_group_modules_and_Core_TF':['receiver log library'],'Mouse_TF_module':['Group','receiver rank library'],'sender_depth_and_epi':'Only sender association models, in sender_stepwise_models.tsv; not added to unrelated group or module endpoints','RNA_method':'Exploratory OLS log1p sample meanCP10k, not edgeR replacement','FDR':'New exploratory BH per species/evidence class across unique data/design tests; identical threshold/complete-case repetitions counted once; original FDR unchanged','max_complete_case_loss':int(r.complete_case_loss.max()),'no_models_selected_by_P':True}
(O/'receiver_stepwise_methods.json').write_text(json.dumps(notes,indent=2))
print('RECEIVER_STEPS_PASS',len(r),len(u),'unique',flush=True)
