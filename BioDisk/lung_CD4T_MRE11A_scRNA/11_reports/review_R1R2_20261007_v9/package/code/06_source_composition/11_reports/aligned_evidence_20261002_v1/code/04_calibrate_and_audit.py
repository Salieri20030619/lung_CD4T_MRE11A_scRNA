from pathlib import Path
import pandas as pd,numpy as np,json,shutil
from scipy import stats
from stat_core import fit,bh
O=Path(__file__).resolve().parents[1];p=O/'results/aligned_all_models.tsv.gz'
r=pd.read_csv(p,sep='\t',low_memory=False);r['covariates']=r.covariates.fillna('None');d=pd.read_csv(O/'inputs/aligned_receiver_candidates_sender.tsv',sep='\t')
first=O/'provenance/aligned_before_conditional_exact.tsv.gz'
if not first.exists():shutil.copy2(p,first)
u=r.drop_duplicates('unique_test_id').copy();cal=[]
for ix,row in u[(u.n<=9)&u.covariates.eq('G')&u.status.eq('FIT')].iterrows():
 z=d[(d.species==row.species)&(d['filter']==row['filter'])&(d.scope==row.scope)&(d.tf==row.tf)&(d.seed==row.seed)&d['sample'].isin(row.sample_ids.split(';'))].sort_values('sample')
 x='sender_mean_CP10k' if row.evidence=='SENDER_TF' else 'mean_auc';y='mean_auc' if row.evidence=='SENDER_TF' else ('expr_' if row.evidence=='TF_CANDIDATE' else 'score_')+row.feature
 res=fit(z,x,y,['G']);assert abs(res['rho']-row.rho)<1e-12
 u.loc[ix,'p']=res['p'];u.loc[ix,'p_method']=res['p_method']
 cal.append(dict(unique_test_id=row.unique_test_id,n=row.n,old_partial_t_p=row.p_t_approx,conditional_exact_p=res['p'],permutation_count=res['permutation_count']))
for _,idx in u.groupby('FDR_family').groups.items():
 ok=u.loc[idx].p.dropna().index;u.loc[ok,'FDR']=bh(u.loc[ok].p)
r=r.drop(columns=['p','p_method','FDR']).merge(u[['unique_test_id','p','p_method','FDR']],on='unique_test_id',validate='many_to_one');r['nominal_positive']=r.p.lt(.05)
r.to_csv(p,sep='\t',index=False);pd.DataFrame(cal).to_csv(O/'provenance/conditional_group_exact_calibration.tsv',sep='\t',index=False)
focus=r[(r.species=='Human')&(r['filter']=='All')&(r.scope=='CD4-IL7R')&(r.sample_policy=='available_cases')&r.covariates.isin(['None','G','G+S','G+S+R','G+S+R+D','G+S+R+D+E'])&((r.branch=='Paired')|((r.branch=='Receiver')&(r.threshold==20)))]
focus.to_csv(O/'results/first_alignment_focus.tsv',sep='\t',index=False)
# Independent checks against scipy for all unique unadjusted estimates, and
# against prior stepwise results for identical pre-existing receiver endpoints.
checks=[]
for ix,row in u[u.covariates.eq('None')&u.status.eq('FIT')].iterrows():
 z=d[(d.species==row.species)&(d['filter']==row['filter'])&(d.scope==row.scope)&(d.tf==row.tf)&(d.seed==row.seed)&d['sample'].isin(row.sample_ids.split(';'))]
 x='sender_mean_CP10k' if row.evidence=='SENDER_TF' else 'mean_auc';y='mean_auc' if row.evidence=='SENDER_TF' else ('expr_' if row.evidence=='TF_CANDIDATE' else 'score_')+row.feature
 rho=stats.spearmanr(z[x],z[y]).statistic;assert abs(rho-row.rho)<1e-12
 checks.append(abs(rho-row.rho))
old=pd.read_csv(O/'inputs/prior_receiver_stepwise_models.tsv',sep='\t',low_memory=False);old['covariates']=old.covariates.fillna('None')
cols=['species','filter','scope','tf','seed','feature','evidence','covariates','threshold','sample_policy','sample_ids']
old=old[old.evidence=='TF_MODULE'];old['seed']=old.seed.astype(int)
cmp=r[(r.evidence=='TF_MODULE')&(r.branch=='Receiver')].merge(old,on=cols,suffixes=('_new','_old'),validate='one_to_one')
assert len(cmp)==len(old)==len(r[(r.evidence=='TF_MODULE')&(r.branch=='Receiver')])
assert np.nanmax(np.abs(cmp.rho-cmp.estimate))<1e-10
unchanged=~cmp.p_method_new.eq('exact_within_group_rank_permutation')
assert np.nanmax(np.abs(cmp.loc[unchanged,'p_new']-cmp.loc[unchanged,'p_old']))<1e-10
summary=dict(status='PASS',model_rows=len(r),unique_tests=len(u),unadjusted_scipy_checked=len(checks),max_scipy_rho_difference=max(checks),prior_receiver_rows_compared=len(cmp),max_prior_rho_difference=float(np.nanmax(np.abs(cmp.rho-cmp.estimate))),conditional_exact_unique_tests=len(cal),max_complete_case_loss=int(r.complete_case_loss.max()),max_endpoint_loss=int(r.endpoint_loss.max()),constant_models=int(r.status.eq('CONSTANT_VARIABLE').sum()),unidentifiable_models=int(r.status.eq('UNIDENTIFIABLE').sum()))
(O/'provenance/aligned_validation.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
notes=json.loads((O/'provenance/aligned_methods.json').read_text());notes['p']='Unadjusted n<=9 exact rank permutation; Group-only n<=9 exact within-group rank permutation; all other partial rank P are exploratory t approximation. Original t approximations retained.'
(O/'provenance/aligned_methods.json').write_text(json.dumps(notes,indent=2))
