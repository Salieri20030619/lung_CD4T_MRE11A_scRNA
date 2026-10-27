from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
from stat_core import fit,bh,subsets,COVCOL
O=Path(__file__).resolve().parents[1];d=pd.read_csv(O/'inputs/aligned_receiver_candidates_sender.tsv',sep='\t')
six=['IL1R2','PENK','TTC39C','FOXP1','TNFRSF9','LAG3'];mods=['OXPHOS','MitoISR','MitoROS','Mitophagy','PyroExec','ApoExec','TAct','TReg','Inhib','Inflam','Naive','Cycle']
rows=[];loo=[];roster=pd.read_csv(O/'results/roster_summary.tsv',sep='\t')
for key,z in d.groupby(['species','filter','scope','tf','seed']):
 sp,arm,scope,tf,seed=key;seed=int(seed)
 rr=roster[(roster.species==sp)&(roster['filter']==arm)&(roster.scope==scope)]
 for ro in rr.itertuples():
  zz=z[z['sample'].isin(ro.sample_ids.split(';'))].copy()
  endpoints=[('TF_CANDIDATE',g,'mean_auc','expr_'+g) for g in six]+[('TF_MODULE',m,'mean_auc','score_'+m) for m in mods]
  if ro.branch=='Paired':endpoints += [('SENDER_TF','IL18','sender_mean_CP10k','mean_auc')]
  for evidence,feature,x,y in endpoints:
   candidates=['G','S','R'] if sp=='Human' else ['G','R']
   if evidence=='SENDER_TF':candidates=['G','S','R','D','E']
   specs=subsets(candidates) if evidence!='SENDER_TF' else list(dict.fromkeys([(),*[(tuple(candidates[:k])) for k in range(1,6)],*[(c,) for c in candidates],*[tuple(a for a in candidates if a!=c) for c in candidates],*[tuple(c for c in candidates if c in candidates[::-1][:k]) for k in range(1,6)]]))
   base=zz.replace([np.inf,-np.inf],np.nan).dropna(subset=[x,y]);complete=base.dropna(subset=[COVCOL[c] for c in candidates])
   for covs in specs:
    av=base.dropna(subset=[COVCOL[c] for c in covs])
    for policy,q in [('available_cases',av),('fixed_complete_cases',complete)]:
     q=q.sort_values('sample');res=fit(q,x,y,covs)
     ident=dict(species=sp,filter=arm,scope=scope,tf=tf,seed=seed,evidence=evidence,feature=feature,covariates='+'.join(covs) or 'None',sample_ids=';'.join(q['sample']))
     uid=hashlib.sha256((json.dumps(ident,sort_keys=True)+q[['sample',x,y]+[COVCOL[c] for c in covs]].to_csv(index=False)).encode()).hexdigest()[:24]
     row=dict(**ident,**res,roster_id=ro.roster_id,branch=ro.branch,threshold=ro.threshold,sample_policy=policy,base_roster_n=len(zz),endpoint_loss=len(zz)-len(base),complete_case_loss=len(av)-len(q),excluded=';'.join(sorted(set(zz['sample'])-set(q['sample']))),n_control=int(q.group.isin(['Control','Saline']).sum()),n_case=int(q.group.isin(['CWP','Silica']).sum()),unique_test_id=uid)
     rows.append(row)
     if policy=='available_cases' and ro.threshold in [1,20] and covs in [(),tuple(candidates),('G','S','R')]:
      for who in q['sample']:
       r=fit(q[q['sample']!=who],x,y,covs,compute_p=False);loo.append(dict(**ident,branch=ro.branch,threshold=ro.threshold,omitted_sample=who,loo_rho=r['rho'],loo_status=r['status'],loo_df=r['df']))
 print('PAIRED_COMPLETED',key,len(rows),flush=True)
r=pd.DataFrame(rows);u=r.drop_duplicates('unique_test_id').copy();u['FDR_family']='Aligned_v1/'+u.species+'/'+u.evidence;u['FDR']=np.nan
for _,idx in u.groupby('FDR_family').groups.items():
 ok=u.loc[idx].p.dropna().index;u.loc[ok,'FDR']=bh(u.loc[ok].p)
r=r.merge(u[['unique_test_id','FDR','FDR_family']],on='unique_test_id',validate='many_to_one');r['nominal_positive']=r.p.lt(.05)
r.to_csv(O/'results/aligned_all_models.tsv.gz',sep='\t',index=False)
pd.DataFrame(loo).to_csv(O/'results/aligned_LOO.tsv.gz',sep='\t',index=False)
focus=r[(r.species=='Human')&(r['filter']=='All')&(r.scope=='CD4-IL7R')&(r.sample_policy=='available_cases')&r.covariates.isin(['None','G','G+S','G+S+R','G+S+R+D','G+S+R+D+E'])&((r.branch=='Paired')|((r.branch=='Receiver')&(r.threshold==20)))]
focus.to_csv(O/'results/first_alignment_focus.tsv',sep='\t',index=False)
methods=dict(rows=len(r),unique_tests=len(u),roster_count=len(roster),fdr='BH per species/evidence over unique data/design tests; old FDR unchanged',covariates=COVCOL,thresholds=[1,5,10,20],scope='Existing four receiver scopes; two TFs; original six genes; 12 frozen modules; fixed human MacAM IL18',p='Exact permutation unadjusted n<=9; otherwise rank t approximation; partial rank t exploratory, low df flagged',constant='NA retained, never zero-filled',primary='Original RNA edgeR retained; correlations do not replace group RNA findings',max_complete_case_loss=int(r.complete_case_loss.max()))
(O/'provenance/aligned_methods.json').write_text(json.dumps(methods,indent=2))
print('ALIGNED_MODELS_PASS',json.dumps(methods),flush=True)
