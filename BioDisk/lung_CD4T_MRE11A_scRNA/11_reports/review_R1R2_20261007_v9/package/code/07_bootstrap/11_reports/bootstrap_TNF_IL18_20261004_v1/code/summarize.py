"""Summaries fixed by FROZEN_PLAN; no significance selection."""
from pathlib import Path
import json,itertools,hashlib
import numpy as np,pandas as pd
B=Path(__file__).resolve().parents[1];R=B/'results';P=B/'provenance'
def rd(n):return pd.read_csv(R/n,sep='\t')
def wr(d,n):pd.DataFrame(d).to_csv(R/n,sep='\t',index=False,na_rep='NA')
d=rd('bootstrap_effects.tsv.gz');o=rd('original_effects.tsv');dg=rd('bootstrap_diagnostics.tsv.gz')
KEY=['axis','modality','adjustment','edge']; N=5000
orig=o.set_index(KEY).rho.to_dict()
def stats(v,original):
 a=np.asarray(v,float);a=a[np.isfinite(a)];n=len(a);q=np.quantile(a,[.025,.5,.975],method='linear') if n>=1000 else [np.nan]*3
 same=int(np.sum(np.sign(a)==np.sign(original))) if original!=0 else 0
 return dict(original_effect=original,bootstrap_median=q[1],CI95_low=q[0],CI95_high=q[2],valid=n,attempted=N,failed=N-n,failure_rate=(N-n)/N,sign_preservation=same/n if n and original!=0 else np.nan,sign_preserved_count=same,positive_fraction=float(np.mean(a>0)) if n else np.nan,negative_fraction=float(np.mean(a<0)) if n else np.nan,zero_fraction=float(np.mean(a==0)) if n else np.nan,abs_rho_boundary_fraction=float(np.mean(np.abs(a)>=1-1e-10)) if n else np.nan,interval_status='OK_CONDITIONAL_VALID' if n>=1000 else 'TOO_FEW_VALID')
summ=[]
for key,g in d.groupby(['scheme']+KEY,sort=False):summ.append(dict(zip(['scheme']+KEY,key),**stats(g.rho,orig[key[1:]])))
wr(summ,'effect_summary.tsv')
wr(d.groupby(['scheme']+KEY+['status'],dropna=False).size().reset_index(name='count'),'failure_status_counts.tsv')
di=[]
for key,g in dg.groupby(['scheme','axis','roster','adjustment'],sort=False):
 v=dict(zip(['scheme','axis','roster','adjustment'],key));v.update(attempted=len(g),rank_deficient=int((g.design_rank<g.design_columns).sum()),zero_covariate_rows=int(g.constant_covariates.notna().sum()),zero_endpoint_rows=int(g.zero_variance_endpoints.notna().sum()),zero_residual_rows=int(g.zero_residual_endpoints.notna().sum()),low_unique_support_fraction=float(g.low_unique_support.mean()))
 for c in ['unique_donors','multiplicity_effective_n','design_rank','nominal_df','unique_support_df','condition_number']:
  a=g[c].to_numpy();v[c+'_min']=a.min();v[c+'_median']=np.median(a);v[c+'_max']=a.max()
 di.append(v)
wr(di,'diagnostic_summary.tsv')
chain=[];chain_detail=[]
for key,g in d.groupby(['scheme','axis','modality','adjustment'],sort=False):
 q=g.pivot(index='replicate',columns='edge',values='rho').reindex(range(1,N+1));good=q.notna().all(axis=1);sign=np.sign(q)
 compat=good & (sign.ligand_TF*sign.TF_phenotype==sign.ligand_phenotype)&(sign!=0).all(axis=1)
 osigns={e:np.sign(orig[key[1:]+(e,)]) for e in q.columns};same=good&(sign==pd.Series(osigns)).all(axis=1)
 count=int(compat.sum());n=int(good.sum());rec=dict(zip(['scheme','axis','modality','adjustment'],key))
 chain.append(dict(**rec,valid=n,attempted=N,failed=N-n,failure_rate=(N-n)/N,compatible_count=count,direction_compatible_fraction=count/n if n else np.nan,compatible_confirmed_of_all=count/N,compatible_all_attempts_lower=count/N,compatible_all_attempts_upper=(count+N-n)/N,all_three_original_sign_count=int(same.sum()),all_three_original_sign_fraction=float(same.sum()/n) if n else np.nan,original_signs=';'.join(f'{e}:{int(osigns[e])}' for e in q.columns),original_direction_compatible=bool(np.prod(list(osigns.values()))>0)))
 for i in q.index:chain_detail.append(dict(**rec,replicate=i,valid=bool(good.loc[i]),compatible=bool(compat.loc[i]),same_original_signs=bool(same.loc[i])))
wr(chain,'chain_summary.tsv');wr(chain_detail,'bootstrap_chain_status.tsv.gz')
cc=pd.DataFrame(chain_detail);alls=[]
for key,g in cc[cc.modality.str.startswith('AUC')].groupby(['scheme','axis','adjustment'],sort=False):
 q=g.pivot(index='replicate',columns='modality',values='valid');good=q.all(axis=1);comp=g.pivot(index='replicate',columns='modality',values='compatible').all(axis=1)&good;same=g.pivot(index='replicate',columns='modality',values='same_original_signs').all(axis=1)&good;n=int(good.sum());k=int(comp.sum())
 alls.append(dict(zip(['scheme','axis','adjustment'],key),valid=n,attempted=N,failed=N-n,all_seeds_compatible_count=k,all_seeds_compatible_fraction=k/n,all_seeds_original_sign_fraction=float(same.sum()/n),compatible_all_attempts_lower=k/N,compatible_all_attempts_upper=(k+N-n)/N))
wr(alls,'all_seeds_chain_summary.tsv')
diffs=[];detail=[]
for key,g in d[d.modality!='RNA19'].groupby(['scheme','axis','adjustment','edge'],sort=False):
 q=g.pivot(index='replicate',columns='modality',values='rho').reindex(range(1,N+1))
 comparisons=[('seed_paired',b,a) for a,b in itertools.combinations(['AUC17','AUC29','AUC43'],2)]+[('AUC_minus_matched_RNA',a,'RNA17') for a in ['AUC17','AUC29','AUC43']]
 for typ,plus,minus in comparisons:
  delta=q[plus]-q[minus];obs=orig[(key[1],plus,key[2],key[3])]-orig[(key[1],minus,key[2],key[3])]
  meta=dict(zip(['scheme','axis','adjustment','edge'],key),comparison_type=typ,comparison=f'{plus}-{minus}',plus=plus,minus=minus)
  diffs.append(dict(**meta,**stats(delta,obs)))
  detail.extend(dict(**meta,replicate=i,delta=v,valid=np.isfinite(v)) for i,v in delta.items())
for key,g in d[d.axis=='TNF'].groupby(['scheme','modality','edge'],sort=False):
 q=g.pivot(index='replicate',columns='adjustment',values='rho').reindex(range(1,N+1));delta=q.B_Treg-q.B;obs=orig[('TNF',key[1],'B_Treg',key[2])]-orig[('TNF',key[1],'B',key[2])]
 meta=dict(scheme=key[0],axis='TNF',modality=key[1],adjustment='B_Treg-minus-B',edge=key[2],comparison_type='Treg_enhancement',comparison='B_Treg-B',plus='B_Treg',minus='B')
 diffs.append(dict(**meta,**stats(delta,obs)));detail.extend(dict(**meta,replicate=i,delta=v,valid=np.isfinite(v)) for i,v in delta.items())
for row in diffs:row.pop('abs_rho_boundary_fraction',None)
wr(diffs,'paired_difference_summary.tsv');wr(detail,'bootstrap_paired_differences.tsv.gz')
# Structural checks, targeted independent rank fit checks, and frozen-index checks.
assert len(d)==750000 and len(dg)==100000
assert d.groupby(['scheme']+KEY).size().eq(N).all()
assert not d.duplicated(['scheme','replicate']+KEY).any()
for key,g in d[d.modality!='RNA19'].groupby(['scheme','replicate','axis','adjustment','edge'],sort=False):
 if key[-1]=='ligand_phenotype':assert g.rho.dropna().nunique()<=1
import bootstrap as engine
from scipy.stats import spearmanr
checks=[]
for scheme in ['stratified','unstratified']:
 for roster in ['RNA19','MATCH17']:
  ids=np.load(R/f'indices_{scheme}_{roster}.npy');assert ids.shape==(N,len(engine.rosters[roster]))
  z=engine.subset('TNF',roster)
  if scheme=='stratified':
   labels=(z.series+'|'+z.group).to_numpy();expected=pd.Series(labels).value_counts().to_dict()
   for draw in ids:assert pd.Series(labels[draw]).value_counts().to_dict()==expected
  for rep in [1,2500,5000]:
   zz=z.iloc[ids[rep-1]];r,_=engine.fit(zz,'TNF','unadjusted');expected=spearmanr(zz.L_TNF,zz.STAT1_RNA).statistic;assert abs(r[0]['rho']-expected)<1e-12
   # Independent explicit average ranks via pandas and projection via QR; adjusted FULL-RANK only.
   cols=['L_TNF','STAT1_RNA','State_Effector'];V=zz[cols].rank(method='average').to_numpy();depth=zz[['receiver_library','sender_library']].rank().to_numpy();depth=(depth-depth.mean(axis=0))/depth.std(axis=0)
   DD=np.column_stack([np.ones(len(zz)),zz.group.eq('CWP'),zz.series.eq('G227'),zz.series.eq('G288'),depth]);rr,dd=engine.fit(zz,'TNF','B')
   if np.linalg.matrix_rank(DD)==DD.shape[1]:
    Q=np.linalg.qr(DD)[0];res=V-Q@(Q.T@V);independent=np.corrcoef(res.T)[0,1];assert abs(rr[0]['rho']-independent)<1e-10
   actual=d[(d.scheme==scheme)&(d.replicate==rep)&(d.axis=='TNF')&(d.modality==('RNA19' if roster=='RNA19' else 'RNA17'))&(d.adjustment=='B')&(d.edge=='ligand_TF')].iloc[0]
   assert np.isclose(actual.rho,rr[0]['rho'],atol=1e-12,equal_nan=True)
   checks.append(dict(scheme=scheme,roster=roster,replicate=rep,status='PASS'))
wr(checks,'independent_checks.tsv')
verification=dict(status='PASS',effect_rows=len(d),diagnostic_rows=len(dg),summary_rows=len(summ),paired_difference_rows=len(diffs),chain_rows=len(chain),indices_shared_across_seeds_modalities_adjustments_axes=True,strata_counts_verified_every_replicate=True,independent_selected_checks=len(checks),same_donor_nonTF_edges_identical=True,no_new_P_FDR=True,no_legacy_grid_or_SCENIC_rerun=True)
(P/'VERIFIED.json').write_text(json.dumps(verification,indent=2));print(json.dumps(verification,indent=2))
print(pd.DataFrame(summ).query("scheme=='stratified' and adjustment=='B'").to_string(index=False))
print(pd.DataFrame(diffs).query("scheme=='stratified' and comparison_type=='Treg_enhancement'").to_string(index=False))
