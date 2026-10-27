from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,itertools,datetime
from scipy import stats
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/missing_links_20261003_v1';R=O/'results'
m=pd.read_csv(R/'all_new_models_with_matched.tsv',sep='\t',keep_default_na=False,na_values=['NA',''])
inputs=pd.concat([pd.read_csv(R/'model_sample_inputs.tsv.gz',sep='\t'),pd.read_csv(R/'matched_modality_sample_inputs.tsv.gz',sep='\t')],ignore_index=True)
assert m.model_id.is_unique and not inputs.duplicated(['model_id','sample']).any()
size=inputs.groupby('model_id').size().reindex(m.model_id).fillna(0).to_numpy();assert np.array_equal(size,m.n)
groups={k:v.sort_values('sample') for k,v in inputs.groupby('model_id',sort=False)}
checks=[]
# Independent recomputation for all unadjusted effects and all exact permutation P,
# plus a deterministic sample and the low-df extreme fits among adjusted models.
chosen=m[m.status.eq('FIT')&(m.covariates.eq('None')|m.P_method.str.startswith('exact')|m.rho.abs().gt(.999999))]
chosen=pd.concat([chosen,m[m.status.eq('FIT')].sample(min(250,int(m.status.eq('FIT').sum())),random_state=733)]).drop_duplicates('test_hash')
for row in chosen.itertuples():
 z=groups[row.model_id];n=len(z);cols=[np.ones(n)];covs=[] if row.covariates=='None' else row.covariates.split('+')
 for c in covs:
  if c=='G':cols.append(z.group.isin(['CWP','Silica']).to_numpy(float))
  elif c=='S':cols.extend([z.series.eq(t).to_numpy(float) for t in sorted(z.series.unique())[1:]])
  else:cols.append(stats.rankdata(z[c]))
 D=np.column_stack(cols);rank=np.linalg.matrix_rank(D);assert n-rank-1==row.df
 # Projection with pseudoinverse is independent of runner's least-squares implementation.
 H=np.eye(n)-D@np.linalg.pinv(D);rx=H@stats.rankdata(z.x_value);ry=H@stats.rankdata(z.y_value)
 rho=float(np.dot(rx,ry)/(np.linalg.norm(rx)*np.linalg.norm(ry)));assert abs(rho-row.rho)<2e-10
 dp=np.nan
 if row.P_method.startswith('exact'):
  if row.P_method=='exact_two_sided_rank_permutation':orders=itertools.permutations(range(n))
  else:
   idx=[np.flatnonzero(z.group.eq(v)) for v in sorted(z.group.unique())]
   orders=[]
   for parts in itertools.product(*[list(itertools.permutations(q)) for q in idx]):
    order=np.arange(n)
    for ii,jj in zip(idx,parts):order[ii]=jj
    orders.append(order)
  vals=[np.dot(rx,ry[list(o)])/(np.linalg.norm(rx)*np.linalg.norm(ry)) for o in orders]
  p=np.mean(np.abs(vals)>=abs(rho)-1e-12);dp=abs(p-row.P);assert dp<1e-12
 checks.append(dict(model_id=row.model_id,rho_abs_difference=abs(rho-row.rho),P_abs_difference=dp,method=row.P_method))
pd.DataFrame(checks).to_csv(R/'independent_statistical_checks.tsv',sep='\t',index=False,na_rep='NA')
for family,z in m[m.P.notna()].groupby('FDR_family'):
 u=z.drop_duplicates('test_hash');qq=stats.false_discovery_control(u.P.to_numpy(),method='bh');assert np.allclose(qq,u.FDR,atol=1e-12)
assert (m.P.dropna().between(0,1)).all() and (m.FDR.dropna().between(0,1)).all()
assert m[m.species.eq('Mouse')].n.max()==6
assert not ((m.kind=='sender')&m.axis.isin(['IL15','TNF'])&m.y.str.contains('AUC')).any()
m['near_perfect_low_df']=m.covariates.ne('None')&m.df.le(3)&m.rho.abs().gt(.999999)
m['small_n_adjusted_t_only']=m.n.lt(10)&m.P_method.eq('partial_rank_t_exploratory')
m['interpretation_flag']=np.where(m.near_perfect_low_df,'UNSTABLE_NEAR_PERFECT_LOW_DF_NOT_SUPPORT',np.where(m.small_n_adjusted_t_only,'SMALL_N_APPROXIMATE_P_DESCRIPTIVE_ONLY',np.where(m.status.ne('FIT'),'NOT_ESTIMABLE','EXPLORATORY_ASSOCIATION_NOT_CAUSAL')))
m.to_csv(R/'all_new_models_final.tsv',sep='\t',index=False,na_rep='NA')
# Derivative table uses final sensitivity-family FDR, not an inherited parent-row FDR.
m[m.model_id.str.startswith('MATCH')].to_csv(R/'matched_modality_models.tsv',sep='\t',index=False,na_rep='NA')
m[m.near_perfect_low_df|m.small_n_adjusted_t_only].to_csv(R/'small_sample_method_flags.tsv',sep='\t',index=False,na_rep='NA')
base=m[m.threshold.eq(1)&m.sample_policy.eq('available')&m.covariates.eq('None')&~m.x.str.endswith(('_29','_43'))&~m.y.str.endswith(('_29','_43'))]
base.to_csv(R/'broad_baseline_all_endpoints.tsv',sep='\t',index=False,na_rep='NA')
focused=m[m.species.eq('Human')&m.scope.eq('CD4-IL7R')&m.x.str.startswith('IKZF2')&m.y.str.startswith('OXPHOS')]
focused.to_csv(R/'human_IL7R_IKZF2_OXPHOS_sensitivity.tsv',sep='\t',index=False,na_rep='NA')
trans=[]
for key,z in m[m.sample_policy.eq('available')&m.threshold.eq(1)].groupby(['species','scope','sender','axis','x','y']):
 q=z[z.covariates.eq('None')]
 if q.empty:continue
 b=q.iloc[0]
 for _,v in z.iterrows():
  trans.append(dict(species=key[0],scope=key[1],sender=key[2],axis=key[3],x=key[4],y=key[5],covariates=v.covariates,baseline_n=b.n,n=v.n,baseline_rho=b.rho,rho=v.rho,baseline_P=b.P,P=v.P,FDR=v.FDR,sample_changed=b.sample_ids!=v.sample_ids,sign_changed=bool(np.isfinite(b.rho) and np.isfinite(v.rho) and np.sign(b.rho)!=np.sign(v.rho)),status=v.status,interpretation_flag=v.interpretation_flag))
pd.DataFrame(trans).to_csv(R/'adjustment_transitions.tsv',sep='\t',index=False,na_rep='NA')
manifest=pd.read_csv(O/'provenance/input_manifest.tsv',sep='\t');rehash=[]
for v in manifest.itertuples():
 h=hashlib.sha256()
 with open(v.path,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 match=h.hexdigest()==v.sha256;assert match,v.path;rehash.append(dict(path=v.path,sha256=v.sha256,unchanged=match))
pd.DataFrame(rehash).to_csv(O/'provenance/input_recheck.tsv',sep='\t',index=False)
summary=dict(status='PASS',time=datetime.datetime.now().isoformat(),model_rows=len(m),primary_unique_tests=int(m[m.FDR_family.eq('primary_new')&m.P.notna()].test_hash.nunique()),sensitivity_unique_tests=int(m[m.FDR_family.eq('sensitivity')&m.P.notna()].test_hash.nunique()),status_counts=m.status.value_counts().to_dict(),independently_checked_unique_tests=len(checks),max_rho_check_difference=float(pd.DataFrame(checks).rho_abs_difference.max()),exact_P_checked=int(pd.DataFrame(checks).P_abs_difference.notna().sum()),BH_independent_check=True,sample_membership_all_models_pass=True,input_files_unchanged=len(rehash),near_perfect_low_df_rows=int(m.near_perfect_low_df.sum()),small_n_adjusted_t_rows=int(m.small_n_adjusted_t_only.sum()),old_models_rerun=False,authoritative_table='results/all_new_models_final.tsv')
(O/'provenance/final_verification.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
