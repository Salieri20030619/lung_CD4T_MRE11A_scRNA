from pathlib import Path
import sys,json,itertools
import numpy as np,pandas as pd
from scipy.stats import rankdata,permutation_test,spearmanr
from statsmodels.stats.multitest import multipletests
O=Path(sys.argv[1]);samples=['Saline_1','Saline_2','Saline_3','Silica_1','Silica_2','Silica_3'];rd=lambda f:pd.read_csv(O/f,sep='\t')
act=rd('03_receiver/TF_all_regulon_sample_activity.tsv.gz');phen=rd('01_phenotype/sample_program_scores.tsv');phen=phen[phen['filter']=='All'];rec=rd('03_receiver/receptor_candidate_sample_expression.tsv.gz');sender=rd('02_sender/all_ligand_sample_expression.tsv.gz');tf=rd('03_receiver/IKZF2_STAT1_program_associations.tsv');rna=rd('03_receiver/candidate_RNA_program_associations.tsv');pairs=rd('04_chains/all_sender_receiver_associations.tsv.gz')
matA={(s,k):a.pivot(index='sample',columns='regulon',values='mean_auc').loc[samples] for (s,k),a in act.groupby(['scope','seed'])};matR={s:a.pivot(index='sample',columns='gene',values='CP10k').loc[samples] for s,a in rec.groupby('scope')};matS={s:a.pivot(index='sample',columns='gene',values='CP10k').loc[samples] for s,a in sender.groupby('scope')};matY={s:a.set_index('sample').loc[samples] for s,a in phen.groupby('scope')}
def independent(x,y,adjust):
 x=rankdata(x).astype(float);y=rankdata(y).astype(float)
 if adjust:
  for sl in [slice(0,3),slice(3,6)]:x[sl]-=x[sl].mean();y[sl]-=y[sl].mean()
 else:x-=x.mean();y-=y.mean()
 denom=np.linalg.norm(x)*np.linalg.norm(y)
 if denom<1e-12:return np.nan,np.nan
 rho=float(x@y/denom)
 if adjust:
  def stat(a,b,axis):return np.abs(np.sum(np.concatenate([a,b],axis=axis)*y,axis=axis)/denom)
  p=permutation_test((x[:3],x[3:]),stat,permutation_type='pairings',vectorized=True,n_resamples=np.inf,alternative='greater').pvalue
 else:
  def stat(a,axis):return np.abs(np.sum(a*y,axis=axis)/denom)
  p=permutation_test((x,),stat,permutation_type='pairings',vectorized=True,n_resamples=np.inf,alternative='greater').pvalue
 return rho,float(p)
checks=[]
for r in tf.itertuples():
 rr,p=independent(matA[r.scope,r.seed][r.x].to_numpy(),matY[r.scope][r.y].to_numpy(),r.model=='ConditionAdjusted')
 assert (np.isnan(p) and pd.isna(r.P)) or abs(p-r.P)<1e-10,(r.scope,r.x,r.y,p,r.P)
 if np.isfinite(rr):assert abs(rr-r.rho)<1e-10
 checks.append(dict(kind='TF_program',scope=r.scope,seed=r.seed,x=r.x,y=r.y,model=r.model,max_P_difference=abs(p-r.P) if np.isfinite(p) else np.nan))
for r in rna.itertuples():
 rr,p=independent(matR[r.scope][r.x].to_numpy(),matY[r.scope][r.y].to_numpy(),r.model=='ConditionAdjusted')
 assert (np.isnan(p) and pd.isna(r.P)) or abs(p-r.P)<1e-10
 checks.append(dict(kind='RNA_program',scope=r.scope,seed=0,x=r.x,y=r.y,model=r.model,max_P_difference=abs(p-r.P) if np.isfinite(p) else np.nan))
for r in pairs[pairs.x.isin(['IL15','TNF','IL18'])].itertuples():
 typ,feature=r.y.split('::');y=matA[r.scope,17][feature] if typ=='AUC' else matR[r.scope][feature];rr,p=independent(matS[r.sender][r.x].to_numpy(),y.to_numpy(),r.model=='ConditionAdjusted')
 assert (np.isnan(p) and pd.isna(r.P)) or abs(p-r.P)<1e-10,(r.x,r.y,r.scope,p,r.P)
 checks.append(dict(kind='sender_receiver',scope=r.scope,seed=17,x=r.x,y=r.y,model=r.model,max_P_difference=abs(p-r.P) if np.isfinite(p) else np.nan))
bhchecks=[]
for name in ['03_receiver/all_TF_program_associations.tsv.gz','03_receiver/candidate_RNA_program_associations.tsv','04_chains/all_sender_receiver_associations.tsv.gz','03_receiver/all_TF_activity_group_effects.tsv','06_support/sender_TF_seed_sensitivity.tsv.gz','06_support/Core_fixed_composition_TF_program.tsv','06_support/Core_fixed_composition_RNA_program.tsv','05_human/human_program_group_tests.tsv','05_human/human_candidate_RNA_group_tests.tsv']:
 a=rd(name)
 for fam,b in a.groupby('family'):
  b=b[b.P.notna()]
  if not len(b):continue
  q=multipletests(b.P.to_numpy(),method='fdr_bh')[1];err=float(np.max(abs(q-b.BH.to_numpy())));assert err<1e-10;bhchecks.append(dict(file=name,family=fam,valid_tests=len(b),max_BH_difference=err))
pd.DataFrame(checks).to_csv(O/'provenance/independent_P_validation.tsv',sep='\t',index=False,na_rep='NA');pd.DataFrame(bhchecks).to_csv(O/'provenance/independent_BH_validation.tsv',sep='\t',index=False)
q={'status':'PASS','independent_focal_P_checks':len(checks),'BH_families_checked':len(bhchecks),'max_BH_difference':max(r['max_BH_difference'] for r in bhchecks),'primary_sample_count':6,'conditional_permutations':36,'pooled_permutations':720}
(O/'provenance/independent_validation.json').write_text(json.dumps(q,indent=2),encoding='utf-8');print(json.dumps(q))
b=rd('04_chains/old_ligand_candidate_evidence_segments.tsv.gz');print('CONDITIONAL_DOUBLE');print(b[b.both_edges_condition_nominal][['sender','scope','ligand','receiver_gene','readout','module','pooled_sender_P','pooled_program_P','condition_sender_P','condition_program_P','receiver_DE_P']].to_string(index=False))