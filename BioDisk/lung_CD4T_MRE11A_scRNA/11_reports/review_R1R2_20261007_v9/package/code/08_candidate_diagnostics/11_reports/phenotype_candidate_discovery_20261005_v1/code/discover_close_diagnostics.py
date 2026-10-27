from pathlib import Path
import os,json
os.environ['OPENBLAS_NUM_THREADS']='1'
import pandas as pd,numpy as np
from scipy import stats
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/phenotype_candidate_discovery_20261005_v1';R=O/'results';rd=lambda n:pd.read_csv(R/n,sep='\t',low_memory=False)
d=rd('complete_direct_phenotype_discovery_models.tsv.gz');ex=rd('TXLNA_exact_small_sample_diagnostic.tsv');d['P_exact_corrected']=d.P;d['P_method']='rank_residual_t_approximation'
for v in ex.itertuples():
 ix=(d.ligand==v.ligand)&(d.receiver==v.receiver)&(d.metric==v.metric)&(d.sender==v.sender)&(d.step==v.step);assert ix.sum()==1;d.loc[ix,'P_exact_corrected']=v.exact_P;d.loc[ix,'P_method']='exhaustive_permutation_diagnostic_exact_rho_caution'
d['FDR_exact_corrected']=np.nan
for _,ii in d.groupby(['kind','receiver','metric','denominator','step']).groups.items():
 ix=[i for i in ii if np.isfinite(d.loc[i,'P_exact_corrected'])]
 if ix:d.loc[ix,'FDR_exact_corrected']=stats.false_discovery_control(d.loc[ix,'P_exact_corrected'])
d['FDR_global_exact_corrected']=np.nan;ix=d.index[d.P_exact_corrected.notna()];d.loc[ix,'FDR_global_exact_corrected']=stats.false_discovery_control(d.loc[ix,'P_exact_corrected']);d.to_csv(R/'complete_direct_phenotype_discovery_models.tsv.gz',sep='\t',index=False,na_rep='NA')
q=d[(d.P_method!='rank_residual_t_approximation')|d.FDR_global.lt(.05)];q.to_csv(R/'all_global_hits_exact_diagnostic_comparison.tsv',sep='\t',index=False,na_rep='NA')
t=rd('all_gene_endpoint_final_decisions.tsv.gz');mask=t.ligand.eq('SERPINA1')&t.receiver.eq('CD4-IL7R')&t.metric.eq('score_Mitophagy');t.loc[mask,'shortlist_status']='strong_mixed_myeloid_association_retained;cluster3_nominal;cluster7_source_df_insufficient;MacAM_not_supported;not_promoted_as_pure_macrophage';mask=t.ligand.eq('TXLNA')&t.receiver.eq('Treg')&t.metric.eq('score_OXPHOS');t.loc[mask,'shortlist_status']='exact_permutation_nominal_but_not_full_pool_FDR;source_model_insufficient;no_promotion_from_approximate_extreme_P';t.to_csv(R/'all_gene_endpoint_final_decisions.tsv.gz',sep='\t',index=False,na_rep='NA')
qa=json.loads((O/'provenance/final_audit_QA.json').read_text());qa['global_BH05_after_exact_rank_diagnostic']=int(d.FDR_global_exact_corrected.lt(.05).sum());qa['exact_rank_diagnostic_rows']=len(ex);qa['SERPINA1_identity_audit']='mixed_source_retained;no_pure_MacAM_replication';(O/'provenance/final_audit_QA.json').write_text(json.dumps(qa,indent=2));print(q[['ligand','sender','step','P','P_exact_corrected','FDR_global','FDR_global_exact_corrected']].to_string(index=False))
