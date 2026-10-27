from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports');O=P/'breakpoint_TNF_IL18_20261004_v1';R=O/'results';F=P/'front_half_20261004_v1';A=P/'aligned_evidence_20261002_v1'
def rd(p):return pd.read_csv(p,sep='\t',low_memory=False)
def wr(d,n):pd.DataFrame(d).to_csv(R/n,sep='\t',index=False,na_rep='NA')
m=rd(R/'aligned_chain_models.tsv');sign=rd(R/'same_donor_chain_sign_audit.tsv');covers=rd(R/'chain_detection_summary.tsv')
donor=[]
new=rd(R/'new_signature_LOO.tsv.gz')
for _,v in new.iterrows():donor.append(dict(signature=v.signature,omitted_sample=v.omitted_sample,rho=v.rho,n=v.n,df=v.df,status=v.status,origin='new_missing_alignment_LOO'))
extra=rd(R/'historical_missing_LOO_completed.tsv')
for _,v in extra.iterrows():donor.append(dict(signature=v.signature,omitted_sample=v.omitted_sample,rho=v.rho,n=v.n,df=v.df,status=v.status,origin='12_historical_missing_LOO'))
oldf=rd(F/'results/new_link_LOO.tsv.gz')
for _,v in m[(m.origin==F.name)&m.calculation.eq('reused_historical')].drop_duplicates('signature').iterrows():
 for _,r in oldf[oldf.model_id.eq(v.legacy_id)&oldf.omitted_sample.ne('legacy_summary_only')].iterrows():donor.append(dict(signature=v.signature,omitted_sample=r.omitted_sample,rho=r.rho,n=r.n,df=r.df,status=r.status,origin='reused_front_half_LOO'))
oldam=rd(A/'results/aligned_all_models.tsv.gz').drop_duplicates('unique_test_id').set_index('unique_test_id');al=rd(A/'results/aligned_LOO.tsv.gz');al=al[(al.species=='Human')&(al['filter']=='All')&(al.evidence=='SENDER_TF')&(al.tf=='IKZF2')].copy();al['covkey']=al.covariates.fillna('None')
for _,v in m[(m.origin=='aligned_evidence')&m.LOO_detail.ne('new_LOO_only_for_missing_historical_diagnostic')].drop_duplicates('signature').iterrows():
 if v.legacy_id not in oldam.index:continue
 r=oldam.loc[v.legacy_id];cov=r.covariates if isinstance(r.covariates,str) else 'None'
 q=al[(al.scope==r.scope)&(al.seed==r.seed)&(al.covkey==cov)&(al.sample_ids==r.sample_ids)].drop_duplicates('omitted_sample')
 for _,a in q.iterrows():donor.append(dict(signature=v.signature,omitted_sample=a.omitted_sample,rho=a.loo_rho,n=v.n-1,df=a.loo_df,status=a.loo_status,origin='reused_aligned_LOO'))
donor=pd.DataFrame(donor).drop_duplicates(['signature','omitted_sample']);wr(donor,'available_per_donor_influence.tsv.gz');influences=[]
for _,v in m.iterrows():
 q=donor[(donor.signature==v.signature)&donor.rho.notna()].copy();q['delta']=q.rho-v.rho
 influences.append(dict(spec_id=v.spec_id,signature=v.signature,chain=v.chain,edge=v.edge,policy=v.policy,individual_LOO_available=len(q),summary_LOO_valid=v.LOO_valid,largest_change_donor=q.loc[q.delta.abs().idxmax(),'omitted_sample'] if len(q) else '',largest_abs_change=q.delta.abs().max() if len(q) else np.nan,sign_flip_donors=';'.join(q.loc[np.sign(q.rho)!=np.sign(v.rho),'omitted_sample']),note='individual omissions available' if len(q) else 'historical range only; omission names not retained; no invented detail'))
wr(influences,'influence_summary.tsv')
matrix=[]
for _,r in sign[(sign.modality=='AUC17')&sign.policy.eq('aligned_steps')&sign.covariates.eq('G+S+receiver_library+sender_library')].iterrows():
 s,scope,ax=r.chain.split('/');axis=ax.split('-')[0];z=m[(m.chain==r.chain)&(m.modality=='AUC17')&(m.policy=='aligned_steps')&(m.covariates=='G+S+receiver_library+sender_library')];c=covers[(covers.chain==r.chain)&(covers.roster=='AUC17')]
 prox=z[z.edge.eq('ligand_proximal')];classes=['uncertain_directional_mechanism','active_ligand_and_receptor_unmeasured']
 if (c.receptor_codetected_subjects<c.n/2).any():classes.append('sparse_receptor_detection_not_absence')
 base=sign[(sign.chain==r.chain)&(sign.modality=='AUC17')&sign.policy.eq('aligned_steps')&sign.covariates.isna()].iloc[0]
 if abs(base.LT_rho-r.LT_rho)>.25:classes.append('method_covariate_sensitive')
 if r.LT_P>=.05 or r.TFP_P>=.05 or r.LP_P>=.05:classes.append('at_least_one_segment_not_nominally_supported')
 if not r.sign_compatible:classes.append('simple_sign_chain_not_supported_in_this_spec')
 matrix.append(dict(chain=r.chain,modality='AUC17_seed17',adjustment='G+S+receiver_depth+sender_depth',n=r.n,LT_rho=r.LT_rho,LT_P=r.LT_P,LT_FDR=r.LT_FDR,TFP_rho=r.TFP_rho,TFP_P=r.TFP_P,TFP_FDR=r.TFP_FDR,LP_rho=r.LP_rho,LP_P=r.LP_P,LP_FDR=r.LP_FDR,sign_compatible=r.sign_compatible,receptor_detected=';'.join(c.receptor+':'+c.receptor_codetected_subjects.astype(str)+'/17'),proximal=';'.join(prox.y+' rho='+prox.rho.round(4).astype(str)+' P='+prox.P.round(5).astype(str)),evidence_flags=';'.join(classes)))
wr(matrix,'chain_breakpoint_matrix.tsv')
wr(m.groupby(['FDR_family','calculation']).agg(rows=('spec_id','size'),unique_signatures=('signature','nunique'),nominal=('nominal_P_lt05','sum'),FDR05=('FDR_lt05','sum')).reset_index(),'inference_status_summary.tsv')
print('CHAIN_MATRIX');print(pd.DataFrame(matrix).drop(columns=['proximal','evidence_flags']).to_string(index=False))
print('ALL_THREE_NOMINAL');print(sign[sign.all_three_nominal].drop(columns=['sample_ids','interpretation']).to_string(index=False))
print('INFERENCE');print(m.groupby(['FDR_family','policy']).agg(rows=('spec_id','size'),nominal=('nominal_P_lt05','sum'),FDR05=('FDR_lt05','sum')).to_string())
print('IL18_17_EXACT');q=m[(m.chain=='MacAM/CD4-IL7R/IL18-IKZF2')&(m.modality=='AUC17')&m.policy.isin(['aligned_steps','cells20'])];print(q[(q.covariates=='G+S+receiver_library+sender_library')][['spec_id','edge','x','y','n','rho','P','FDR','CI95_low','CI95_high','df','LOO_min','LOO_max','endpoint_VIF','max_covariate_VIF']].to_string(index=False))
print('TNF_CORE_SENS');q=m[(m.chain=='MacAM/Core/TNF-STAT1')&(m.modality=='AUC17')&m.policy.isin(['Treg_composition','TNK_competition','IFN_environment','cells20'])&m.edge.isin(['ligand_TF','TF_phenotype','ligand_phenotype'])];print(q[['policy','edge','y','n','rho','P','FDR','CI95_low','CI95_high','LOO_min','LOO_max']].to_string(index=False))
print('FOCAL_INFLUENCE');print(pd.DataFrame(influences)[lambda d:d.chain.str.startswith('MacAM/')&d.edge.eq('ligand_TF')&d.policy.eq('aligned_steps')].to_string(index=False))
print('SEEDS');print(m[(m.sender=='MacAM')&m.policy.eq('seed_sensitivity')&m.covariates.eq('G+S+receiver_library+sender_library')][['chain','seed','edge','y','rho','P','FDR']].to_string(index=False))
