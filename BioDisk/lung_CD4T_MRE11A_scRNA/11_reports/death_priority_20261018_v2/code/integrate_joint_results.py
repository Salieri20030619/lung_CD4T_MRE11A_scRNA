from pathlib import Path
import pandas as pd
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/death_priority_20261018_v2')
s=pd.read_csv(O/'tables/candidate_evidence_summary.tsv',sep='\t');j=pd.read_csv(O/'tables/joint_network_candidate_summary.tsv',sep='\t');j=j[j.variant=='Raw'][['gene','scope','tested','nominal','FDR']].rename(columns={'tested':'jointTF_tested','nominal':'jointTF_nominal','FDR':'jointTF_FDR'})
s=s.drop(columns=[c for c in ['jointTF_tested','jointTF_nominal','jointTF_FDR'] if c in s]);s=s.merge(j,on=['gene','scope'],how='left');s.to_csv(O/'tables/candidate_evidence_summary.tsv',sep='\t',index=False)
u=pd.read_csv(O/'tables/joint_upstream_associations.tsv',sep='\t');t=pd.read_csv(O/'tables/joint_network_upstream_TF_associations.tsv',sep='\t');rows=[]
for lig in ['IL18','SIRPB1C','SPP1','IL6','IGF1']:
 a=u[(u.variant=='Raw')&(u.model=='JointStratum')&(u.sender=='AllMac')&(u.scope=='Treg')&(u.gene=='IL2RA')&(u.ligand==lig)&u.kind.eq('ligand_candidate_primary')]
 b=t[(t.variant=='Raw')&(t.model=='JointStratum')&(t.sender=='AllMac')&(t.scope=='Treg')&(t.tf=='NFKB2')&(t.ligand==lig)&t.kind.eq('jointnet_ligand_TF_primary')]
 r={'ligand':lig}
 for prefix,z in [('ligand_IL2RA',a),('ligand_NFKB2',b)]:
  r[prefix+'_status']='evaluated' if len(z) else 'not_in_candidate_prior_route'
  for col in ['rho','P','FDR']:r[prefix+'_'+col]=z[col].iloc[0] if len(z) else float('nan')
 rows.append(r)
pd.DataFrame(rows).to_csv(O/'tables/IL2RA_NFKB2_matched_source_bridge.tsv',sep='\t',index=False)
print(pd.DataFrame(rows).to_string(index=False))
