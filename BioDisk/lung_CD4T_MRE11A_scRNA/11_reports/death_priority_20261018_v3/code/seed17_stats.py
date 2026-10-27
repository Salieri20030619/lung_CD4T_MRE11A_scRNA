import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
from pathlib import Path
import sys,json
import pandas as pd,numpy as np
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/death_priority_20261018_v3';V=P/'11_reports/death_priority_20261018_v2'
sys.path.insert(0,str(V/'code'));import statslib as S;S.O=O
a=pd.read_csv(O/'tables/network_seed17_reference_AUC.tsv',sep='\t');s=pd.read_csv(V/'tables/sample_scores.tsv',sep='\t');rts=pd.read_csv(V/'tables/joint_network_candidate_TF_routes.tsv',sep='\t')
for (arm,var,sc),aa in a.groupby(['arm','variant','scope']):
 z=s[(s.variant==var)&(s.scope==sc)&s.dataset.ne('GSE274372')].merge(aa.pivot(index=['dataset','sample'],columns='tf',values='AUC').add_prefix('TF_'),on=['dataset','sample']).sort_values(['dataset','sample']).reset_index(drop=True)
 assert len(z)==10
 for model in ['JointStratum','JointStratumCondition']:
  for r in rts[(rts.arm==arm)&(rts.scope==sc)].itertuples():S.assoc(z,'TF_'+r.tf,'G_'+r.gene,'FullSeed17',var,'candidate_'+r.branch,sc,model,extra={'arm':arm,'gene':r.gene,'tf':r.tf,'branch':r.branch})
ref=S.finalize(S.results,'network_seed17_reference_associations.tsv')
su=pd.read_csv(O/'tables/network_LOO_association_summary.tsv',sep='\t');lo=pd.read_csv(O/'tables/network_LOO_candidate_associations.tsv',sep='\t')
keys=['arm','variant','scope','model','gene','tf','branch']
rr=ref[keys+['status','rho','P','FDR']].rename(columns={c:'full_seed17_'+c for c in ['status','rho','P','FDR']})
out=su.merge(rr,on=keys,how='left',validate='one_to_one')
counts=[]
for r in out.itertuples():
 q=lo[(lo.arm==r.arm)&(lo.scope==r.scope)&(lo.model==r.model)&(lo.gene==r.gene)&(lo.tf==r.tf)&lo.status.eq('OK')]
 counts.append(int((np.sign(q.rho)==np.sign(r.full_seed17_rho)).sum()) if np.isfinite(r.full_seed17_rho) else np.nan)
out['same_direction_to_full_seed17']=counts;out.to_csv(O/'tables/network_LOO_matched_summary.tsv',sep='\t',index=False)
(O/'network_loo/matched_reference_complete.json').write_text(json.dumps({'status':'success','reference_tests':len(ref),'purpose':'separate signature construction from donor omission; original consensus results and LOO statistics unchanged'},indent=2),encoding='utf-8')
print(out[out.gene.eq('IL2RA')&out.tf.eq('NFKB2')].to_string(index=False),flush=True)
