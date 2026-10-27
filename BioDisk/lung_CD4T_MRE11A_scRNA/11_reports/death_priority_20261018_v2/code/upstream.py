import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
import pandas as pd,numpy as np,json
import statslib as S
O=S.O
links=pd.read_csv(O/'tables/joint_NicheNet_candidate_links.tsv',sep='\t')
cand=pd.read_csv(O/'tables/frozen_candidate_pairs.tsv',sep='\t')
links=links.merge(cand,on=['gene','scope'],how='inner')
eligible=links[links.prior_top250].copy();eligible.to_csv(O/'tables/candidate_NicheNet_prior_routes.tsv',sep='\t',index=False)
routes=eligible[['source','scope','gene','ligand','branch']].drop_duplicates()
scores=pd.read_csv(O/'tables/sample_scores.tsv',sep='\t');send=pd.read_csv(O/'tables/sender_expression.tsv',sep='\t');auc=pd.read_csv(O/'tables/TF_sample_AUC.tsv',sep='\t');tfr=pd.read_csv(O/'tables/frozen_candidate_TF_routes.tsv',sep='\t')
for (source,scope),rr in routes.groupby(['source','scope']):
 ss=send[send.source==source].pivot(index=['dataset','sample'],columns='gene',values='CP10k')
 for var in ['Raw','BG1']:
  d=scores[(scores.scope==scope)&(scores.variant==var)&scores.dataset.ne('GSE274372')].copy()
  au=auc[(auc.scope==scope)&(auc.variant==var)].pivot(index=['dataset','sample'],columns='tf',values='AUC').add_prefix('TF_')
  d=d.merge(au,on=['dataset','sample'],how='left')
  for lig in rr.ligand.unique():d=d.merge(ss[[lig]].rename(columns={lig:'L_'+lig}),on=['dataset','sample'],how='left');d['L_'+lig]=np.log1p(d['L_'+lig])
  d=d.sort_values(['dataset','sample']).reset_index(drop=True);assert len(d)==10
  for r in rr.itertuples():
   for model in ['JointStratum','JointStratumCondition']:
    S.assoc(d,'L_'+r.ligand,'G_'+r.gene,'Joint',var,'ligand_candidate_'+r.branch,scope,model,source,{'gene':r.gene,'ligand':r.ligand,'branch':r.branch})
    for tf in tfr.loc[(tfr.gene==r.gene)&(tfr.scope==scope),'tf']:
     S.assoc(d,'L_'+r.ligand,'TF_'+tf,'Joint',var,'ligand_TF_'+r.branch,scope,model,source,{'gene':r.gene,'ligand':r.ligand,'tf':tf,'branch':r.branch})
if S.results:S.finalize(S.results,'joint_upstream_associations.tsv')
else:pd.DataFrame(columns=['P','FDR','family']).to_csv(O/'tables/joint_upstream_associations.tsv',sep='\t',index=False)
counts=cand.merge(eligible.groupby(['gene','scope']).size().rename('prior_rows'),on=['gene','scope'],how='left').fillna({'prior_rows':0});counts.to_csv(O/'tables/candidate_upstream_coverage.tsv',sep='\t',index=False)
(O/'logs/upstream_associations_complete.json').write_text(json.dumps({'status':'success','prior_routes':len(routes),'associations':len(S.results)},indent=2),encoding='utf-8')
print('UPSTREAM',len(routes),len(S.results))
