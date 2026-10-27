import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
import pandas as pd,numpy as np,json
import statslib as S
O=S.O
scores=pd.read_csv(O/'tables/sample_scores.tsv',sep='\t');auc=pd.read_csv(O/'tables/joint_network_TF_sample_AUC.tsv',sep='\t');cand=pd.read_csv(O/'tables/frozen_candidate_pairs.tsv',sep='\t');edge=pd.read_csv(O/'tables/joint_network_candidate_edges.tsv',sep='\t');cov=pd.read_csv(O/'tables/joint_network_TF_coverage.tsv',sep='\t')
focal={'HIF1A','FOXP1','TRPS1','STAT1','STAT3','STAT5A','STAT5B','NFKB1','NFKB2','RELA','RELB','JUN','JUNB','FOS','FOSB','TRP53','NFATC1','NFATC2','FOXO1','FOXO3','ATF4','TFEB'}
mods=['ApoExec','PyroExec','MitoISR','MitoROS','Autophagy','Mitophagy','TAct','TReg','Inhib','Naive','Cycle','OXPHOS','Glycolysis']
routes=[]
for arm,var in [('All','Raw'),('BG1','BG1')]:
 valid=set(cov.loc[(cov.arm==arm)&(cov.status=='SCORED'),'tf'])
 for r in cand.itertuples():
  es=edge[(edge.arm==arm)&(edge.target==r.gene)&edge.stable_2of3];tfs=set(es.tf)|{r.gene}
  for tf in sorted(tfs&valid):routes.append(dict(variant=var,arm=arm,gene=r.gene,scope=r.scope,branch=r.branch,tf=tf,source='self_TF' if tf==r.gene else 'joint_2of3_motif_edge'))
routes=pd.DataFrame(routes);routes.to_csv(O/'tables/joint_network_candidate_TF_routes.tsv',sep='\t',index=False)
for var in ['Raw','BG1','FixedComp','BG1FixedComp']:
 for scope in ['Core','Naive','ActCD4','Treg']:
  if var in ['FixedComp','BG1FixedComp'] and scope!='Core':continue
  z=scores[(scores.variant==var)&(scores.scope==scope)&scores.dataset.ne('GSE274372')].copy();aa=auc[(auc.variant==var)&(auc.scope==scope)].pivot(index=['dataset','sample'],columns='tf',values='AUC')
  valid=set(aa.columns);aa=aa.add_prefix('TF_');z=z.merge(aa,on=['dataset','sample'],how='inner').sort_values(['dataset','sample']).reset_index(drop=True);assert len(z)==10
  for tf in sorted(valid):
   zz=z.copy();zz['M_'+tf]=zz['TF_'+tf];S.group_effect(zz,'M_'+tf,'Joint',var,scope)
  for model in ['JointStratum','JointStratumCondition']:
   for r in routes[(routes.variant==var)&(routes.scope==scope)].itertuples():S.assoc(z,'TF_'+r.tf,'G_'+r.gene,'Joint',var,'jointnet_TF_candidate_'+r.branch,scope,model,extra={'tf':r.tf,'gene':r.gene,'branch':r.branch})
   for tf in sorted(focal&valid):
    for mod in mods:S.assoc(z,'TF_'+tf,'M_'+mod,'Joint',var,'jointnet_focalTF_program',scope,model,extra={'tf':tf,'program':mod})
for r in S.effectrows:r['family']='joint_network|'+r['family'].replace('programme_disease_effect','TF_disease_effect')
S.finalize(S.results,'joint_network_TF_associations.tsv');S.finalize(S.effectrows,'joint_network_TF_disease_effects.tsv')
# Ligand-to-TF edges reuse the fixed NicheNet candidate routes, never join cross-source/state.
S.results=[];S.testids=set()
links=pd.read_csv(O/'tables/candidate_NicheNet_prior_routes.tsv',sep='\t');lr=links[['source','scope','gene','ligand','branch']].drop_duplicates();sender=pd.read_csv(O/'tables/sender_expression.tsv',sep='\t')
for (source,scope),rr in lr.groupby(['source','scope']):
 ss=sender[sender.source==source].pivot(index=['dataset','sample'],columns='gene',values='CP10k')
 for var in ['Raw','BG1']:
  zz=scores[(scores.scope==scope)&(scores.variant==var)&scores.dataset.ne('GSE274372')].copy();aa=auc[(auc.scope==scope)&(auc.variant==var)].pivot(index=['dataset','sample'],columns='tf',values='AUC').add_prefix('TF_');zz=zz.merge(aa,on=['dataset','sample'])
  ls=np.log1p(ss[sorted(rr.ligand.unique())]).add_prefix('L_');zz=zz.merge(ls,on=['dataset','sample']).sort_values(['dataset','sample']).reset_index(drop=True)
  for r in rr.itertuples():
   for tf in routes.loc[(routes.variant==var)&(routes.scope==scope)&(routes.gene==r.gene),'tf']:
    for model in ['JointStratum','JointStratumCondition']:S.assoc(zz,'L_'+r.ligand,'TF_'+tf,'Joint',var,'jointnet_ligand_TF_'+r.branch,scope,model,source,{'gene':r.gene,'ligand':r.ligand,'tf':tf,'branch':r.branch})
S.finalize(S.results,'joint_network_upstream_TF_associations.tsv')
(O/'joint_network/statistics_complete.json').write_text(json.dumps({'status':'success','candidate_TF_routes':len(routes),'upstream_TF_tests':len(S.results)},indent=2),encoding='utf-8');print('JOINT_NETWORK_STATS_COMPLETE',len(routes),len(S.results))
