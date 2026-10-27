import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
import json
import pandas as pd,numpy as np
import statslib as S
O=S.O
scores=pd.read_csv(O/'tables/sample_scores.tsv',sep='\t');auc=pd.read_csv(O/'tables/TF_sample_AUC.tsv',sep='\t');cand=pd.read_csv(O/'tables/frozen_candidate_pairs.tsv',sep='\t');mem=pd.read_csv(O/'tables/TF_regulon_members.tsv',sep='\t')
focal={'HIF1A','FOXP1','TRPS1','STAT1','STAT3','STAT5A','STAT5B','NFKB1','NFKB2','RELA','RELB','JUN','JUNB','FOS','FOSB','TRP53','NFATC1','NFATC2','FOXO1','FOXO3','ATF4','TFEB'}
mods=['ApoExec','PyroExec','MitoISR','MitoROS','Autophagy','Mitophagy','TAct','TReg','Inhib','Naive','Cycle','OXPHOS','Glycolysis']
routes=[]
for r in cand.itertuples():
 tfs=set(mem.loc[mem.gene.eq(r.gene),'tf'])|{r.gene}
 for tf in sorted(tfs & set(auc.tf)):routes.append(dict(gene=r.gene,scope=r.scope,branch=r.branch,tf=tf,source='self_TF' if tf==r.gene else 'frozen_motif_regulon_membership'))
routes=pd.DataFrame(routes);routes.to_csv(O/'tables/frozen_candidate_TF_routes.tsv',sep='\t',index=False)
for var in ['Raw','BG1','FixedComp','BG1FixedComp']:
 for scope in ['Core','Naive','ActCD4','Treg']:
  if var in ['FixedComp','BG1FixedComp'] and scope!='Core':continue
  z=scores[(scores.variant==var)&(scores.scope==scope)&scores.dataset.ne('GSE274372')].copy()
  a=auc[(auc.variant==var)&(auc.scope==scope)].pivot(index=['dataset','sample'],columns='tf',values='AUC').add_prefix('M_')
  z=z.merge(a,on=['dataset','sample'],how='inner').sort_values(['dataset','sample']).reset_index(drop=True)
  if len(z)!=10:continue
  for tf in sorted(auc.tf.unique()):S.group_effect(z,'M_'+tf,'Joint',var,scope)
  for model in ['JointStratum','JointStratumCondition']:
   for r in routes[routes.scope==scope].itertuples():S.assoc(z,'M_'+r.tf,'G_'+r.gene,'Joint',var,'TF_candidate_'+r.branch,scope,model,extra={'tf':r.tf,'gene':r.gene,'branch':r.branch,'route_source':r.source})
   for tf in sorted(focal & set(auc.tf)):
    for mod in mods:S.assoc(z,'M_'+tf,'M_'+mod,'Joint',var,'focalTF_program',scope,model,extra={'tf':tf,'program':mod})
a=S.finalize(S.results,'TF_associations.tsv');e=S.finalize(S.effectrows,'TF_disease_effects.tsv')
(O/'logs/TF_stats_complete.json').write_text(json.dumps({'status':'success','associations':len(a),'effects':len(e),'candidate_TF_routes':len(routes)},indent=2),encoding='utf-8')
print('TF_STATS',len(a),len(e),len(routes))
