import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
from pathlib import Path
import json
import pandas as pd,numpy as np
import statslib as S
O=S.O
d=pd.read_csv(O/'tables/sample_scores.tsv',sep='\t');c=pd.read_csv(O/'tables/frozen_candidate_pairs.tsv',sep='\t')
ph=['TAct','TReg','Inhib','Naive','Cycle'];damage=['MitoROS','MitoISR','ApoExec','PyroExec'];aux=['OXPHOS','Inflam','Mitophagy','Glycolysis','Autophagy']
for var in ['Raw','BG1','FixedComp','BG1FixedComp']:
 for cohort in ['Joint','314','274','274_min10']:
  for scope in (['TLike'] if cohort.startswith('274') else ['Core','Naive','ActCD4','Treg']):
   if var in ['FixedComp','BG1FixedComp'] and scope!='Core':continue
   ds=['GSE314943','GSE183682'] if cohort=='Joint' else ['GSE314943'] if cohort=='314' else ['GSE274372']
   z=d[d.dataset.isin(ds)&d.variant.eq(var)&d.scope.eq(scope)].sort_values(['dataset','sample']).reset_index(drop=True)
   if cohort=='274_min10':z=z[z.cells>=10].reset_index(drop=True)
   if len(z)<4:continue
   models=['JointStratum','JointStratumCondition'] if cohort=='Joint' else ['Pooled','Condition'] if cohort=='314' else ['TimePooled','TimeDayAdjusted']
   cc=c if scope in ['TLike','Core'] else c[c.scope==scope]
   cc=cc.sort_values('branch').drop_duplicates('gene')
   for model in models:
    for r in cc.itertuples():
     for ms,kind in [(ph,'phenotype'),(damage,'damage'),(aux,'auxiliary')]:
      for mod in ms:S.assoc(z,'G_'+r.gene,'M_'+mod,cohort,var,'candidate_'+kind+'_'+r.branch,scope,model,extra={'gene':r.gene,'program':mod,'branch':r.branch})
    for a in ph:
     for b in damage:S.assoc(z,'PAIR_'+a+'__'+b+'__'+a,'PAIR_'+a+'__'+b+'__'+b,cohort,var,'phenotype_damage_disjoint',scope,model,extra={'phenotype':a,'program':b})
   for mod in ph+damage+aux:S.group_effect(z,'M_'+mod,cohort,var,scope)
a=S.finalize(S.results,'all_associations.tsv');e=S.finalize(S.effectrows,'program_effects.tsv')
desc=[]
for (scope,var,day),z in d[d.dataset=='GSE183682'].groupby(['scope','variant','day']):
 if len(z)!=2 or set(z.group)!=set(['Silica','Saline']):continue
 for col in [x for x in z if x.startswith(('G_','M_'))]:
  desc.append(dict(scope=scope,variant=var,day=day,variable=col,delta=float(z.loc[z.group=='Silica',col].iloc[0]-z.loc[z.group=='Saline',col].iloc[0]),P=np.nan,status='DESCRIPTIVE_n1+1'))
pd.DataFrame(desc).to_csv(O/'tables/second_cohort_directions.tsv',sep='\t',index=False)
# Exact stratified composition tests, unchanged full cohort.
comp=d[(d.variant=='Raw')&d.scope.isin(['Naive','ActCD4','Treg'])&d.dataset.ne('GSE274372')].copy()
comp['M_proportion']=comp.cells/comp.groupby(['dataset','sample']).cells.transform('sum')
S.effectrows=[]
for sc,z in comp.groupby('scope'):S.group_effect(z.reset_index(drop=True),'M_proportion','Joint','Raw',sc)
S.finalize(S.effectrows,'joint_composition.tsv')
summary={'status':'success','associations':len(a),'effects':len(e),'primary_joint_raw_nominal':int(((a.cohort=='Joint')&(a.variant=='Raw')&(a.model=='JointStratum')&a.kind.str.contains('primary')&(a.P<.05)).sum())}
(O/'logs/associations_complete.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(summary)
