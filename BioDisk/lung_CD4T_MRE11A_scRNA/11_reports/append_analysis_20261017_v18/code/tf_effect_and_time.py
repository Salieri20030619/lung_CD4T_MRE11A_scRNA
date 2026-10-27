from pathlib import Path
import pandas as pd,numpy as np,itertools,json
from scipy import stats
O=Path(r'C:\Users\31937\lung_append_analysis_20261017_v18');a=pd.read_csv(O/'tables/clean_TF_sample_AUC.tsv',sep='\t');plan=json.loads((O/'protocol_frozen.json').read_text());pairs=set()
for g,sc,l,tfs in plan['routes']:
 if tfs not in ['none','decoy']:
  for tf in tfs.split(';'):pairs.add((sc,tf))
pairs|={('Core','FOXP1'),('Treg','FOXP1'),('Treg','TRP53'),('Naive','NFATC1'),('Naive','NFATC2')};rows=[]
for variant in ['Raw','BG1','FixedComp','BG1FixedComp']:
 for scope,tf in sorted(pairs):
  if variant in ['FixedComp','BG1FixedComp'] and scope!='Core':continue
  z=a[(a.variant==variant)&(a.scope==scope)&(a.tf==tf)].sort_values('sample');r=dict(variant=variant,scope=scope,tf=tf,family=variant+'|focal_TF_activity_disease',n=len(z))
  if len(z)!=6:rows.append({**r,'status':'NOT_EVALUABLE'});continue
  x=z.AUC.to_numpy();case=z.group.eq('Silica').to_numpy();delta=x[case].mean()-x[~case].mean();null=[]
  for comb in itertools.combinations(range(6),3):
   mask=np.zeros(6,bool);mask[list(comb)]=True;null.append(x[mask].mean()-x[~mask].mean())
  rows.append({**r,'status':'OK','effect':delta,'P':float(np.mean(np.abs(null)>=abs(delta)-1e-12)),'assignments':20})
d=pd.DataFrame(rows);d['FDR']=np.nan
for fam,z in d.groupby('family'):
 ii=z.index[z.P.notna()];d.loc[ii,'FDR']=stats.false_discovery_control(d.loc[ii,'P'].values)
d.to_csv(O/'tables/focal_TF_disease_effects.tsv',sep='\t',index=False)
r=pd.read_csv(O/'tables/all_associations.tsv',sep='\t');print('TIME FOCUS')
print(r[r.cohort.isin(['274','274_min10'])&r.variant.eq('Raw')&r.kind.eq('candidate_damage')&r.gene.isin(['FOXP1','TTC39C','MDFIC','DDIT4','IL1R2','SKIL'])][['cohort','model','gene','program','rho','P','FDR']].to_string(index=False))
print('TF EFFECTS',d[(d.variant=='Raw')][['scope','tf','effect','P','FDR']].to_string(index=False))