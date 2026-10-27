from pathlib import Path
import pandas as pd,numpy as np,json,sys
O=Path(sys.argv[1])
pt=pd.read_csv(O/'03A_state_pseudotime/genomewide_trajectory_tests.tsv.gz',sep='\t');tm=pd.read_csv(O/'03B_real_time/Omnibus.tsv',sep='\t');tm['canonical']=tm.symbol.str.upper()
print('TIME_DUPLICATES',int(tm.canonical.duplicated().sum()))
time={}
for name in ['D7_D0','D28_D0','D56_D0']:
 t=pd.read_csv(O/'03B_real_time'/f'{name}.tsv',sep='\t');t['canonical']=t.symbol.str.upper();time[name]=t.set_index('canonical')
rows=[]
for state in ['Naive','ActCD4']:
 de=pd.read_csv(O/'04_candidate_intersection'/f'{state}_discovery_DE_frozen.tsv',sep='\t');a=de.merge(pt[(pt.state==state)&(pt.run=='K15_R1')],left_on='canonical',right_on='feature',how='left',suffixes=('_DE','_PT'))
 a=a.rename(columns={'PValue':'DE_P','FDR':'DE_BH','logFC':'DE_logFC','P':'PT_P','BH':'PT_BH'})
 a=a.merge(tm[['canonical','PValue','FDR']],on='canonical',how='left').rename(columns={'PValue':'time_omnibus_P','FDR':'time_omnibus_BH'})
 a['time_nominal_direction_match']=False;a['matched_days']=''
 for name,t in time.items():
  for c in ['logFC','PValue','FDR']:a[name+'_'+c]=a.canonical.map(t[c])
  ok=(a[name+'_PValue']<.05)&(np.sign(a[name+'_logFC'])==np.sign(a.DE_logFC))
  a.loc[ok,'time_nominal_direction_match']=True;a.loc[ok,'matched_days']+=name+';'
 a['pass_DE']=a.DE_P<.05;a['pass_PT']=a.PT_P<.05;a['PT_direction_match']=np.sign(a.mean_slope)==np.sign(a.DE_logFC);a['pass_time_omnibus']=a.time_omnibus_P<.05
 a['selected']=a.pass_DE&a.pass_PT&a.PT_direction_match&a.pass_time_omnibus&a.time_nominal_direction_match
 a['state']=state
 a.to_csv(O/'04_candidate_intersection'/f'{state}_all_gene_decisions.tsv.gz',sep='\t',index=False,na_rep='NA')
 stages=[('tested_disease',np.ones(len(a),bool)),('DE_P05',a.pass_DE),('DE_and_PT_P05',a.pass_DE&a.pass_PT),('DE_and_PT_direction',a.pass_DE&a.pass_PT&a.PT_direction_match),('plus_time_omnibus',a.pass_DE&a.pass_PT&a.PT_direction_match&a.pass_time_omnibus),('plus_time_direction',a.selected)]
 print(state,[(s,int(v.sum())) for s,v in stages]);print(a[a.selected][['canonical','DE_logFC','DE_P','DE_BH','mean_slope','PT_P','PT_BH','time_omnibus_P','matched_days']].to_string(index=False))
 for stage,mask in stages:
  rows.append(dict(state=state,stage=stage,genes=int(mask.sum())))
  a.loc[mask,['canonical']].to_csv(O/'04_candidate_intersection'/f'{state}_{stage}_members.tsv',sep='\t',index=False)
 a[a.selected].to_csv(O/'04_candidate_intersection'/f'{state}_primary_intersection.tsv',sep='\t',index=False)
pd.DataFrame(rows).to_csv(O/'04_candidate_intersection/stepwise_gene_counts.tsv',sep='\t',index=False)