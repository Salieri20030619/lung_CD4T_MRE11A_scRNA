from pathlib import Path
import pandas as pd,json
O=Path(r'C:\Users\31937\lung_append_analysis_20261017_v18');B=Path(r'C:\Users\31937\lung_three_cohort_20261017_v16');V=Path(r'C:\Users\31937\lung_review_20261017_v14\revision_v14');p=json.loads((O/'protocol_frozen.json').read_text());v16=pd.read_csv(B/'data/tables/nominal_DE_axis_membership.tsv',sep='\t');v14=pd.read_csv(V/'supplement/S07_additional_candidate_evidence.tsv',sep='\t');v18=pd.read_csv(O/'tables/183682_timepoint_descriptive.tsv',sep='\t');rows=[]
for gene,scope in p['focus_pairs']:
 a=v16[(v16.gene==gene)&(v16.scope==scope)];b=v14[(v14.gene==gene)&(v14.mouse_state==scope)]
 for day in [7,56]:
  c=v18[(v18.variant=='Raw')&(v18.scope==scope)&(v18.day==day)&(v18.variable=='G_'+gene)];rows.append(dict(gene=gene,scope=scope,day=day,v14_cell_aggregate_difference=b['external_D'+str(day)+'_difference'].iloc[0] if len(b) else None,v16_TMM_logCPM_difference=a['183_D'+str(day)+'_descriptive_delta'].iloc[0] if len(a) else None,v18_log1p_CP10k_difference=c.delta.iloc[0] if len(c) else None,interpretation='各口径单位/聚合不同，比较方向而非数值大小；n1+1不作独立推断'))
d=pd.DataFrame(rows);d.to_csv(O/'review_tables/第二队列方向_归一化口径对照.tsv',sep='\t',index=False);print(d[d.gene=='TTC39C'].to_string(index=False))