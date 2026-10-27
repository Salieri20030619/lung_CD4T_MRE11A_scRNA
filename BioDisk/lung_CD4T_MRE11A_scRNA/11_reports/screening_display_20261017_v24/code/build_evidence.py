from pathlib import Path
import pandas as pd,json,numpy as np
O=Path(__file__).parent;H=O.parent;B=H/'lung_three_cohort_20261017_v16';W=H/'lung_feedback_revision_20261017_v19';S=H/'lung_dual_benchmark_20261017_v17';T=O/'tables';T.mkdir(exist_ok=True)
r=lambda p:pd.read_csv(p,sep='\t',low_memory=False)
flow=r(W/'tables/R01_selection_flow.tsv');flow.to_csv(T/'screening_flow.tsv',sep='\t',index=False)
base=r(B/'final_tables/四轴候选总表_55基因.tsv');supp=r(B/'final_tables/历史候选逐项回溯.tsv');all90=r(W/'tables/R05_current_all90.tsv');x=r(S/'assessment/逐状态双基准评估.tsv');fn=r(B/'data/tables/function_membership_extended.tsv');impact=r(W/'tables/R02_candidate_entry_impact.tsv')
assert len(set(base.gene))==55 and len(set(supp.gene))==42 and len(set(base.gene)&set(supp.gene))==7 and set(base.gene)|set(supp.gene)==set(all90.gene)
stages={'F05':x};stages['F06']=x[x.gene.isin(fn.gene)];z=stages['F06'].merge(impact[['gene','scope','family','has_other_entry']],on=['gene','scope','family'],validate='one_to_one');stages['F07']=z[z.has_other_entry];stages['F08']=stages['F07'][~stages['F07'].gene.isin(['GSDMD','S100A8'])];stages['F09']=stages['F08'][~stages['F08'].gene.isin(['GZMB','TNFSF10','IL10'])]
for key,z in stages.items():
 q=flow[flow.step==key].iloc[0];assert z.gene.nunique()==q.retained_genes and len(z)==q.retained_records
assert set(stages['F09'].gene)==set(all90[all90.A_status_current.str.startswith('A1')].gene)
status={'A1':'符合计算筛选条件','A2':'需补充证据','A3':'未达到当前名义差异条件'};out=[]
for a in all90.itertuples():
 row={'gene':a.gene,'in_three_axis_initial55':a.gene in set(base.gene),'in_supplemental42':a.gene in set(supp.gene),'current_status':status[a.A_status_current[:2]]}
 for key,z in stages.items():row[key+'_retained']=a.gene in set(z.gene)
 out.append(row)
pd.DataFrame(out).to_csv(T/'candidate90_stage_membership.tsv',sep='\t',index=False)
rows=[];prev=set(all90.gene)
for key,z in stages.items():
 now=set(z.gene);lost=sorted(prev-now);rows.append({'step':key,'retained_genes':len(now),'retained_records':len(z),'moved_out_at_step':len(lost),'genes':' ; '.join(lost)});prev=now
pd.DataFrame(rows).to_csv(T/'stage_transition_members.tsv',sep='\t',index=False)
focus=['IL1R2','DDIT4','SKIL','FOXP1','GSN','AHNAK','TTC39C','MDFIC','ST6GALNAC3'];assert sum(all90.set_index('gene').loc[focus].A_status_current.str.startswith('A1'))==5
(O/'screening_QA.json').write_text(json.dumps({'initial55':55,'supplemental42':42,'overlap':7,'union':90,'sequential_genes':[89,57,56,54,52],'sequential_records':[344,172,168,164,160],'status_counts':[52,37,1],'focus9_supported':5,'focus9_pending':4,'screening_results_changed':False},indent=2),encoding='utf-8');print(pd.DataFrame(rows).to_string(index=False))