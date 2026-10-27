from pathlib import Path
import pandas as pd,numpy as np,json,itertools
O=Path(r'C:\Users\31937\lung_feedback_revision_20261017_v19');B=Path(r'C:\Users\31937\lung_three_cohort_20261017_v16');V=Path(r'C:\Users\31937\lung_append_analysis_20261017_v18');T=O/'tables';r=pd.read_csv(V/'tables/all_associations.tsv',sep='\t');de=pd.read_csv(B/'data/tables/DE_all_models.tsv',sep='\t',low_memory=False);plan=json.loads((V/'protocol_frozen.json').read_text());dirs=pd.read_csv(V/'review_tables/第二队列方向_归一化口径对照.tsv',sep='\t');upd=pd.read_csv(V/'review_tables/相对v17的逐基因证据更新.tsv',sep='\t').set_index('gene')
positions={'TTC39C':'凋亡相关转录状态候选','MDFIC':'损伤转录关联备选','SKIL':'Treg活化/状态调控候选','IL1R2':'代谢与反馈调节候选','FOXP1':'状态重塑候选；损伤关联须解释组成','DDIT4':'文献驱动的应激/损伤假说候选','GSN':'备选；PROS1功能路线支持撤回','AHNAK':'T细胞表型/OXPHOS转录线索','ST6GALNAC3':'MitoROS响应转录线索'}
ph={'TTC39C':'TAct','MDFIC':'TAct','SKIL':'TAct','IL1R2':'Inhib','FOXP1':'Inhib','DDIT4':'Inhib','GSN':'Inhib','AHNAK':'Naive','ST6GALNAC3':'Naive'};dam={'TTC39C':'ApoExec','MDFIC':'ApoExec','SKIL':'ApoExec','IL1R2':'PyroExec','FOXP1':'PyroExec','DDIT4':'PyroExec','GSN':'ApoExec','AHNAK':'MitoROS','ST6GALNAC3':'MitoROS'}
def assoc(g,sc,y,co='Joint',model='JointStratum',var='Raw'):
 z=r[r.cohort.eq(co)&r.variant.eq(var)&r.scope.eq(sc)&r.model.eq(model)&r.x.eq('G_'+g)&r.y.eq(y)];return z.iloc[0] if len(z)==1 else None
rows=[]
for g,sc in plan['focus_pairs']:
 d=dict(gene=g,scope=sc,current_position=positions[g],not_final_innovation_or_experiment_investment_rank=True,phenotype_program=ph[g],injury_program=dam[g])
 for model in ['314_disease','joint_stratified']:
  z=de[(de.gene==g)&de.scope.eq(sc)&de.model.eq(model)];d[model+'_tested']=bool(len(z))
  for col in ['logFC','PValue','FDR']:d[model+'_'+col]=z[col].iloc[0] if len(z) else np.nan
 for label,y,co,model,var,ss in [('state','M_'+ph[g],'Joint','JointStratum','Raw',sc),('injury','M_'+dam[g],'Joint','JointStratum','Raw',sc),('injury_adjusted','M_'+dam[g],'Joint','JointStratumCondition','Raw',sc),('injury_original6','M_'+dam[g],'314','Pooled','Raw',sc),('time_day_adjusted','M_'+('MitoISR' if g=='TTC39C' else dam[g]),'274','TimeDayAdjusted','Raw','TLike')]:
  z=assoc(g,ss,y,co,model,var)
  for col in ['rho','P','FDR']:d[label+'_'+col]=z[col] if z is not None else np.nan
  d[label+'_source']='all_associations.tsv; '+co+'|'+var+'|'+ss+'|'+model+'|G_'+g+'|'+y
 for day in [7,56]:
  z=dirs[(dirs.gene==g)&dirs.scope.eq(sc)&dirs.day.eq(day)]
  for col in ['v14_cell_aggregate_difference','v16_TMM_logCPM_difference','v18_log1p_CP10k_difference']:d['D'+str(day)+'_'+col]=z[col].iloc[0] if len(z) else np.nan
 d['upstream_current']='PROS1-AXL适用性支持撤回；TYRO3/MERTK未通过覆盖，其他既有模型入口见R02_GSN_other_entries.tsv' if g=='GSN' else upd.loc[g,'limitation_or_contradiction']
 d['unresolved']=upd.loc[g,'limitation_or_contradiction'];d['B_current']='尚未达到完整机制支持标准；不是无选题价值';d['DE_source']='DE_all_models.tsv; see R01_candidate_DE_all_tested_records.tsv';rows.append(d)
d=pd.DataFrame(rows);d.to_csv(T/'R05_current_focus_conclusions.tsv',sep='\t',index=False)
y=pd.read_csv(T/'R05_current_all90.tsv',sep='\t');y=y.rename(columns={'v16_priority':'historical_v16_priority_not_current_rank','B_missing':'mechanism_support_gap','A_status':'historical_A_status'}).drop(columns=['B_status']);y['current_focus_position']=y.gene.map(positions).fillna('保留原A资格/待补定位；未纳入v18重点追加');y.to_csv(T/'R05_current_all90.tsv',sep='\t',index=False)
rr=pd.read_csv(B/'data/tables/candidate_complete_routes.tsv',sep='\t');rr[(rr.gene=='GSN')&(rr.scope=='ActCD4')&(rr.ligand!='PROS1')].drop_duplicates().to_csv(T/'R02_GSN_other_entries.tsv',sep='\t',index=False)
q=pd.read_csv(T/'TAM_receiver_group_qualification.tsv',sep='\t');q=q[(q.dataset=='GSE314943')&q.scope.eq('ActCD4')&q.gene.isin(['TYRO3','MERTK'])].copy();q['ligand']='PROS1';q['functional_source']='PMID25265470 mouse receptor experiments';q['current_qualification']=np.where(q.group_expression_pass,'expression compatible only','not supported by current expression coverage');q.to_csv(T/'R02_corrected_PROS1_alternatives.tsv',sep='\t',index=False)
# Explicit effect design matrices and case-label assignments, without rerunning any outcome.
s=pd.read_csv(V/'tables/sample_gene_and_program_scores.tsv',sep='\t');ef=pd.read_csv(V/'tables/programme_group_effects.tsv',sep='\t');rows=[];assign=[]
for (co,var,sc),zz in ef.groupby(['cohort','variant','scope']):
 ds=['GSE314943'] if co=='314' else ['GSE314943','GSE183682'] if co=='Joint' else ['GSE274372'];z=s[s.dataset.isin(ds)&s.variant.eq(var)&s.scope.eq(sc)].sort_values(['dataset','sample']).reset_index(drop=True)
 if co=='274_min10':z=z[z.cells>=10].reset_index(drop=True)
 id='|'.join([str(co),var,sc]);C={'Intercept':np.ones(len(z))}
 if str(co).startswith('274'):
  for col,vv in pd.get_dummies(z.day.astype(str),drop_first=True).items():C['day_'+col]=vv.to_numpy(float)
 else:
  lab=z.dataset.astype(str)+'_'+z.day.astype(str) if co=='Joint' else pd.Series(['all']*len(z))
  for col,vv in pd.get_dummies(lab,drop_first=True).items():C['stratum_'+col]=vv.to_numpy(float)
  case=z.group.eq('Silica').to_numpy(int);C['Silica']=case
  groups=[np.flatnonzero(lab.to_numpy()==l) for l in sorted(lab.unique())];choices=[list(itertools.combinations(g,int(case[g].sum()))) for g in groups]
  for j,parts in enumerate(itertools.product(*choices)):
   v=np.zeros(len(z),int)
   for pick in parts:v[list(pick)]=1
   for k,rr0 in z.iterrows():assign.append(dict(design_id=id,allocation=j,sample=rr0['sample'],Silica=int(v[k])))
 for i,rr0 in z.iterrows():
  for col,vals in C.items():rows.append(dict(design_id=id,sample=rr0['sample'],column=col,value=float(vals[i])))
pd.DataFrame(rows).to_csv(T/'R04_effect_design_matrices.tsv',sep='\t',index=False);pd.DataFrame(assign).to_csv(T/'R04_all_group_assignments.tsv',sep='\t',index=False)
cs=pd.read_csv(V/'tables/common_gene_sensitivity.tsv',sep='\t');cs['original_design_id']=cs.cohort.astype(str)+'|'+cs.original_variant+'|'+cs.scope+'|'+cs.model;cs[['family','original_design_id']].drop_duplicates().to_csv(T/'R04_common_gene_design_reference.tsv',sep='\t',index=False)
print('CURRENT_FOCUS',len(d),'GENES',d.gene.nunique(),'score_reconstruction_ready')