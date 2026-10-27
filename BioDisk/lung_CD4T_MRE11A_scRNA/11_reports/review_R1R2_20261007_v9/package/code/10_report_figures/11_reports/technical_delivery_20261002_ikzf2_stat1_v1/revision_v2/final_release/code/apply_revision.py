from pathlib import Path
import json,pandas as pd
P=Path('technical_delivery')
caps=json.loads((P/'figure_manifest.json').read_text(encoding='utf-8'))
for c in caps:
 if c['figure']=='F09_patient_fixed_axis_robustness':
  c['caption_zh']='沿用既报固定关联：MacAM IL18与CD4-IL7R IKZF2 DropOverlap；未择最小P。完整保留模型、seed和敏感性。All/Matched n=11：SourceDep名义P=0.00726–0.01169；BothEpi P=0.05256–0.07621、df降至3；全部FDR未通过。实心=名义P<0.05，空心=未达，叉号=无有效P；线为LOO范围，不是CI。受体限制：独立的RNA受体审计中，人All CD4-IL7R有16个合格样本，IL18R1/IL18RAP共检出中位0%，范围0–2.13%。RNA层受体支持有限，不等同蛋白缺失；该轴仍为关联候选，非已证实通信机制。'
  c['source'] = 'data/fixed_IL18_IKZF2_all_sensitivities.tsv; prior_review_repaired/scenic3_review/results/reviewed_sender_TF_models.tsv; source_tables/S3_CD4_receptor_codetection.tsv'
  c['caption_zh']+=' 接收范围是CD4-IL7R，不能与Core的TF–模块关联拼成同一已检验链条。'
 if c['figure']=='F07_human_module_associations':
  c['caption_zh']='人Core，原SourceDep模型，DropOverlap；All/BG1及三seed逐项展示12个冻结模块。颜色为原偏秩相关rho；N=有效名义P<0.05，F=原BH FDR<0.05；不重算。IKZF2–PyroExec六配置rho为−0.838至−0.476，ApoExec也全部负向；方向一致不等于每项显著，P/FDR逐项保留。不能写成IKZF2活性升高促进焦亡的正向因果链。该图接收范围为Core，IL18发送端关联为CD4-IL7R，不能跨范围拼成同一已检验链。MitoROS反映抗氧化应答，PyroExec/ApoExec不是蛋白裂解或细胞死亡实验证据。'
(P/'figure_manifest.json').write_text(json.dumps(caps,ensure_ascii=False,indent=2),encoding='utf-8')
r=pd.read_csv(P/'source_tables/S3_CD4_receptor_codetection.tsv',sep='\t')
r=r[(r.species=='Human')&(r['filter']=='All')&(r.scope=='CD4-IL7R')&(r.complex=='IL18_rec')&r.eligible20]
assert len(r)==16 and r.pct.median()==0 and abs(r.pct.max()-2.127659574468085)<1e-9
r.to_csv(P/'data/IL18_receptor_support_audit.tsv',sep='\t',index=False)
m=pd.read_csv(P/'source_tables/M2_human_source_models.tsv',sep='\t')
m=m[m.metric.str.startswith('score_')]
assert len(m)==504 and m.p_nominal.notna().sum()==396 and (m.p_nominal<.05).sum()==18
m.to_csv(P/'data/M2_all_module_group_models.tsv',sep='\t',index=False)
core=m[(m.scope=='Core')&(m['filter']=='All')&m.model.isin(['Source','SourceDep'])]
assert (core.p_nominal<.05).sum()==0 and len(core)==24
core.to_csv(P/'data/M2_Core_All_module_group_models.tsv',sep='\t',index=False)
a=pd.read_csv(P/'data/S2_TF_associations.tsv',sep='\t')
a=a[(a.species=='Human')&(a.scope=='Core')&(a.tf=='IKZF2')&(a.variant=='DropOverlap')&(a.model=='SourceDep')&a.feature.isin(['score_PyroExec','score_ApoExec'])]
assert len(a)==12 and (a.rho<0).all()
a.to_csv(P/'data/IKZF2_Core_death_module_direction_audit.tsv',sep='\t',index=False)
print('Revision evidence verified: receptor 16 samples; M2 504 models, 396 valid P, 18 nominal, Core All 0/12 in each model.')
