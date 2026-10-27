from pathlib import Path
import pandas as pd,numpy as np,json
O=Path(r'C:\Users\31937\lung_append_analysis_20261017_v18');r=pd.read_csv(O/'tables/all_associations.tsv',sep='\t');p=json.loads((O/'protocol_frozen.json').read_text());out=O/'review_tables';out.mkdir(exist_ok=True)
def get(cohort,var,model,scope,x,y,kind=None,sender=None):
 z=r[r.cohort.eq(cohort)&r.variant.eq(var)&r.model.eq(model)&r.scope.eq(scope)&r.x.eq(x)&r.y.eq(y)]
 if kind is not None:z=z[z.kind==kind]
 if sender is not None:z=z[z.sender==sender]
 assert len(z)<=1,(cohort,var,model,scope,x,y,len(z))
 return z.iloc[0] if len(z) else None
primary=r[(r.cohort=='Joint')&(r.variant=='Raw')&(r.model=='JointStratum')&r.kind.str.startswith('candidate_')];rows=[]
for a in primary.itertuples():
 d=dict(gene=a.gene,scope=a.scope,program=a.program,kind=a.kind,primary_rho=a.rho,primary_P=a.P,primary_FDR=a.FDR,primary_loo_same_sign=a.loo_same_sign)
 for label,cohort,var,model,sc in [('condition','Joint','Raw','JointStratumCondition',a.scope),('original314','314','Raw','Pooled',a.scope),('BG1','314','BG1','Pooled',a.scope),('fixedComposition','Joint','FixedComp','JointStratum',a.scope),('time274','274','Raw','TimePooled','TLike'),('time274_dayAdjusted','274','Raw','TimeDayAdjusted','TLike')]:
  z=get(cohort,var,model,sc,a.x,a.y)
  for k in ['rho','P','FDR']:d[label+'_'+k]=z[k] if z is not None else np.nan
 rows.append(d)
pd.DataFrame(rows).to_csv(out/'候选表型损伤_主分析与敏感性对照.tsv',sep='\t',index=False)
rows=[]
for g,sc,l,ts in p['routes']:
 for sender in ['AllMac','AMlike']:
  for tf in (ts.split(';') if ts not in ['none','decoy'] else [ts]):
   d=dict(gene=g,scope=sc,ligand=l,sender=sender,proximal_TF=tf,proximal_interpretation='不设激活TF：诱饵受体' if tf=='decoy' else '本轮未预设特异近端TF' if tf=='none' else '共享转录响应兼容性，非唯一受体路径')
   for model in ['Pooled','Condition']:
    for label,x,y,kind,send in [('ligand_candidate','L_'+l,'G_'+g,'ligand_candidate',sender),('ligand_TF','L_'+l,'TF_'+tf,'ligand_proximalTF',sender),('TF_candidate','TF_'+tf,'G_'+g,'proximalTF_candidate','none')]:
     z=get('314','Raw',model,sc,x,y,kind,send)
     for k in ['rho','P','FDR']:d[model+'_'+label+'_'+k]=z[k] if z is not None else np.nan
   rows.append(d)
pd.DataFrame(rows).to_csv(out/'上游连接逐段证据_不跨来源拼接.tsv',sep='\t',index=False)
updates=[
['IL1R2','Treg','新增代谢/抑制表型相关支持','ENO1 RNA联合rho=0.912，P=0.000694，FDR=0.0285；控制组别rho=0.610，P=0.111','RNA正相关不证明抑制ENO1酶活；尚缺IL1A/B驱动候选的实际作用','保留代谢反馈方向；不要以RNA相关证明抗焦亡机制'],
['DDIT4','ActCD4','表型支持增加，焦亡连接未增强','Inhib/Naive/Cycle有联合名义关联；PyroExec rho=0.126，P=0.751；TNF-DDIT4 rho=-0.257，P=0.658','MitoISR调整后rho=0.888，置换P=0.0556；RELB-DDIT4仅一段名义支持','保留文献支撑的应激候选，不把TNF-DDIT4焦亡链作为已支持主线'],
['SKIL','Treg','Treg活化联系增强，损伤主张需收窄','TAct联合rho=0.965，P=0.00208；控制组别rho=0.897，P=0.0278','ApoExec未达名义关联；SKIL-p53 AUC pooled P=0.919；CCL2关联来自AMlike且调整后减弱','优先定位活化/状态调节问题，损伤连接另行验证'],
['FOXP1','Core/Treg','状态重塑支持较强，Core焦亡关联组成敏感','Core-PyroExec联合rho=-0.822，P=0.0219；固定组成后rho=-0.266，P=0.494','独立TLike时间有同向线索，但非纯CD4验证；SEMA3E/CD200来源竞争仍在','保留T细胞状态/损伤易感性问题，不用原Core关联证明内源焦亡调控'],
['GSN','ActCD4','维持备选，新增损伤连接不足','Inhib联合名义关联；ApoExec rho=-0.459，P=0.264；PROS1关联仅AMlike pooled名义支持','PROS1非巨噬特异，来源/蛋白形式和对照异质性未消除','保留备选，不升级机制主线'],
['AHNAK','Naive','氧化磷酸化转录联系，非损伤证据','OXPHOS联合rho=0.811，P=0.0417；NFAT活动关联未达名义支持','OXPHOS RNA程序不是呼吸功能或线粒体损伤实测','保留T细胞表型候选，损伤缺口仍需补'],
['TTC39C','ActCD4','新增较强凋亡转录关联，值得重新关注','ApoExec联合rho=0.943，P=0.00347，FDR=0.0556；控制组别rho=0.922，P=0.0278；时间TLike-MitoISR日调整P=0.0185','第二队列D7方向对归一化/聚合口径敏感，D56反向；原TTC39C表型与上游链未闭合','提升为损伤转录关联的第二梯队备选，提交独立评审'],
['MDFIC','ActCD4','新增较强凋亡转录关联','ApoExec联合rho=0.947，P=0.00278，FDR=0.0556；控制组别rho=0.950，P=0.0833','效应方向强但条件置换分辨率低；第二队列疾病方向相反，上游未闭合','提升为损伤转录关联的第二梯队备选，不要求每个敏感性P显著'],
['ST6GALNAC3','ActCD4','新增ROS响应转录关联，仍可能受疾病状态驱动','MitoROS联合rho=0.872，P=0.00972；控制组别rho=0.569，P=0.222','MitoROS为应激响应，不是ROS实测；直接功能依据不足','继续保留条件待补，不用泛应激关联升级损伤机制']]
pd.DataFrame(updates,columns=['gene','scope','new_conclusion','support','limitation_or_contradiction','updated_position']).to_csv(out/'相对v17的逐基因证据更新.tsv',sep='\t',index=False)
summary={'focus_genes':9,'gene_state_pairs':10,'upstream_routes':18,'all_associations':len(r),'planned_association_families':r.family.nunique(),'joint_primary_candidate_tests':len(primary),'joint_primary_nominal':int((primary.P<.05).sum()),'joint_primary_FDR':int((primary.FDR<.05).sum()),'new_secondary_candidates':['TTC39C','MDFIC'],'key_risk':'FOXP1 Core pyroptosis-score association attenuates under fixed composition; neither original five nor additions establish complete macrophage->T phenotype->injury causality','delivery':'append after unchanged v17 baseline, full negative and sensitivity tables included'};(out/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(summary,ensure_ascii=False,indent=2))