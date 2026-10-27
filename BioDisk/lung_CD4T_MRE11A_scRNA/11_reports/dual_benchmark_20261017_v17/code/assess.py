from pathlib import Path
import pandas as pd,numpy as np,json,re
B=Path(r'C:\Users\31937\lung_three_cohort_20261017_v16');O=Path(r'C:\Users\31937\lung_dual_benchmark_20261017_v17');V14=Path(r'C:\Users\31937\lung_review_20261017_v14\revision_v14');T=B/'data/tables';out=O/'assessment';out.mkdir(exist_ok=True)
rd=lambda n:pd.read_csv(T/n,sep='\t',low_memory=False)
genes=json.loads((O/'genes.json').read_text());a=rd('nominal_DE_axis_membership.tsv');a=a[a.gene.isin(genes)].copy();fn=rd('function_membership_extended.tsv');q=pd.read_csv(O/'tables/90_gene_sample_expression.tsv',sep='\t');net=pd.read_csv(O/'tables/network_links_no_rank_exclusion.tsv.gz',sep='\t');lr=rd('expressed_LR_with_decoy_literature.tsv');old=pd.read_csv(B/'final_tables/四轴候选总表_55基因.tsv',sep='\t').set_index('gene')
manual={
'AHNAK':('一般凋亡/明确T细胞表型','CD4钙流、NFAT和激活有直接小鼠研究；人癌细胞53BP1-p53/AnnexinV研究未直接证明T细胞线粒体损伤或焦亡','https://pmc.ncbi.nlm.nih.gov/articles/PMC2350190/;https://pmc.ncbi.nlm.nih.gov/articles/PMC8221568/'),
'TTC39C':('一般凋亡','A549/H1299敲低增加AnnexinV；能量/p53通路分析不是线粒体损伤或焦亡的直接终点','https://pmc.ncbi.nlm.nih.gov/articles/PMC9615393/'),
'MDFIC':('跨物种凋亡启动线索','鸡成肌细胞MDFIC干预及cleaved-caspase9结果；需补哺乳T细胞作用与线粒体终点，不作为已确证T细胞机制','https://pmc.ncbi.nlm.nih.gov/articles/PMC6678156/'),
'ST6GALNAC3':('目标功能尚缺证据','当前未核实到目标T细胞线粒体/焦亡直接证据；同名lncRNA与蛋白编码基因不能混用','https://www.ncbi.nlm.nih.gov/nuccore/NM_152996'),
'GZMB':('作用细胞角色待澄清','杀伤细胞GZMB可切割靶细胞GSDME；不等于表达GZMB的T细胞自身焦亡','https://www.nature.com/articles/s41586-020-2071-9'),
'TNFSF10':('作用细胞角色待澄清','作为死亡配体的旁分泌/自分泌作用需区分；T细胞表达不能直接说明自身死亡','v16/function_membership_extended.tsv'),
'IL10':('作用细胞角色待澄清','分泌性抗炎/死亡调控不能自动证明IL10表达T细胞自身损伤','v16/function_membership_extended.tsv')}
identity_flags={'S100A8':'髓系标记/背景来源待核实','GSDMD':'v14低UMI及严格过滤丢失对照样本，内源身份未解决','TFF1':'上皮来源待核实','SPP1':'髓系/其他来源待核实','LCN2':'髓系/上皮来源待核实','ESM1':'内皮来源待核实'}
rows=[]
for r in a.itertuples():
 ds='GSE274372' if r.model=='274_time' else 'GSE314943';qq=q[(q.gene==r.gene)&(q.scope==r.scope)&(q.dataset==ds)];express=bool((qq.positive>0).any());ident=identity_flags.get(r.gene,'免疫球蛋白/浆细胞来源待核实' if re.match(r'^(IGH|IGK|IGL|JCHAIN)',r.gene) else '')
 gg=fn[fn.gene==r.gene];function=len(gg)>0;role=r.gene in ['GZMB','TNFSF10','IL10'];n=net[(net.gene==r.gene)&(net.scope==r.scope)&(net.family==r.family)] if r.contrast!='time_omnibus' else net[(net.gene==r.gene)&(net.scope==r.scope)&net.family.str.startswith('274_time__')]
 ll=lr[(lr.dataset==ds)&(lr.scope==r.scope)&(lr.sender=='AllMac')&lr.receptor_parts.str.split(';').map(lambda z:r.gene in z)]
 # Nonzero matrix weight = prior coverage; it is explicitly not an inferred treatment effect.
 hasroute=len(n)>0 or len(ll)>0;decoy_only=len(ll)>0 and ll.annotation.str.contains('decoy').all() and len(n)==0
 gaps=[]
 if not express:gaps.append('无当前范围表达支持')
 if ident:gaps.append(ident)
 if not function:gaps.append(manual.get(r.gene,('目标功能尚缺证据',))[0])
 if role:gaps.append('自身损伤与作用其他细胞尚未区分')
 if not hasroute:gaps.append('本次资源未找到匹配上游入口，属于待补而非生物学否定')
 if decoy_only:gaps.append('仅诱饵结合，缺候选内源功能耦合')
 standard='A1_计算支持入池' if not gaps else 'A2_条件保留待补'
 rows.append(dict(gene=r.gene,scope=r.scope,model=r.model,family=r.family,contrast=r.contrast,logFC=r.logFC,P=r.PValue,FDR=r.FDR,A_status=standard,A_gaps='；'.join(gaps),expression_samples=int((qq.positive>0).sum()),expression_available_samples=len(qq),positive_cells=int(qq.positive.sum()),multiUMI_cells=int(qq.multiUMI.sum()),Tmarker_fraction=float(qq.Tmarker_positive.sum()/qq.positive.sum()) if qq.positive.sum() else np.nan,identity_note=ident or '冻结T细胞身份及共标记支持；不等于蛋白来源确证',function_status=manual[r.gene][0] if r.gene in manual else '特定通路/文献调控依据' if function else '目标功能尚缺证据',function_branches=';'.join(sorted(gg.mechanism.unique())),direct_receptor=bool(len(ll)),network_prior_links=len(n),best_target_rank=float(n.target_rank.min()) if len(n) else np.nan,best_activity_rank=float(n.activity_rank.min()) if len(n) else np.nan,within_v16_reference=bool(n.v16_reference_rank.any()) if len(n) else False,network_grade='直接受体/结合入口（具体功能另评）' if len(ll) else '较靠前的模型支持，非因果' if len(n) and n.v16_reference_rank.any() else '弱模型先验，保留不因排名淘汰' if len(n) else '入口待补',B_status='B_未达到机制主候选',B_gap='本项目T细胞自身损伤/死亡执行终点及同源同状态连续连接未建立；RNA变化/模块评分为代理证据',temporal_causal_order='未建立',raw_P_gate=True,FDR_gate=False,rank_gate=False))
x=pd.DataFrame(rows);x.to_csv(out/'逐状态双基准评估.tsv',sep='\t',index=False)
allrows=[]
for gene in genes:
 z=x[x.gene==gene];good=z[z.A_status=='A1_计算支持入池'];rep=(good if len(good) else z);rep=rep.sort_values(['model','P']) if len(rep) else rep
 status='A1_计算支持入池' if len(good) else 'A2_条件保留待补' if len(z) else 'A3_本轮未达名义差异或不可评估'
 scopes=';'.join(sorted(good.scope.unique()));note='；'.join(sorted(set(z.A_gaps.dropna())-{''})) if not len(good) else '状态限定：'+scopes+'；其余状态见明细'
 allrows.append(dict(gene=gene,v16_in_55=gene in old.index,A_status=status,A_supported_scopes=scopes,A_notes=note,reference_best_P=float(z.P.min()) if len(z) else np.nan,any_old_reference_route=bool(z.within_v16_reference.any()) if len(z) else False,any_network_prior=bool((z.network_prior_links>0).any()) if len(z) else False,function_or_role_note=manual[gene][1] if gene in manual else '',manual_source=manual[gene][2] if gene in manual else '',v16_priority=old.loc[gene,'review_priority'] if gene in old.index else '历史补评',B_status='B_未达到机制主候选',B_missing='T细胞自身线粒体损伤/焦亡表型尚未被充分确证；巨噬到表型及损伤的连续作用连接尚缺',novelty_hard_gate=False))
y=pd.DataFrame(allrows);y.to_csv(out/'90基因双基准总表.tsv',sep='\t',index=False)
# Retain the complete pre-existing phenotype proxies and proximal tests, including negative and NA.
for folder,names in [('supplement',['S04_death_program_comparisons.tsv','S05_mito_autophagy_comparisons.tsv','S07_additional_candidate_evidence.tsv']),('supplement_v14',['S23_targeted_proximal_associations.tsv','S24_complex_codetection_raw.tsv','S27_control_overlap_associations.tsv'])]:
 for n in names:(out/('v14_'+n)).write_bytes((V14/folder/n).read_bytes())
pd.DataFrame([dict(gene=g,evidence_class=v[0],interpretation=v[1],source=v[2]) for g,v in manual.items()]).to_csv(out/'新增文献及细胞作用角色复核.tsv',sep='\t',index=False)
assert y.gene.is_unique and len(y)==90;assert (x.P<.05).all();assert not x.FDR_gate.any() and not x.rank_gate.any();assert set(x.loc[x.A_status=='A1_计算支持入池','A_gaps'])=={''}
counts=y.A_status.value_counts().to_dict();j={'genes':90,'scope_contrasts':len(x),'A_counts':counts,'B_qualified':0,'A_new_vs_v16':y.loc[(y.A_status=='A1_计算支持入池')&~y.v16_in_55,'gene'].tolist(),'v16_moved_to_pending':y.loc[y.v16_in_55&(y.A_status!='A1_计算支持入池'),'gene'].tolist(),'interpretation':'A1 is computational candidate eligibility only, not proof of macrophage causality or T-cell death. A2 retains missing/weak evidence. B=0 means evidence not yet adequate, not absence of biology.'};(out/'summary.json').write_text(json.dumps(j,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(j,ensure_ascii=False,indent=2));print(y[y.gene.isin(['AHNAK','TTC39C','MDFIC','ST6GALNAC3','IL1R2','DDIT4','SKIL','FOXP1','GSN','GSDMD','IL10','GZMB','TNFSF10'])][['gene','A_status','A_supported_scopes','A_notes']].to_string(index=False))