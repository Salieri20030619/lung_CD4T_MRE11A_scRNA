from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,datetime
R=Path('C:/codex_test');T=R/'tmp/silicosis_v11';O=R/'output/silicosis_v11';S=R/'output/silicosis_v10/source_tables/03_receiver'
for x in ['supplement','provenance','code']:(O/x).mkdir(parents=True,exist_ok=True)
SHARE='https://chatgpt.com/share/6ad1bfd1-2c50-83ea-9c27-6df48e857982'
CAND=['IL1R2','PENK','TNFRSF9','LAG3','TTC39C','FOXP1','IKZF2','STAT1','AHNAK','MDFIC','ST6GALNAC3']
REF={}
def ref(id,title,year,url,context,finding):REF[id]=dict(id=id,title=title,year=year,url=url,context=context,finding=finding)
ref('L01','An atlas of the aging lung mapped by single cell transcriptomics and deep tissue proteomics',2019,'https://www.nature.com/articles/s41467-019-08831-9','小鼠肺单细胞图谱','肺图谱区分B细胞、浆细胞与T细胞；用于解释免疫球蛋白及JCHAIN的来源核查需要。')
ref('L02','Reciprocal expression of trefoil factor-1 and thyroid transcription factor-1 in lung adenocarcinomas',2020,'https://pmc.ncbi.nlm.nih.gov/articles/PMC7293082/','人肺腺癌上皮','TFF1在部分具有胃肠分化特征的肺腺癌中表达。该证据不等于正常CD4表达或矽肺CD4功能证据。')
ref('L03','The short-chain fatty acid acetate coordinates with CD30 to modulate T-cell survival',2023,'https://pmc.ncbi.nlm.nih.gov/articles/PMC10398883/','小鼠T细胞','研究CD30/TNFRSF8、HDAC6、微管乙酰化与T细胞生存/凋亡关系。')
ref('L04','Targeted Disruption of the Transition Protein 2 Gene Affects Sperm Chromatin Structure and Reduces Fertility in Mice',2001,'https://pmc.ncbi.nlm.nih.gov/articles/PMC99899/','小鼠生精细胞','Tnp2缺失影响精子染色质凝集及生育力。')
ref('L05','Prohormone convertase 1/3 deficiency causes obesity due to impaired proinsulin processing',2022,'https://pmc.ncbi.nlm.nih.gov/articles/PMC9376086/','神经内分泌及胰岛细胞','PCSK1/PC1/3参与前激素加工，缺失影响前胰岛素处理。')
ref('L06','Defective development of secretory neurones in the hypothalamus of Arnt2-knockout mice',2001,'https://doi.org/10.1046/j.1365-2443.2001.00421.x','小鼠下丘脑神经元','Arnt2敲除导致分泌神经元发育异常。')
ref('L07','Knockout of STriatal Enriched Protein Tyrosine Phosphatase in Mice Results in Increased ERK1/2 Phosphorylation',2008,'https://pmc.ncbi.nlm.nih.gov/articles/PMC2706508/','小鼠脑组织','PTPN5/STEP敲除与脑内ERK1/2磷酸化升高有关。')
ref('L08','ESM-1 is a novel human endothelial cell-specific molecule expressed in lung and regulated by cytokines',1996,'https://pubmed.ncbi.nlm.nih.gov/8702785/','人内皮细胞及肺组织','从HUVEC克隆ESM1，并验证其表达和细胞因子调节。')
ref('L09','Role of IL-10-producing regulatory B cells in modulating T-helper cell immune responses during silica-induced lung inflammation and fibrosis',2016,'https://pmc.ncbi.nlm.nih.gov/articles/PMC4926212/','硅尘小鼠、调节B细胞及Th免疫反应','已有硅尘模型中IL10/B10调节Th反应及炎症纤维化的实验研究。')
ref('L10','Maintenance of sarcomeric integrity in adult muscle cells crucially depends on Z-disc anchored titin',2020,'https://pmc.ncbi.nlm.nih.gov/articles/PMC7478974/','成年小鼠骨骼肌','Ttn条件失活破坏肌节完整性并影响肌肉功能。')
ref('L11','SPP1+ macrophages contribute to silica-induced pulmonary inflammation and fibrosis via the ERK/JNK signaling pathway',2026,'https://academic.oup.com/ajrcmb/article/74/7/936/8538001','矽肺、巨噬细胞','单细胞/空间及干预研究显示SPP1阳性巨噬细胞与硅尘炎症纤维化有关，抑制SPP1减轻损伤。')
ref('L12','Nucleoside diphosphate kinase B knock-out mice have impaired activation of the K+ channel KCa3.1, resulting in defective T cell activation',2010,'https://pubmed.ncbi.nlm.nih.gov/20884616/','小鼠及人CD4 T细胞','NME2/NDPK-B对KCa3.1活性和Th1/Th2功能已有遗传学支持。')
ref('L13','CCR8+FOXp3+ Treg cells as master drivers of immune regulation',2017,'https://pubmed.ncbi.nlm.nih.gov/28533380/','人/小鼠Treg与自身免疫模型','CCL1-CCR8可增强Treg免疫抑制功能；不等同于矽肺ActCD4机制。')
ref('L14','Alveolar fibroblast lineage orchestrates lung inflammation and fibrosis',2024,'https://doi.org/10.1038/s41586-024-07660-1','肺损伤模型、肺泡成纤维细胞','炎性成纤维细胞表达Lcn2，提示需区分结构细胞与T细胞来源。')
ref('L15','Selective attraction of monocytes and T lymphocytes of the memory phenotype by cytokine RANTES',1990,'https://pubmed.ncbi.nlm.nih.gov/1699135/','人单核细胞与记忆表型CD4 T细胞','RANTES/CCL5对单核细胞和记忆表型CD4 T细胞的趋化作用已有直接实验。')
ref('L16','Gelsolin expression is necessary for the development of modelled pulmonary inflammation and fibrosis',2009,'https://thorax.bmj.com/content/64/6/467','肺纤维化组织和博来霉素小鼠','GSN与肺炎症、细胞迁移及上皮凋亡/纤维化已有模型研究。')
ref('L17','Transcriptional repressor ZEB2 promotes terminal differentiation of CD8+ effector and memory T cell populations during infection',2015,'https://pmc.ncbi.nlm.nih.gov/articles/PMC4647262/','感染模型、CD8 T细胞','ZEB2参与效应CD8 T细胞终末分化；该结果不能直接外推为矽肺CD4机制。')
ref('L18','Neurotrophin-3 attenuates human peripheral blood T cell and monocyte activation status and cytokine production post stroke',2022,'https://doi.org/10.1016/j.expneurol.2021.113901','卒中患者/对照PBMC离体实验','检测TrkC并观察NT3对T细胞活化及细胞因子的影响；特定受体因果机制仍需区分。')
ref('D01','NCBI Gene: Gm2663, predicted gene 2663',2026,'https://www.ncbi.nlm.nih.gov/gene/100040208','小鼠基因注释','蛋白编码基因，RefSeq为INFERRED，优选蛋白名为trypsinogen 4；功能多为预测。')
ref('D02','MGI: D830030K20Rik, RIKEN cDNA D830030K20 gene',2026,'https://www.informatics.jax.org/marker/MGI%3A2443830','小鼠基因注释','数据库记录支持基因身份；当前记录不足以支持矽肺CD4机制成熟这一判断。')

ANN={}
def ann(g,r,reason,status,check,hist='未找到逐基因历史排除决定；以下为本次补充评估'):
 ANN[g]=dict(ref=r,reason=reason,status=status,next_check=check,history=hist)
for g in ['IGHG2C','IGHG1','IGHG2B','IGLV1','IGKC','IGHA','IGLC1','JCHAIN']:
 h='2026-09-14明确列为需核查细胞来源；不是已确认污染或最终排除' if g in ['IGKC','IGHG2B','IGHG2C','JCHAIN'] else '历史提出免疫球蛋白来源核查；本基因的具体处置为本次补充评估'
 ann(g,'L01','免疫球蛋白/分泌抗体相关信号须先确认是否来自真实CD4表达。不能因log2FC大直接作为CD4内机制靶点，也不能仅凭检出率高认定污染。','先核实细胞来源','逐样本共检出、UMI负担、双细胞/环境RNA敏感性及纯化T细胞验证',h)
ann('TFF1','L02','已知肺上皮肿瘤研究不直接支持矽肺CD4功能。高效应需要先核查本队列的细胞来源和表达分布。','先核实细胞来源','分样本表达、上皮标记共检出及纯化CD4验证')
ann('TNFRSF8','L03','CD30与T细胞生存/凋亡已有机制基础，单纯重复该关系的创新空间有限；但与本主题相关，尚无充分依据将其排除。','保留比较候选','比较CD4亚群内表达、STAT1/死亡关联及矽肺特异上下游关系')
ann('TNP2','L04','原始功能证据集中于生精细胞染色质，不能用该背景解释CD4差异。应先验证检测与比对，再决定机制实验。','先核实细胞来源','基因唯一比对、样本分布和独立CD4表达验证')
ann('PCSK1','L05','既有证据主要为前激素加工。暂缺将本队列PCSK1差异连接到CD4死亡或巨噬细胞信号的证据。','先核实表达与作用位置','确认CD4表达后再检验相关底物和功能联系')
ann('ARNT2','L06','神经发育机制较明确，但不足以判定矽肺CD4中的作用；不能仅以其他组织研究成熟作为淘汰依据。','先核实表达与作用位置','CD4来源验证及亚群内调控证据')
ann('PTPN5','L07','脑内STEP-ERK机制已有研究，本项目尚缺CD4来源和对应功能连接；先验证而非直接进入主轴。','先核实表达与作用位置','唯一比对、纯化CD4表达及ERK相关功能')
ann('ESM1','L08','内皮表达与细胞因子调节已有直接实验。CD4分群中的强差异需排查内皮/上皮来源贡献，暂不作为CD4内驱动。','先核实细胞来源','与内皮及上皮标记共检出、独立CD4验证')
ann('IL10','L09','同疾病的免疫调节研究已有覆盖，泛化的IL10抗炎叙事新颖性有限；仍可作为状态参照或细胞特异机制候选。','保留比较候选','区分B10/Treg/ActCD4来源，并验证与当前调控轴的关系')
ann('TTN','L10','主要遗传机制证据来自肌节结构。CD4中的差异首先需要基因覆盖、比对和来源核查，不能直接作为死亡执行分子。','先核实细胞来源','检查比对覆盖、异源标记、样本分布及纯化细胞表达')
ann('SPP1','L11','矽肺巨噬细胞SPP1已有同疾病、同来源及干预研究；作为新主靶点需证明额外关系。ActCD4中的SPP1还需独立核实来源。','成熟机制参照并核实来源','先区分巨噬来源与真实CD4表达，再比较新的接收端关系')
ann('NME2','L12','NME2-KCa3.1-CD4活化已有机制；本队列为下调，应与下调FOXP1比较，不能因其不是上调基因而忽略。尚无充分排除依据。','保留比较候选','检验下调与CD4状态、钙信号及死亡程序的关系')
ann('CCR8','L13','CCL1-CCR8-Treg免疫调节已有研究，但当前ActCD4差异与该场景不完全相同，不能直接写成机制已成熟所以排除。','保留比较候选','复核亚群身份、受体检出及矽肺中的配体来源')
ann('LCN2','L14','肺损伤中成纤维细胞LCN2表达已有证据。CD4分群差异需确认来源，本轮不能把其他细胞证据写成CD4直接作用。','先核实细胞来源','结构细胞/髓系共检出、纯化T细胞表达与样本一致性')
ann('CCL5','L15','CCL5趋化功能已有长期实验基础。该队列ActCD4中下调，应作为状态或迁移读数保留，不能套用上调配体促病的叙事。','保留比较候选','同方向比较FOXP1，并核查接收端与来源特异性')
ann('GSN','L16','肺纤维化模型已有GSN功能证据，且涉及多种细胞与胞内/分泌形式；本队列CD4下调仍需细胞特异解释。','保留比较候选','区分异构体/来源，验证亚群内效应和死亡关系')
ann('ZEB2','L17','已有T细胞分化研究，但引用研究主要是CD8感染模型。不能据此排除矽肺CD4研究价值。','保留比较候选','比较CD4亚群内活性、靶基因及细胞状态关联')
ann('NTRK3','L18','TrkC/NT3在人体T细胞中已有研究，因此不能简单归为神经元污染。当前尚缺矽肺内配体覆盖与功能证据。','保留比较候选','核实CD4中的受体蛋白及NT3来源，再检验干预反应')
ann('GM2663','D01','功能注释多为预测，不能归类为机制成熟。先核对基因身份、人鼠对应及检测可靠性，再评估研究可行性。','注释与功能证据待补','基因标识、唯一比对、人鼠对应和独立表达验证')
ann('D830030K20RIK','D02','功能机制证据不足，不宜仅凭强差异直接设为主靶点，也不能称为研究过于成熟。','注释与功能证据待补','身份/转录本、人鼠对应及可验证功能的确认')

frames={};selected=[];bench=[];full=[]
for scope in ['Core','ActCD4','Naive','Treg']:
 d=pd.read_csv(S/f'{scope}_DE_all_genes.tsv.gz',sep='\t');assert not d.canonical.duplicated().any();d['scope']=scope;d['abs_logFC']=d.logFC.abs();d=d.sort_values(['abs_logFC','canonical'],ascending=[False,True]);sig=d[d.PValue<.05].copy();sig['nominal_abs_rank']=range(1,len(sig)+1);frames[scope]=d
 for g in CAND:
  q=d[d.canonical==g]
  bench.append({'状态':scope,'候选基因':g,'log2FC':None if q.empty else float(q.iloc[0].logFC),'原始P':None if q.empty else float(q.iloc[0].PValue),'原BH':None if q.empty else float(q.iloc[0].FDR),'检验状态':'未受检' if q.empty else '名义显著' if q.iloc[0].PValue<.05 else '未达名义显著'})
 c=d[d.canonical.isin(CAND)&d.PValue.lt(.05)]
 for _,r in sig.iterrows():
  if r.canonical in CAND:continue
  match=c[(np.sign(c.logFC)==np.sign(r.logFC))&(c.abs_logFC<abs(r.logFC))]
  anydir=c[c.abs_logFC<abs(r.logFC)]
  if not len(anydir):continue
  full.append({'状态':scope,'基因':r.canonical,'名义基因绝对效应排名':int(r.nominal_abs_rank),'log2FC':r.logFC,'logCPM':r.logCPM,'原始P':r.PValue,'原BH':r.FDR,'同方向且效应更低的候选':';'.join(match.canonical),'不分方向效应更低的候选':';'.join(anydir.canonical),'详细文献评估':'本表S17逐项评估' if r.canonical in ANN else '未逐项评估，不能视作已排除'})
 if scope in ['Core','ActCD4']:
  selected.extend(sig.head(20).loc[lambda x:~x.canonical.isin(CAND)].to_dict('records'))
assert set(r['canonical'] for r in selected)==set(ANN)
detail=[]
for g in sorted(ANN,key=lambda g:-max(abs(r['logFC']) for r in selected if r['canonical']==g)):
 a=ANN[g];ref=REF[a['ref']];rr=max([r for r in selected if r['canonical']==g],key=lambda x:abs(x['logFC']));f=next(x for x in full if x['基因']==g and x['状态']==rr['scope'])
 detail.append({'基因':g,'代表状态':rr['scope'],'log2FC':rr['logFC'],'原始P':rr['PValue'],'原BH':rr['FDR'],'同方向效应更低的候选':f['同方向且效应更低的候选'],'历史记录':a['history'],'文献具体发现':ref['finding'],'研究场景':ref['context'],'当前未优先/保留依据':a['reason'],'当前处置建议':a['status'],'再次考虑的条件':a['next_check'],'文献':f"{ref['year']} | {ref['title']}",'来源链接':ref['url']})
main=pd.DataFrame(detail);stats=[]
for scope in ['Core','ActCD4']:
 for g in ANN:
  q=frames[scope].loc[lambda d:d.canonical==g]
  if q.empty:
   stats.append({'状态':scope,'基因':g,'log2FC':None,'logCPM':None,'原始P':None,'原BH':None,'该状态进入前20展示':False,'检验状态':'未受检'})
   continue
  r=q.iloc[0];stats.append({'状态':scope,'基因':g,'log2FC':r.logFC,'logCPM':r.logCPM,'原始P':r.PValue,'原BH':r.FDR,'该状态进入前20展示':any(x['canonical']==g and x['scope']==scope for x in selected),'检验状态':'名义显著' if r.PValue<.05 else '未达名义显著'})

messages=json.loads((T/'shared_messages.json').read_text(encoding='utf-8'));ids=['f7ed6168-4c7e-5344-bea2-d89dfc500b06','79b3da88-d668-57b3-93b6-3abbadb88e55','b67aff77-3c97-5902-8215-974353b4dae9','47532a35-fce4-5231-810f-46d5a499650f','6dc3890e-58f3-5f67-9b0d-889b5ff69260','48168bce-3f64-5d11-9d10-b0db31216249']
history=[]
for i,mid in enumerate(ids,1):
 m=next(x for x in messages if x['id']==mid);history.append({'record':f'H{i:02d}','message_id':mid,'date_utc':datetime.datetime.fromtimestamp(m['create_time'],datetime.timezone.utc).isoformat(),'text':m['text'],'source':SHARE})
(O/'provenance/screening_history_excerpts.json').write_text(json.dumps(history,ensure_ascii=False,indent=2),encoding='utf-8')
histrows=[
 ['初筛导出','每张差异表按P<0.05及|log2FC|≥0.5，分别保留正/负方向前20个','H01','历史导出规则；不等于所有基因逐一完成文献排除'],
 ['IL1R2、PENK','优先表达核查，背景敏感性后仍优先','H02;H04;H05','本表不改变现有主候选'],
 ['TTC39C','探索候选，须确认真实T细胞表达','H03;H04','不能称为全新肺纤维化基因'],
 ['FOXP1','已有Treg机制基础，作为辅助/状态变化指标','H03;H04;H06','保留下降方向，避免仅选上调基因'],
 ['TNFRSF9、LAG3','已有矽肺/肺纤维化直接研究，降低创新优先级，作为机制参照','H03;H06','属于降级而非从统计候选中删除'],
 ['IGKC、IGHG2B、IGHG2C、JCHAIN','前列差异结果需要核查细胞来源','H02','历史未判定为已确认污染或逐基因排除'],
 ['背景转录本','检出率高不等于UMI负担高；NoBG不等于环境RNA去污染','H02;H05','本次不能用高检出率自动否定这些基因']]
hist=pd.DataFrame(histrows,columns=['对象','历史明确记录','记录编号','解释边界']);hist['原始对话']=SHARE
guide=[['用途','S17补充说明高效应基因未优先或仍应保留的依据，分开记录历史决定与本次文献评估。'],['基准候选','; '.join(CAND)],['比较口径','同一scope、冻结全基因差异表、原始P<0.05；效应为log2FC。先按绝对值比较，另列同方向的更低效应候选。'],['详细覆盖','Core与ActCD4名义差异按|log2FC|排名前20，去除11个既有重点基因后，35个基因-状态条目、27个独立基因。此为本次展示范围，不声称恢复历史逐基因淘汰清单。'],['完整范围','完整效应比较表包含4个CD4范围内高于至少一个名义候选的非重点基因；未逐项评估者明确标注，不能视作已排除。'],['历史来源','已读取用户提供的“数据集质量比较分析”分享页，提取筛选相关记录H01-H06。历史未记载的原因按本次补充评估表述。'],['成熟度解释','已有研究只影响重复创新空间；其他疾病/细胞中的机制不等于矽肺CD4问题已解决，不单独构成排除理由。'],['统计解释','不以FDR为准入硬门槛，不新增推断检验。每个状态为3对照+3硅尘小鼠；Core与亚群不作为独立队列。'],['方向说明','NME2、CCL5、GSN在ActCD4下调，按同方向与FOXP1比较，不把较大绝对值写成更强上调。'],['资料日期','文献及数据库核查截至2026-10-16。'],['评分限制','未编造综合创新性分数；未取得逐基因历史排除决定时不作追溯性断言。']]
frames_out={'逐基因说明':main,'历史筛选依据':hist,'详细统计':pd.DataFrame(stats),'候选基准':pd.DataFrame(bench),'完整效应比较':pd.DataFrame(full),'使用说明':pd.DataFrame(guide,columns=['项目','说明'])}
for name,d in frames_out.items():d.to_csv(O/'supplement'/f'S17_{name}.tsv',sep='\t',index=False,encoding='utf-8')
payload={'sheets':[{ 'name':k,'columns':list(v.columns),'rows':v.astype(object).where(pd.notna(v),None).values.tolist()} for k,v in frames_out.items()],'references':list(REF.values())}
(T/'workbook_data.json').write_text(json.dumps(payload,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
(O/'provenance/literature_sources.json').write_text(json.dumps(list(REF.values()),ensure_ascii=False,indent=2),encoding='utf-8')
counts={'detailed_unique_genes':len(main),'top20_gene_scope_entries':len(selected),'complete_comparison_rows':len(full),'candidate_benchmark_rows':len(bench),'history_records':len(history),'references':len(REF),'historical_share_sha256':hashlib.sha256((T/'shared_conversation.html').read_bytes()).hexdigest()}
(O/'provenance/S17_counts.json').write_text(json.dumps(counts,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(counts,ensure_ascii=False));print(main[['基因','当前处置建议']].to_string(index=False))
