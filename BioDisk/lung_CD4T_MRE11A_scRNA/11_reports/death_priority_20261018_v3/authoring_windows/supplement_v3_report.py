# Executed inside report.py. Frozen v2 results remain unchanged; v3 tables are read separately.
XT=P/'tables'
def xr(name):return pd.read_csv(XT/name,sep='\t')
def exfig(fid,title,insight,caption):
 global y
 f=P/'supplement_figures'/(fid+'.png')
 if not f.exists():return
 new('图'+str(len(figure_entries)+1)+'｜'+title);para(insight)
 figure_entries.append(dict(figure=len(figure_entries)+1,page=page,title=title,id=fid,source='v3补充分析统计表；'+fid))
 im=PILImage.open(f);iw,ih=im.size
 cp=Paragraph(caption,styles['small']);_,ch=cp.wrap(W-96,H)
 scale=min((W-96)/iw,(y-60-ch-12)/ih)
 c.drawImage(str(f),48+(W-96-iw*scale)/2,y-ih*scale,width=iw*scale,height=ih*scale);y-=ih*scale+10;para(caption,'small')

new('追加分析引入真实敲除证据，并检验跨队列复现')
done=(XT/'network_LOO_association_summary.tsv').exists()
table(['追加内容','实际执行与证据层级'],[
 ['独立同亚群复验','独立人肺Treg及CCR7阳性CD4亚群；供者级评分、同来源疾病差异、亚群细胞数敏感性。'],
 ['小鼠留一网络','All/BG1各10次逐小鼠排除、固定seed17重建；'+('已完成边稳定性与活性关联。' if done else '计算尚未结束，本页为内部排版草稿。')],
 ['跨队列方向','原始两个疾病队列的三个时间层分别计算；第二队列1＋1仅报告效应方向。'],
 ['公开真实干预','IL2RA小鼠Treg与HIF1A人T细胞计数重拟合；小鼠HIF1A完整作者差异统计量重新富集。'],
 ['计算扰动适用性','IL2RA/ARL6IP1不是SCENIC中的TF干预节点；HIF1A去重叠后仅All有6个目标、BG1无可评分自身regulon。本轮采用真实干预证据，未强行构造虚拟敲除P值。']
],[112,387])
para('独立人肺数据来自煤工尘肺背景，提供跨物种、相关肺病的外部检验；不是新的小鼠矽尘重复。原5＋5和主筛选不随补充结果改变。')
new('公开干预队列的设计和输入均已核查')
table(['队列','干预与样本','本轮分析'],[
 ['GSE121883','IL2RA/CD25条件敲除；脾Treg为诱导删除4＋4只小鼠；胸腺Treg靶向删除5＋5只','各组织独立edgeR拟合；程序效应、GSEA和关键基因。'],
 ['GSE253090','由Naive CD4激活48 h；常氧和低氧各WT 3库＋KO 3库','作者完整DESeq2结果重新做GSEA；不是原始计数重拟合。'],
 ['GSE324038','人原代混合T细胞，低氧HIF1A CRISPR；4名供者各KO和非靶向对照','8个文库、4个独立供者；供者配对edgeR模型。']
],[90,200,209])
para('GSE253090的Count.txt实际属于CUT&amp;RUN，未作为RNA计数使用；作者表中的重复Ensembl注释先核对统计量一致，再以唯一映射进入新富集。')
para('GSE324038计数矩阵覆盖第一批4名供者；第二批仅有元信息而不在该矩阵中的样本未计入。无CRISPR处理对照也未混作主比较的非靶向对照。')
para('来源：<link href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121883">GSE121883</link>；<link href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253090">GSE253090</link>；<link href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE324038">GSE324038</link>。','small')
new('IL2RA缺失使外周Treg的凋亡相关转录变化增强')
k=xr('external_key_genes.tsv');z=k[k.contrast=='GSE121883_spleen'].set_index('gene')
table(['基因','KO-WT log2FC','原始P','全基因FDR'],[[g,fmt(z.loc[g,'logFC']),fmt(z.loc[g,'PValue']),fmt(z.loc[g,'FDR'])] for g in ['IL2RA','BCL2','CASP3','BAX','PYCARD','HIF1A'] if g in z.index],[100,130,130,139])
para('重分析显示IL2RA下降，同时抗凋亡BCL2下降、CASP3上升。结合该研究的真实干预设计，结果支持IL2RA参与维持外周Treg存活，不能将原矽尘队列中的IL2RA上调直接写成促凋亡。')
para('原文另有Annexin-V、caspase-3/7和线粒体相关实验，支持CD25缺失后的存活障碍；这些是外部研究的实验证据，不是本轮新增湿实验，也不是对矽尘肺内IL2RA-NFKB2通路的救援验证。')
para('原文：<link href="https://www.nature.com/articles/s41467-019-08960-1">Toomer等，Nature Communications 2019，10:1037</link>。','small')
exfig('E01','IL2RA敲除的外周与胸腺Treg响应不同','程序效应必须按组织解释；外周Treg凋亡执行评分升高，胸腺没有相同响应。','独立小鼠为统计单位，外周n=8、胸腺n=10。每组织14类程序分别构成FDR家族。评分为过滤后成员logCPM均值；箭头不表示蛋白执行。')
new('外周Treg死亡评分升高，但不据此把死亡方式定为焦亡')
ef=xr('external_module_effects.tsv')
table(['外周Treg程序','KO-WT分数差','P','FDR'],[[r.program,fmt(r.effect),fmt(r.P),fmt(r.FDR)] for r in ef[(ef.contrast=='GSE121883_spleen')&ef.program.isin(['ApoExec','PyroExec','MitoROS','Glycolysis','OXPHOS'])].itertuples()],[135,125,115,124])
para('ApoExec与PyroExec均是RNA成员集合。PyroExec升高应完整保留，但没有新增gasdermin切割、膜孔形成或通路阻断数据，不能用这个评分重新认定原文的死亡方式。')
para('外周Treg糖酵解和OXPHOS平均评分未达到P＜0.05；胸腺的相应评分下降。两种组织与实验背景下的响应不同，不能将所有组织合并来“验证”原肺Treg的代谢关联。')
new('外周Treg敲除增强凋亡程序，并伴随ROS及焦亡转录信号')
pg=xr('external_module_GSEA.tsv');pg=pg[pg.contrast=='GSE121883_spleen'].set_index('pathway')
kg=xr('GSE121883_spleen_KEGG_GSEA.tsv').set_index('pathway')
rr=[]
for mod in ['ApoExec','PyroExec','MitoROS','Glycolysis','OXPHOS']:
 row=pg.loc[mod];rr.append(['预设程序 / '+mod,fmt(row.NES),fmt(row.pval),fmt(row.padj)])
for code,label in [('mmu04210','凋亡'),('mmu00100','类固醇生物合成'),('mmu00010','糖酵解/糖异生'),('mmu00190','氧化磷酸化'),('mmu04979','胆固醇代谢')]:
 row=kg.loc[code];rr.append(['KEGG / '+label,fmt(row.NES),fmt(row.pval),fmt(row.padj)])
table(['检验家族/条目','NES','P','FDR'],rr,[220,85,95,99])
para('凋亡执行程序与KEGG凋亡均上移并通过各自家族FDR，类固醇生物合成下降。胆固醇代谢总通路未获得同等支持，不能把局部生物合成结果扩大为全部脂代谢下降。')
para('提高数值精度后，PyroExec与MitoROS的全基因排序GSEA均达到原始P＜0.05，但FDR未通过，属于名义上移；平均评分则通过其程序家族FDR。两种统计回答不同问题，证据强度分别呈现。凋亡执行程序同时获得平均评分和排序富集的FDR支持。')
new('外周与胸腺敲除实验的NFKB2靶程序方向相反')
rg=xr('external_joint_regulon_GSEA.tsv')
q=rg[rg.pathway.eq('NFKB2')&rg.contrast.str.startswith('GSE121883')]
table(['组织/网络','覆盖靶基因','NES','P','完整TF家族FDR'],[[r.contrast.replace('GSE121883_','')+' / '+r.arm,int(r['size']),fmt(r.NES),fmt(r.pval),fmt(r.padj)] for _,r in q.iterrows()],[125,85,80,100,109])
para('外周Treg中，IL2RA缺失后NFKB2靶基因集合上移；胸腺Treg中则下降，两个网络范围方向一致。该结果说明真实干预能够改变相关转录程序，但不支持将原肺内的正相关简化为“IL2RA单向正向激活NFKB2”的通路结论。')
para('这些签名全部来自当前5＋5联合网络，且已移除候选、程序成员和TF自身。富集测量靶基因的排序变化，不直接测量NFKB2蛋白活化，也不是NFKB2阻断或救援实验。')
para('All/BG1、5个外部对比分别评价完整可评估TF家族，共1,815项；完整结果和靶覆盖均保留。HIF1A去重叠小签名没有在两套真实HIF1A干预中获得名义下降支持，未据此强行构造虚拟敲除结论。')
new('小鼠HIF1A敲除支持代谢调控，焦亡评分方向与原相关性不一致')
hf=xr('mouse_HIF1A_perturbation_GSEA.tsv');z=hf[hf.database=='Modules']
table(['氧条件/程序','NES','P','FDR'],[[r.contrast.replace('GSE253090_','')+' / '+r.pathway,fmt(r.NES),fmt(r.pval),fmt(r.padj)] for r in z[z.pathway.isin(['Glycolysis','OXPHOS','ApoExec','PyroExec'])].itertuples()],[205,85,100,109])
para('低氧HIF1A敲除后，糖酵解排序下降、OXPHOS排序上升；PyroExec则上升。后者不支持把原Naive中HIF1A与PyroExec的正相关直接解释为“HIF1A促焦亡”。凋亡RNA集合未获得同等富集支持。')
para('RNA采于激活48 h，原文死亡表型包含后续培养时间；不能要求RNA执行评分等同于后续死亡数量。原文报告HIF1A-糖酵解与低氧激活诱导死亡的实验关系，外推范围是该激活体系。')
para('来源：<link href="https://www.nature.com/articles/s41467-024-53593-8">HIF1α-regulated glycolysis promotes activation-induced cell death，Nature Communications 2024</link>。','small')
exfig('E02','真实HIF1A干预揭示氧条件依赖的程序变化','同一候选的干预效应随细胞状态和氧条件改变，原Naive关联不能替代干预效应。','小鼠图用作者全基因Wald统计量重新GSEA；人图用4名配对供者重新拟合的signed sqrt(F)。两者的NES与P来自不同检验家族，不作绝对强弱排名。')
new('人HIF1A敲除提供配对供者证据，但不复现全部小鼠结果')
z=ef[ef.contrast=='GSE324038_HumanPanT_Hypoxia']
table(['程序','KO-WT分数差','P','FDR'],[[r.program,fmt(r.effect),fmt(r.P),fmt(r.FDR)] for r in z[z.program.isin(['ApoExec','PyroExec','Glycolysis','OXPHOS','MitoROS'])].itertuples()],[135,125,115,124])
r=k[(k.contrast=='GSE324038_HumanPanT_Hypoxia')&k.gene.eq('HIF1A')].iloc[0]
para(f'HIF1A基因表达变化为log2FC={fmt(r.logFC)}，P={fmt(r.PValue)}，FDR={fmt(r.FDR)}。8个文库来自4名供者，推断样本量为4。')
para('OXPHOS平均评分为名义上升，未通过程序家族FDR；凋亡和焦亡评分没有名义变化。这是混合原代T细胞，不等于肺Naive亚群。')
eg=xr('external_module_GSEA.tsv');eg=eg[eg.contrast=='GSE324038_HumanPanT_Hypoxia'].set_index('pathway')
para('全基因排序GSEA则支持糖酵解下降与OXPHOS上升：'+ '；'.join(mod+' NES='+fmt(eg.loc[mod,'NES'])+'，P='+fmt(eg.loc[mod,'pval'])+'，FDR='+fmt(eg.loc[mod,'padj']) for mod in ['Glycolysis','OXPHOS'])+'。这与小鼠低氧干预方向一致；排序富集与成员平均评分检验的问题不同，两种统计均保留。')
para('小鼠Hif1a敲除的作者表中，低氧Hif1a整基因RNA反而升高；整基因RNA丰度不等同于完整功能蛋白。本轮保留该值，敲除身份依据原实验设计，而不是强行要求整基因RNA下降。')
new('独立人肺亚群未复现IL2RA的疾病上调')
h=xr('human_G227_candidate_DE.tsv');hq=xr('human_G227_DE_QC.tsv')
table(['亚群','对照/病例','用途'],[[r.scope,f'{r.n_control} / {r.n_case}','同一G227来源内比较'] for r in hq.itertuples()],[130,115,254])
pairs=[('IL2RA','Treg'),('HIF1A','CD4-CCR7'),('ARL6IP1','CD4-CCR7')]
rr=[]
for g,sc in pairs:
 z=h[(h.gene==g)&(h.scope==sc)]
 if len(z):
  r=z.iloc[0];rr.append([g+' / '+sc,fmt(r.logFC),fmt(r.PValue),fmt(r.FDR)])
table(['候选/亚群','病例-对照log2FC','P','FDR'],rr,[185,115,95,104])
para('人肺Treg的IL2RA疾病效应方向与小鼠联合结果相反，且未达到名义显著。HIF1A和ARL6IP1在CCR7阳性CD4中也未复现名义疾病差异。CCR7阳性CD4仅作为Naive相关亚群对照，不宣称与小鼠Naive完全等价。')
para('仅使用同一来源内同时存在的病例与对照进行主疾病比较，避免把G288纯病例和其他研究的对照直接混成独立复验。全部来源的关联另作来源调整敏感性。')
exfig('E03','独立人肺的候选程序关联呈现部分方向支持与不一致','IL2RA-糖酵解为正向但未达到P＜0.05；HIF1A-焦亡评分没有复现。','图为同来源G227，亚群至少10细胞/供者。各候选14类程序为独立FDR家族；另保留至少20细胞及来源/疾病调整敏感性。')
new('多来源人肺Treg支持IL2RA与代谢程序的正相关')
ha=xr('human_candidate_associations.tsv')
hh=ha[(ha.cohort=='AllHuman')&(ha.min_cells==10)&ha.gene.eq('IL2RA')&ha.program.isin(['Glycolysis','OXPHOS','PyroExec'])]
table(['模型/程序','n','ρ','P','FDR'],[[('来源＋疾病' if r.model=='SourceCondition' else '来源')+' / '+r.program,int(r.n),fmt(r.rho),fmt(r.P),fmt(r.FDR)] for r in hh.itertuples()],[200,35,80,85,99])
para('13例多来源Treg的糖酵解与OXPHOS正相关在来源＋疾病调整后的14程序检验家族内均通过FDR，支持候选的代谢关联。G227同来源8例的效应同向，但未达到名义P；两种分析共享8名供者，不计作两次独立复验。')
para('多来源模型评价样本间的程序关联，不提供额外的同来源疾病病例-对照验证。HIF1A-焦亡及ARL6IP1-ROS在13例来源调整后仍未达到名义P；对应完整结果保留。')
new('人肺Treg的IL2RA-焦亡评分关联与小鼠主分析不同')

ha=xr('human_candidate_associations.tsv');z=ha[(ha.cohort=='G227')&(ha.min_cells==10)&(ha.model=='Source')]
targets=[('IL2RA','Glycolysis'),('IL2RA','OXPHOS'),('IL2RA','PyroExec'),('HIF1A','PyroExec'),('ARL6IP1','MitoROS')]
rr=[]
for g,mod in targets:
 r=z[(z.gene==g)&(z.program==mod)].iloc[0];rr.append([g+' / '+mod,int(r.n),fmt(r.rho),fmt(r.P),fmt(r.FDR)])
table(['关联','独立供者','ρ','P','FDR'],rr,[175,75,75,85,89])
para('IL2RA-PyroExec在人肺Treg中达到名义关联，但原小鼠Treg主分析不支持这一关系，G227主分析未通过14程序家族FDR，而13例多来源调整分析通过FDR。至少20细胞的多来源10例敏感性中，来源＋疾病调整后的关联仍通过FDR。该跨物种差异完整保留，不将RNA关联写成死亡执行。')
new('同来源关键评分成员复验未达到名义显著')
para('本页使用G227同来源子集，分别为9例CCR7阳性CD4与8例Treg；与13例多来源程序关联分开解释。','small')
km=xr('human_targeted_member_associations.tsv')
table(['候选/成员','亚群','n','ρ','P / FDR'],[[r.gene+' / '+r.member,r.scope,int(r.n),fmt(r.rho),fmt(r.P)+' / '+fmt(r.FDR)] for r in km[km.model=='Unadjusted'].itertuples()],[135,100,35,75,154])
para('HIF1A-PYCARD未显示关联；ARL6IP1-GCLM及IL2RA-CYC1保留正向，但P未达到0.05。不能把方向一致描述为已复现的名义关联。疾病调整后的完整结果保留在表中。')
new('独立矽肺BAL的亚群覆盖不足，只保留描述性结果')
bal=xr('BAL174725_sample_scores.tsv')
table(['供者','组别','NaiveCD4细胞','Treg细胞'],[[sam,z.group.iloc[0],int(z.loc[z.subtype=='NaiveCD4','cells'].sum()),int(z.loc[z.subtype=='Treg','cells'].sum())] for sam,z in bal.groupby('sample')],[140,115,122,122])
para('GSE174725包含2名暴露对照与3名矽肺患者。要求亚群至少10细胞后，NaiveCD4和Treg均只剩1名暴露对照与2名患者；不作新的组间或相关性推断检验。全部5名供者及10/20细胞门槛下的效应表均保留。')
para('因此，已完成公开独立队列的适用性检查与可用数据分析，但没有获得新的充分重复、同亚群小鼠矽尘验证。')
para('来源：<link href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174725">GSE174725</link>。','small')
new('三个时间层的候选均值差同向，疾病效应仍以模型为准')
dc=xr('cross_cohort_direction.tsv')
rr=[]
for g,sc in [('HIF1A','Naive'),('IL2RA','Treg'),('ARL6IP1','Naive')]:
 for r in dc[(dc.variant=='Raw')&(dc.scope==sc)&dc.feature.eq('G_'+g)].itertuples():
  rr.append([g+' / '+sc,r.dataset+' D'+str(r.day),fmt(r.effect),f'{r.n_control}＋{r.n_case}'])
table(['候选/亚群','队列/时间','病例-对照分数差','样本'],rr,[135,155,135,74])
para('三个主候选在各层的逐细胞log1p均值差均为正向；这一指标与伪bulk总计数及TMM归一化的edgeR log2FC不等价，尤其ARL6IP1不能据此改判疾病效应。第二队列每时间点仅1＋1，只报告方向；原5＋5中的两个队列属于发现数据内部一致性检验，不是新的独立复验。')
new('逐小鼠重建检验了调控边对独立样本的敏感性')
if done:
 la=xr('network_LOO_nonself_fit_audit.tsv')
 table(['范围','完成网络','motif边数范围','与完整seed17的Jaccard范围'],[[arm,len(z),str(int(z.edges.min()))+' - '+str(int(z.edges.max())),fmt(z.jaccard.min())+' - '+fmt(z.jaccard.max())] for arm,z in la.groupby('arm')],[65,85,145,204])
 para('剔除TF自身成员后，全图边集合的Jaccard约为0.19至0.26，显示不少边随留一样本发生改变；局部候选边的保留不能泛化为整个网络高度稳定。')
 para('每个范围逐一去除10只小鼠，保持原seed17的细胞抽取和基因背景，去除零方差基因后重建GRNBoost2并重新motif剪枝。这样把供者删除效应与重新抽取细胞的随机变化分开。')
 ls=xr('network_LOO_association_summary.tsv');rr=[]
 for arm in ['All','BG1']:
  q=ls[(ls.arm==arm)&ls.gene.eq('IL2RA')&ls.tf.eq('NFKB2')&ls.model.eq('JointStratum')]
  if len(q):
   r=q.iloc[0];rr.append([arm,int(r.folds_evaluable),f'{int(r.same_direction)}/10',f'{int(r.nominal)}/10',fmt(r.rho_min)+' - '+fmt(r.rho_max)])
 table(['IL2RA-NFKB2','可评估','同向','P＜0.05','ρ范围'],rr,[100,70,80,85,164])
 para('每折活性签名重新去除全部候选、程序成员与TF自身，并只对剩余9只小鼠检验；另外对完整10样本的seed17网络采用相同评分和统计口径建立匹配参考。各折的P与FDR保留完整可评估候选连接家族；重复留一不增加独立样本量。')
else:para('网络仍在运行；最终交付前必须检查20个网络的GRN与motif步骤退出码、完成标记及输出完整性。')
new('匹配单种子参考后，IL2RA-NFKB2证据强度仍呈敏感性')
mr=xr('network_LOO_matched_summary.tsv')
mr=mr[(mr.gene=='IL2RA')&(mr.tf=='NFKB2')&(mr.model=='JointStratum')]
table(['范围','共识ρ / P / FDR','seed17 ρ / P / FDR','9样本同向 / 名义折数'],[[r.arm,fmt(r.baseline_rho)+' / '+fmt(r.baseline_P)+' / '+fmt(r.baseline_FDR),fmt(r.full_seed17_rho)+' / '+fmt(r.full_seed17_P)+' / '+fmt(r.full_seed17_FDR),f'{int(r.same_direction_to_full_seed17)}/10 / {int(r.nominal)}/10'] for r in mr.itertuples()],[60,150,150,139])
para('All完整seed17本身P=0.0566，未达名义阈值；其10次留一中9次同向、0次名义显著，不能把与共识结果的差异全部归因于删去一只小鼠。BG1完整seed17为名义关联，FDR=0.1038；10次留一全部同向、7次达到名义P。')
para('三种子共识仍是原主分析；单种子参考与留一结果分别评价签名构建和供者影响。两范围的留一关联均未通过各折完整候选连接家族FDR；该结果不否定多数折保持的方向。')
edges=xr('network_LOO_candidate_edges_annotated.tsv')
ee=edges[edges.tf.eq('NFKB2')&edges.target.eq('IL2RA')]
para('NFKB2-IL2RA的motif边在All为9/10、BG1为7/10折出现，两者都存在于完整seed17参考。结构先验的保留与RNA-活性关联的显著性属于不同证据层级，不能合并为已确立的通路。')
para('疾病条件调整后，留一置换空间为12或36种；其中12种的最小可达P为1/12，无法达到P＜0.05。此类敏感性须结合ρ及其范围，不以未过P/FDR阈值单独否定方向。','small')
exfig('E04','候选调控边的供者留一稳定性','稳定边须同时审视原三种子支持和逐小鼠剔除后的保留次数。','展示原All三种子稳定候选边在10折中的出现次数，剔除TF自身成员；颜色区分完整seed17基准中是否存在。部分2/3边不在seed17中，不能将其低频直接归因于供者删除。多数折出现也不等于直接结合或因果。')
new('综合定位：IL2RA偏存活调控，HIF1A偏状态依赖代谢调控')
table(['候选','本轮新增证据后的定位'],[
 ['IL2RA / Treg','原小鼠肺疾病上调及代谢关联保留；多来源人肺Treg的代谢正相关获FDR支持，同来源疾病上调未复现。外部条件敲除支持外周Treg存活与抗凋亡作用，NFKB2通路因果关系未建立。'],
 ['HIF1A / Naive','真实敲除支持CD4激活背景下的代谢重编程；原Naive焦亡评分关联没有得到同方向干预复现。不能据此将该候选定为促焦亡因子。'],
 ['ARL6IP1 / Naive','维持ROS候选线索定位；本轮未找到适配CD4/Treg且可重分析的公开敲除转录组。其他细胞类型的实验不能代替本项目亚群验证。']
],[125,374])
para('“凋亡相关”包括维持存活或抑制凋亡，不预设所有候选都促进死亡。完成生信分析不以实验机制闭环为门槛；机制因果表述则必须对应实际干预、通路及死亡执行证据。')
new('方法：外部复验与真实干预分析')
para('外部计数采用Ensembl到物种对应org.*.eg.db的唯一符号映射；歧义映射不强行选择，保留映射审计。过滤后TMM、edgeR robust quasi-likelihood；GSE121883按组织分别拟合，GSE324038使用供者＋干预设计。RSEM expected counts保留原数值，不人为取整。')
para('程序评分为预设且去除7个候选的成员logCPM均值，以独立小鼠/供者拟合线性模型。成员覆盖随队列过滤而变化，完整成员数和列表随表保留，不跨数据类型比较绝对分数。每对比14程序效应、预设程序GSEA和KEGG GSEA分别校正。')
para('GSE253090复用作者全基因Wald统计量，不重新拟合基因P值；新GSEA采用fgseaMultilevel、eps=0。重复注释核对统计量一致；多个Ensembl映射到同一符号时按baseMean优先，避免按显著性挑选。')
para('人肺沿用冻结细胞注释，亚群每供者至少10细胞，至少20细胞作为敏感性。程序采用冻结成员的同名符号投射，并保留实际成员覆盖。相关性调整来源或来源＋疾病；可枚举时精确排列，否则19,999次排列并使用加一修正。')
para('真实干预中的TF靶程序GSEA采用当前All/BG1稳定去重叠签名，至少5个被测试靶基因才评价；每个外部对比、每个网络范围的完整可评估TF为独立校正家族，不与其他程序或KEGG合并。')
para('新增GSEA初次P在0.04至0.06之间的153项统一进行一次高精度计算：nPermSimple=100,000、sampleSize=1,001、eps=0、固定seed=20261018。使用高精度结果重算原完整家族FDR，保留初次表与冻结清单。')
new('补充材料：新增分析、输入核查与边界')
table(['内容','v3材料'],[
 ['真实干预','external_fit_QC.tsv；external_key_genes.tsv；external_module_effects.tsv；external_module_GSEA.tsv；mouse_HIF1A_perturbation_GSEA.tsv；external_joint_regulon_GSEA.tsv'],
 ['独立人肺','human_coverage.tsv；human_sample_scores.tsv；human_G227_candidate_DE.tsv；human_candidate_associations.tsv；human_targeted_member_associations.tsv；BAL174725_*.tsv'],
 ['跨队列','cross_cohort_direction.tsv；独立验证适用性与检索记录'],
 ['供者网络稳定性','network_LOO_nonself_fit_audit.tsv；network_LOO_candidate_edges.tsv；network_LOO_association_summary.tsv；network_LOO_matched_summary.tsv；network_seed17_reference_associations.tsv；computational_perturbation_suitability.tsv'],
 ['复核及存储','supplement_QA.json；storage_migration.json；system_tmp_migration.json；原始数据SHA256清单']
],[110,389])
para('未检索到适配数据不代表不存在相应生物学作用。ARL6IP1已检索GEO及相关原始研究；检索到的神经/成纤维等细胞研究仅提供背景，未替代CD4干预复验。')
para('公共数据保存在/mnt/BioData/bio/wangmingyuan/lung_CD4T_MRE11A_scRNA下对应目录，旧分析路径以链接兼容。','small')
