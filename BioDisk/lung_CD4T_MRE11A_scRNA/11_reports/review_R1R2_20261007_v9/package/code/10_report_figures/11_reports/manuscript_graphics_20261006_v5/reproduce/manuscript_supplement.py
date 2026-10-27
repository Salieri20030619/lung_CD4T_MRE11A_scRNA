doc.add_page_break();h('4　统计与数值补充',1)
p('以下补充表由原始聚合结果直接生成。主图所示正负效应均保留；“不可估计”只用于原结果缺失，不等同不显著。完整候选、各门槛组合、配体–靶成员和背景成员另存服务器可复现档案；本报告不复制个体供者行或细胞表达矩阵。')
scope_map={'Primary':'全体CD4','Core':'核心CD4','CD4-CCR7':'CCR7','CD4-GZMA':'GZMA','CD4-IL7R':'IL7R','Treg':'Treg','CD4state':'CD4state','ActCD4':'ActCD4','Naive':'Naive'}
sc=lambda v:scope_map.get(str(v).replace('Subtype:',''),str(v))
model_map={'Base':'基础','Source':'来源校正','SourceDep':'来源与深度校正','WithinG227':'G227来源内','SourceAdj':'来源校正','QCAdj':'来源与质量校正'}
arm={'All':'完整细胞','BG1':'背景筛选'}
variants={'Original':'原始调控子','DropTarget':'去候选基因','DropOverlap':'去候选及程序重叠'}
def triple(r,a,b,c):return f'{n(r[a])} / {n(r[b])} / {n(r[c])}'
h('4.1　统计模型及校正范围')
table('读数、模型与多重检验家族',['分析对象','供者模型或统计量','BH家族与特殊限制'],[
['组成/表达匹配程序','Y~G；Y~G+S；Y~G+S+log1p(U)','12程序按范围×背景×模型校正；窄面板跨范围家族另列'],
['新增RNA/AUC疾病模型','Y~G；Y~G+S；Y~G+S+log10(U)','读数×版本×种子×背景×模型内，跨范围×两个TF'],
['历史AUC组间','小鼠Y~G或+log1p(U)；人类Y~G+S或+log1p(U)','旧规则下每组不足10未算组间FDR；保留原P'],
['调控子–程序相关','对AUC/评分取秩；残差化G、S和秩(U)；Pearson相关残差','物种×背景×种子×版本×范围×模型内全部受检TF×程序/候选表达特征；n≥10才算原P/BH'],
['人类UCell','Y~G+S；再加z[log10(U中位数+1)]、z[线粒体比例中位数]','每范围×模型跨9模块；亚群模型每亚群跨9模块'],
['人类GSEA/GSVA','edgeR~S+G排序；GSVA后limma~S+G','每范围×数据库；敏感性各自保留原BH'],
['NicheNet疾病DE','伪批量edgeR~G或~G+S；完整文库×TMM偏置','受检基因内BH；靶集和背景门槛独立'],
['配体同供者关联','两端秩读数残差化G、S、两端秩(U)','范围×背景×DE模型×靶集阈值×方向内跨配体–发送端'],
['分层Bootstrap','整行供者，按S×G分层，5000次；百分位95%CI','无新P/BH；失败不补抽，配对差保持共享索引'],
['小鼠UCell','20种3+3标签分配的双侧精确置换','8模块×3主要范围；Treg另作探索家族']],[3.3,7.1,6.8],note='G为疾病，S为研究来源，U为相应范围平均原始UMI；z为标准化。来源与深度的具体变换依读数而定，不互换。背景筛选为CD4固定背景面板UMI比例≤1%，巨噬细胞不套用该筛选。',supp=True)
h('4.2　固定程序成员与评分覆盖')
d=df('module_definitions');rows=[]
cn={'OXPHOS':'核编码氧化磷酸化','MitoISR':'线粒体整合应激','MitoROS':'线粒体氧化还原反应','Mitophagy':'线粒体自噬','Inflam':'炎症小体','PyroExec':'焦亡执行相关','ApoExec':'凋亡执行相关','TAct':'T细胞活化','TReg':'调节性相关','Inhib':'抑制相关','Naive':'初始状态','Cycle':'细胞周期'}
for mod,g in d.groupby('module',sort=False):rows.append([mod,cn[mod],len(g),', '.join(g.gene)])
table('表达匹配评分的完整签名成员',['缩写','程序','成员数','基因'],rows,[2.0,3.2,1.3,10.7],note='共12个固定手工转录面板，不称为官方GO通路。控制基因按方法1.6固定；原始可测、表达与控制基因覆盖全部保存在服务器。',supp=True)
rows=[]
for _,r in df('human_UCell_module_gene_coverage').iterrows():rows.append([r.module,int(r.expected_n),int(r.present_n),r.present_genes])
table('人类UCell面板与实际成员覆盖',['面板','定义数','可测数','成员'],rows,[2,1.2,1.2,12.8],note='可测表示存在于对象特征空间，不能据此声称所有供者或所有细胞均检出。MitoResp与PyroCore为相应子面板并集，非独立验证；UCell无独立凋亡执行面板。',supp=True)
h('4.3　通路富集与评分的对应数值')
terms=['KEGG_OXIDATIVE_PHOSPHORYLATION','GOBP_OXIDATIVE_PHOSPHORYLATION','GOBP_REGULATION_OF_EXECUTION_PHASE_OF_APOPTOSIS','GOBP_POSITIVE_REGULATION_OF_EXECUTION_PHASE_OF_APOPTOSIS','GOBP_EXECUTION_PHASE_OF_APOPTOSIS','KEGG_APOPTOSIS','KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY']
tn=['KEGG氧化磷酸化','GO氧化磷酸化','GO凋亡执行调控','GO凋亡执行正调控','GO凋亡执行','KEGG凋亡','KEGG TCR']
g=df('human_GSEA_focused');v=df('human_GSVA_focused');rows=[]
for t,name in zip(terms,tn):
 for scope in ['Primary','Core']:
  r=g[(g.ID==t)&(g.analysis==scope)].iloc[0];s=v[(v.pathway_id==t)&(v.analysis==scope)].iloc[0];rows.append([name,sc(scope),triple(r,'NES','pvalue','p.adjust'),triple(s,'logFC','P.Value','adj.P.Val')])
table('人类GSEA与GSVA对应比较',['通路','范围','GSEA：NES / P / BH','GSVA：评分差 / P / BH'],rows,[4.1,2,5.6,5.5],note='全部为来源校正、n=17（8对照/9病例）；GO与KEGG分别校正。GSEA的负向凋亡条目不等于死亡率下降。',supp=True)
# The next table number is referenced explicitly in the text.
h('4.4　IKZF2死亡相关程序的完整种子与去重叠敏感性')
d=df('all_TF_program_associations');d=d[(d.species=='Human')&(d.scope=='Core')&(d.model=='SourceDep')&(d.tf=='IKZF2')&d.feature.isin(['score_PyroExec','score_ApoExec'])];rows=[]
for ar in ['All','BG1']:
 for va in ['Original','DropTarget','DropOverlap']:
  for feature in ['score_PyroExec','score_ApoExec']:
   rrs=[d[(d['filter']==ar)&(d.variant==va)&(d.feature==feature)&(d.seed==seed)].iloc[0] for seed in [17,29,43]]
   rows.append([arm[ar]+'\n'+variants[va],feature.replace('score_',''),str(int(rrs[0]['n']))+'/'+str(int(rrs[0]['df']))]+[triple(r,'rho','p','FDR') for r in rrs])
table('IKZF2与凋亡/焦亡评分：所有预定调控子版本',['背景/调控子','程序','n/df','种子17：ρ / P / BH','种子29：ρ / P / BH','种子43：ρ / P / BH'],rows,[3.2,1.7,1.1,3.75,3.75,3.7],note='核心CD4，均控制疾病、来源和接收者秩深度。三个种子不作合并检验。原BH在每物种×背景×种子×版本×范围×模型内，跨受检TF与12程序及7候选表达特征，取有限P值校正。完整细胞三个种子的有限检验数为918/990/990，受检TF为51/55/55；背景筛选后为972/990/972及54/55/54，三个调控子版本各自使用相同对应家族大小。受检TF为候选上游非自身TF及八个既有TF，不是仅本表IKZF2。去候选版本移除既有七基因及TF本身，去重叠再移除12程序成员。背景筛选同时改变供者数，不能把差异只归因于评分。',supp=True)
h('4.5　完整细胞中两个候选的全部程序关联')
d=df('all_TF_program_associations').query("species=='Human' and filter=='All' and variant=='DropOverlap' and scope=='Core' and model=='SourceDep'");rows=[]
for tf in ['IKZF2','STAT1']:
 for mod in cn:
  rr=[d[(d.tf==tf)&(d.feature=='score_'+mod)&(d.seed==seed)].iloc[0] for seed in [17,29,43]];rows.append([tf,mod]+[triple(r,'rho','p','FDR') for r in rr])
table('核心CD4的12程序关联，包括未支持结果',['TF','程序','种子17：ρ / P / BH','种子29：ρ / P / BH','种子43：ρ / P / BH'],rows,[1.7,2.2,4.4,4.4,4.5],note='所有行n=17、残差df=11；完整细胞、去候选及程序重叠调控子、来源与深度校正。小鼠相同相关分析n=6，原研究未计算相关P/BH，方向结果保存于服务器，不用空FDR列表示阴性。',supp=True)
h('4.6　线粒体及细胞死亡程序的疾病比较')
d=df('human_all_12_programs');rows=[]
for scope in ['Primary','Core','Subtype:CD4-CCR7','Subtype:CD4-GZMA','Subtype:CD4-IL7R','Subtype:Treg']:
 for mod in ['OXPHOS','MitoISR','MitoROS','Mitophagy','Inflam','PyroExec','ApoExec']:
  rr=[d[(d.scope==scope)&(d['filter']=='All')&(d.model==mo)&(d.metric=='score_'+mod)].iloc[0] for mo in ['Source','SourceDep']]
  rows.append([sc(scope),mod,f'{int(rr[0].n_control)}/{int(rr[0].n_case)}']+[triple(r,'coefficient','p_nominal','FDR') for r in rr])
table('表达匹配评分的来源及深度模型',['范围','程序','对照/病例','来源校正：β / P / BH','来源与深度：β / P / BH'],rows,[2,2.1,1.7,5.7,5.7],note='完整细胞、每供者≥20细胞；BH跨同范围/背景/模型的全部12程序，不仅本表7程序。全套504项模型（含免疫程序、背景筛选及来源内比较）保存服务器；不可估计组合保留其原因。',supp=True)
d=df('human_UCell_module_group_tests');rows=[]
for scope in ['Primary','Core']:
 for mod in df('human_UCell_module_gene_coverage').module:
  rr=[d[(d.analysis==scope)&(d.module==mod)&(d.model==mo)].iloc[0] for mo in ['SourceAdj','QCAdj']];rows.append([sc(scope),mod]+[triple(r,'estimate','p_value','fdr_bh') for r in rr])
table('人类UCell的全部整体疾病比较',['范围','面板','来源校正：β / P / BH','来源与质量：β / P / BH'],rows,[2.2,2.5,6.2,6.3],note='所有行8对照/9病例，共36个模型结果；未有P<0.05。亚群≥10细胞的45项结果中，IL7R线粒体自噬及替代Gasdermin为名义提示，BH均0.073532，其余未达名义阈值或不可估计。',supp=True)
h('4.7　接收者与发送端的同供者配对')
d=de('pairing');rows=[]
for rec in ['Primary','Core','CD4-CCR7','CD4-GZMA','CD4-IL7R','Treg','CD4state']:
 for sender in ['MacAM','MonoMac','MacWide']:
  z=d[(d.receiver==rec)&(d.sender==sender)];cc=z[z.group=='Control'];ww=z[z.group=='CWP']
  if len(cc) and len(ww):rows.append([sc(rec),sender,int(cc.iloc[0].paired_donors),int(ww.iloc[0].paired_donors),f'{int(cc.iloc[0].receiver_cells)}/{int(ww.iloc[0].receiver_cells)}',f'{int(cc.iloc[0].sender_cells)}/{int(ww.iloc[0].sender_cells)}'])
table('两端细胞覆盖与合格配对供者',['接收范围','发送端','对照配对','病例配对','接收细胞 对照/病例','发送细胞 对照/病例'],rows,[2.2,2.2,1.8,1.8,4.6,4.6],note='配对供者要求两端均≥20细胞；细胞数列为该组范围内的聚合覆盖，未限于达到配对门槛者，不能替代配对数。MacWide为并集，各发送端重复出现的供者不能累加。',supp=True)
h('4.8　配体活性与同供者关联')
d=de('stage2_targetset_status');rows=[]
for rec in ['Primary','Core','CD4-CCR7','CD4-GZMA','CD4-IL7R','Treg','CD4state']:
 z=d[(d.receiver==rec)&(d['filter']=='All')&(d.model=='Source')];vals=[]
 for test,dr in [('FDR05','Up'),('FDR05','Down'),('Nominal05','Up'),('Nominal05','Down')]:
  rr=z[(z.test==test)&(z.direction==dr)];r=rr.iloc[0] if len(rr) else None;vals.append('不足/不可估计' if r is None or pd.isna(r.geneset_n) else str(int(r.geneset_n)))
 row0=z.iloc[0];rows.append([sc(rec),'不可估计' if pd.isna(row0.background_n) else int(row0.background_n)]+vals)
table('来源校正后人类各接收者的靶集与背景',['范围','背景基因','FDR上调','FDR下调','名义上调','名义下调'],rows,[2.4,2.4,3.1,3.1,3.1,3.1],note='所有靶集另需|log₂FC|≥0.5；FDR及名义靶集不可合并。至少10靶基因、500背景基因等门槛均通过后才有活性结果。',supp=True)
d=de('human_NN_BH');rows=[]
for _,r in d.iterrows():rows.append([sc(r.receiver),r.sender+'/'+r.ligand,arm[r['filter']],model_map[r.model]+'\n'+('FDR' if r.test=='FDR05' else '名义')+('上调' if r.direction=='Up' else '下调'),str(int(r['n']))+'/'+str(int(r['df'])),n(r.rho),n(r.P),n(r.FDR)])
table('全部九条通过原BH的同供者关联',['接收者','发送端/配体','背景','疾病靶集','n/df','ρ','P','BH'],rows,[1.8,3.1,1.8,3,1.3,1.9,2.1,2.2],note='这些是同供者关联，不能等同NicheNet配体活性通过。CCR7相关及GZMA–VTN的删一供者可估计性极低；CXCL2所配FDR靶集不足10。全部未通过、活性为负和不可计算候选也保存在完整源表。',supp=True)
h('4.9　供者重采样的全部关系与区间')
d=de('bootstrap_effect_summary');rows=[]
for axis,adj,label in [('IL18','B','IL18–IKZF2'),('TNF','B','TNF–STAT1'),('TNF','B_Treg','TNF–STAT1+Treg比例')]:
 for mod in ['RNA19','RNA17','AUC17','AUC29','AUC43']:
  rr=[d[(d.scheme=='stratified')&(d.axis==axis)&(d.adjustment==adj)&(d.modality==mod)&(d.edge==edge)].iloc[0] for edge in ['ligand_TF','TF_phenotype','ligand_phenotype']]
  rows.append([label,mod]+[f'{n(r.original_effect,3)} [{n(r.CI95_low,3)}, {n(r.CI95_high,3)}]' for r in rr]+[str(int(rr[0].valid))+'/5000'])
table('分层重采样的原效应与95%区间',['关系/模型','读数','配体–TF','TF–效应程序','配体–效应程序','有效次数'],rows,[3.0,1.5,3.55,3.55,3.6,2],note='表内为ρ[95%CI]；效应程序为PRF1/IFNG/NKG7。基础关联均控制疾病、来源及两端秩深度，标注Treg者额外控制其比例。AUC数字为种子。未生成新P/BH，不设置空FDR列。模型配对差和未分层失败记录单独保存在原始重采样表。',supp=True)
h('4.10　小鼠评分、调控子转移及历史 MRE11 对照')
d=df('mouse_UCell_score_group_tests');rows=[]
for _,r in d.iterrows():rows.append([sc(r['set']),r.module,n(r.difference),n(r.PValue),n(r.FDR),'探索性低细胞' if r['set']=='Treg' else '主要范围'])
table('小鼠UCell的全部疾病效应与精确置换结果',['范围','面板','硅尘−对照','P','BH','层级'],rows,[2.0,3,3.3,2.8,2.8,3.3],note='每组3文库，20种3+3分配；所有方向保留。各面板成员不同于人类手工面板，不能直接比较绝对评分。',supp=True)
fig('figS1_transfer.png','冻结调控子在独立小鼠时间点的评分方向','A、B分别为ActCD4的IKZF2与STAT1；每条线对应主队列的一个训练种子，横轴为独立7天/56天取材时间。每条件每时间点仅1文库，曲线仅连接两个描述性差值，不代表同一个体轨迹；不计算推断P或FDR。全部成员在接收数据特征空间可映射，不等同所有基因在每个细胞中检出。',supp=True)
rows=[]
for _,r in df('MRE11_human_historical').iterrows():rows.append([sc(r.analysis),r.module,int(r.subject_n),f'{n(r.spearman_rho)} / {n(r.spearman_p)}',triple(r,'adjusted_estimate','adjusted_p','fdr_bh')])
table('人类历史MRE11相关结果，独立于候选TF分析',['范围','面板','n','Spearman ρ / P','调整系数 / P / BH'],rows,[2.3,2.3,1,5.2,6.4],note='调整模型为面板评分~来源+疾病+MRE11检出比例+标准化log10(中位UMI+1)。MRE11结果未被重新命名为IKZF2或STAT1；小鼠历史MRE11的128项聚合相关记录保存在服务器，未据其筛选本报告候选。',supp=True)
h('4.11　数据资源与可复现范围')
textblock('''正式资源包含人类/小鼠基序数据库、NicheNet人类经典及小鼠v2先验、MSigDB来源的GO/KEGG成员和小鼠GO/KEGG映射。各资源源文件、哈希、物种、成员及映射以服务器档案为准。报告制作读取聚合结果，新增的程序只执行筛选、汇总、作图及排版，没有新增科学模型或显著性搜索。
完整档案包含：人类供者与亚群覆盖、504项程序疾病模型、全部人类通讯候选及状态、完整靶集/背景/配体–靶成员、小鼠方向与通讯源表、12个网络的统计及完整成员、供者重采样的效应/差值/失败诊断、GSEA/GSVA完整输出与敏感性，以及图表重建代码。个体数据与原始矩阵留在项目服务器；Library交付仅为报告PDF和可编辑文档。''')
