from pathlib import Path
import json,html,shutil,hashlib,math,zipfile
import pandas as pd
from PIL import Image as PILImage
from reportlab.pdfgen import canvas
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Image,Table,TableStyle,PageBreak,KeepTogether
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.enums import TA_CENTER
L=Path(r'C:\Users\31937\lung_cd4_plan1_20261015');O=L/'complete_report_20261015_v7';R=lambda f:pd.read_csv(L/f,sep='\t');pages=[];cur=None
pdfmetrics.registerFont(TTFont('CN',r'C:\Windows\Fonts\simsun.ttc',subfontIndex=0));pdfmetrics.registerFont(TTFont('CNB',r'C:\Windows\Fonts\simhei.ttf'))
def page(title):
 global cur
 cur={'title':title,'blocks':[]};pages.append(cur)
def p(text):cur['blocks'].append({'type':'p','text':text})
def sub(text):cur['blocks'].append({'type':'sub','text':text})
def tab(headers,rows,widths=None):cur['blocks'].append({'type':'table','headers':headers,'rows':[[str(x) for x in row] for row in rows],'widths':widths})
def fig(src,caption,height=300):
 src=L/src;name=src.parent.parent.name+'_'+src.name;shutil.copy2(src,O/'assets'/name);cur['blocks'].append({'type':'fig','src':'assets/'+name,'caption':caption,'height':height})
def num(x):
 if pd.isna(x):return '未检验'
 return f'{x:.3g}'
page('矽肺研究中巨噬细胞与CD4 T细胞表型关联分析')
p('完整修订报告  2026年10月15日  v7')
p('基准为正文核对截至2026年10月6日、于10月7日发布的复审R1R2报告。本版保留巨噬细胞信号与CD4表型的研究框架，重建候选选择依据，调整人类证据的疾病范围，并整合最新的真实时间、调控、通讯与患者定位结果。')
sub('摘要')
p('目的  在矽肺背景下筛选可供基础实验验证的巨噬细胞信号、CD4调控候选及功能读数，并回答候选从何而来、如何选择实验靶点以及煤工尘肺证据如何使用。')
p('方法  以GSE314943硅尘小鼠为发现队列，复用冻结差异表并恢复全部非空样本；分别比较CD4亚群的死亡、线粒体和自噬相关转录程序。以GSE274372提供真实时间支持，GSE183682提供方向补充。GSE174725用于有限的人类矽肺病例对照支持，GSE296682仅用于患者细胞定位。统计单位为独立生物学样本。')
p('结果  历史六候选在ActCD4中仍均为名义差异基因，其中IL1R2通过原BH校正。STAT1与GSDMD转录表达在Core、ActCD4、Treg中呈名义正相关，固定组成分析保留该联系；IKZF2与Core凋亡程序的关联则受亚群组成影响。CXCL16–CXCR6获得小鼠来源差异和两位患者的细胞定位支持，但尚未连接为直接STAT1驱动机制。患者Core中STAT1/GSDMD较低的方向在留一供者及质控敏感性后仍保留。')
p('结论  当前最适合分层验证的是CXCL16–CXCR6募集或滞留相关线索，以及CD4的STAT1–GSDMD相关转录状态。历史六基因继续作为差异响应候选，拟时序和真实时间作为辅助证据。完整的巨噬细胞配体至STAT1至焦亡因果链仍需功能实验验证。')
sub('阅读顺序')
p('第1至3节交代客户问题、数据和方法；第4至16节呈现筛选及新证据；第17至18节给出实验决策；第19至21节集中说明边界、旧稿变更及补充材料；末附参考文献。')
page('1 客户问题已有可追溯的对应答案')
sub('基因不是由火山图自动产生')
p('火山图展示每个基因的效应量和统计证据，本身不是选靶算法。历史六候选IL1R2、PENK、TNFRSF9、LAG3、TTC39C、FOXP1来自冻结差异筛选及后续综合判断。原有资料未提供一个能唯一自动生成这六个基因的全交集规则，本版不追溯改写为自动交集。第4节补充筛选数量、火山图、六基因效应与原始P及BH。')
p('IKZF2和STAT1是后续调控假设候选，不属于历史六基因，不能描述为从六基因自动筛出的两个。新版将其与其他可评估TF置于完整检验家族中评价，并分别保留RNA、活性和表型证据。')
sub('实验对象按作用位置选择')
p('配体或受体决定研究哪种细胞间作用，TF决定验证哪种接收端调控，GSDMD等执行分子决定测量哪类细胞死亡，历史六基因主要作为疾病响应指标。各角色不能凭一张排名表互相替代。第17至18节给出优先级、进入实验的依据和结果分支。')
sub('旧煤工尘肺资料不是纯SiO2患者亚组')
p('旧稿的人类主分析包含10名对照和9名煤工尘肺患者，来自GSE135893、GSE227136及GSE288887。现有核实记录不能把其中病例重新认定为经个体暴露确认的纯SiO2矽肺亚组。旧稿是煤工尘肺患者与硅尘小鼠的跨队列综合解释，并非已建立纯矽肺患者主发现。')
p('新版不将这些煤工尘肺细胞与GSE174725矽肺患者合并拟合疾病效应；GSE288887退出主要发现和选靶。旧人类SCENIC、GSEA等结果归入历史证据，不再改称矽肺患者验证。原报告保留在baseline目录，供核对。')
page('2 数据分工与样本单位明确区分')
tab(['数据集','实际分析范围','本版角色'],[
['GSE314943','硅尘20周；3 saline与3 silica小鼠文库','主发现；Core与Naive、ActCD4、Treg'],['GSE274372','D0、D7、D28、D56各3文库；共12','真实时间支持；TLike不等于纯CD4'],['GSE183682','D7、D56每条件每时间1只小鼠','仅描述方向，不作生物学重复推断'],['GSE174725','2暴露对照、3矽肺病例','有限病例对照支持；按复核注释'],['GSE296682','2位矽肺患者；4579个所选细胞','患者定位；无疾病对照检验'],['GSE135893 / GSE227136 / GSE288887','旧稿人类组合队列；10对照、9 CWP','仅历史结果说明，不参与当前选靶']], [110,180,190])
p('所有效应估计以独立小鼠或人类供者为单位。细胞数表示测量覆盖，不作为统计样本量；三个SCENIC随机种子仅用于算法敏感性。[4,6]Core是所选CD4总集合，亚群分别分析，避免把组成变化直接解释为亚群内变化。')
p('本版采用所有非空且可评估的样本和状态，不设至少20或10个细胞的硬门槛。小鼠Treg的13细胞对照保留；GSE274372的8细胞D0文库保留。零细胞或缺失特征保持不可评估，不补零表达或伪造P值。')
p('人类GSE174725的Core各供者细胞数为21、151、50、5、2725；各亚群的有效供者数按实际非空覆盖记录。BALF_AM病例仅1位有细胞，不用单核细胞替代以补足巨噬细胞病例样本量。')
page('3 分析路径以疾病变化和证据分层为主')
sub('可复现的分析顺序')
p('第一步：冻结输入、样本名单和历史差异表，恢复历史筛选数量及候选位置。第二步：按供者汇总各CD4亚群程序评分，比较疾病相对对照的变化、效应及一致性。第三步：评价差异基因、同状态轨迹及真实时间支持；三类证据分别输出。第四步：连接发送端表达、接收端受体、配体靶基因预测、TF活性及执行分子。第五步：在人类矽肺队列中分别检查疾病方向和细胞定位。')
sub('差异和相关检验分别定义')
p('差异表达复用经核对的伪批量edgeR结果及原全基因BH家族。程序疾病比较使用样本均值和精确组别置换；3对3的双侧精确均值差检验仅有20种分配，最小可达P通常为0.10。此类P不能与借助全基因信息估计离散度的edgeR P值简单等同。')
p('主要样本关联采用6个小鼠的Spearman精确置换，枚举720种配对；控制疾病组别的敏感性使用组内36种置换。另记录留一供者方向和固定亚群组成结果。TF主分析使用种子17的去重叠调控子活性，29和43为敏感性。')
p('默认以原始P<0.05筛选名义显著结果，同时保留效应量、完整P及BH。BH在预先定义的同类假设家族内计算，不混合差异基因、通路、相关性或外部验证，也不只对选中的两个TF重新校正。配体活性预测的AUPR及排名不是生物学效应P值。[5]')
p('历史六候选、调控候选和补充交集候选分别记录来源。创新性影响实验优先级，但不用于改变统计阈值或强制保留某个基因。')
page('4 历史筛选可以恢复但不能改写为自动交集')
p('使用20260914冻结差异表恢复历史流程。Full为原完整分析，NoBG为当时排除背景基因后的差异分析分支；NoBG不等于本轮逐细胞环境RNA校正。历史效应门槛|log2FC|>0.5仅用于复原旧流程，不作为新增候选统一准入条件。')
c=R('results/tables/historical_screening_counts.tsv');c=c[(c.arm=='Full')&c.scope.isin(['Core','ActCD4'])]
tab(['范围','受检基因','P<0.05','再加|logFC|>0.5','原BH<0.05'],[[z.scope,z.tested,z.nominal_P05,z.nominal_and_absLFC05,z.original_BH05] for z in c.itertuples()],[70,90,95,130,95])
fig('results/figures/historical_Core_ActCD4_volcano.png','图1  冻结Core和ActCD4差异结果。图中标出历史六候选及Full/NoBG分支；火山图展示来源及统计位置，不把图形本身当成候选选择算法。完整输入和原校正家族见补充表。',350)
p('复原得到的是差异候选池，而不是唯一的六基因生成规则。六候选的最终保留属于历史综合判断；本版将这一点明示，并给出足以独立重新排序的完整统计证据。')
page('5 历史六候选保留差异证据并明确角色')
a=R('results/tables/historical_six_DE_all_scopes.tsv');a=a[(a.arm=='Full')&a.scope.isin(['Core','ActCD4'])]
tab(['范围','基因','log2FC','原始P','原BH'],[[z.scope,z.gene_upper,num(z.logFC),num(z.PValue),num(z.FDR)] for z in a.itertuples()],[70,100,90,110,110])
p('ActCD4中六候选均达到原始P<0.05，IL1R2通过原BH<0.05；Core中五个上调候选通过原BH，FOXP1为名义下调。统计证据没有因为重新组织叙事而被删除。')
fig('results/figures/historical_six_evidence.png','图2  六候选的差异效应及统计证据。完整亚群、原家族和未受检记录随补充材料保留。',215)
p('六候选首先作为疾病响应和状态验证指标。IL1R2适合优先复核表达；PENK、TNFRSF9、LAG3反映不同响应线索；FOXP1为下调对照候选；TTC39C同时进入后续时间交集，但其阶段依赖方向仍需解决。以上角色不等同于已确立的致病驱动。')
page('6 亚群表型以疾病变化而非绝对最高评分选择')
p('对Naive、ActCD4和Treg分别计算供者平均程序，再比较疾病相对对照的变化。Core PyroExec差值+0.0564，P=0.20、BH=0.30；ActCD4差值+0.0516，P=0.50、BH=1；Treg差值+0.0396，P=0.60、BH=1。当前不支持把任一亚群宣布为疾病中显著焦亡最强的亚群。')
fig('customer_sequence_v2/figures/01_death_programs.png','图3  死亡相关程序的样本级疾病变化。每个点代表样本，非细胞重复。完整结果包含凋亡、炎症小体和焦亡相关评分及背景敏感性。',420)
p('PyroExec由PYCARD、CASP1、GSDMD、IL1B、IL18构成，是包含执行成分与炎症因子的转录程序，不是裂解型蛋白或细胞死亡测量。后续逐基因分析补充GSDMD证据，并未替换该原始评分。')
page('7 线粒体与自噬程序提示状态差异')
p('Core MitoROS变化+0.0579，P=0.10、BH=0.50；Treg变化+0.0637，P=0.10、BH=1。Naive Mitophagy变化-0.0485，P=0.60；Treg变化+0.0371，P=0.40。亚群间方向不同，不能通过绝对分数最高直接选择机制细胞。')
fig('customer_sequence_v2/figures/02_mito_autophagy_programs.png','图4  OXPHOS、线粒体应激、MitoROS、线粒体自噬及GO自噬的样本级比较。不同程序单位和基因组成不作为同一绝对强度尺度。',310)
p('广义Autophagy_GO已单独纳入表型比较，但未将其混入原去重叠TF程序关联的主要终点；原调控子排除列表不完全覆盖新增GO集合。程序转录变化与线粒体膜电位、ROS生成、自噬通量及死亡执行在实验上分别检测。')
p('亚群优先级综合疾病表达、受体覆盖、调控关联及实验可分选性。ActCD4适合作为CXCR6及差异响应重点，Treg保留STAT1/GSDMD关联线索；Naive可作为状态对照。此选择不声称已证明某亚群焦亡最强。')
page('8 交集分析已补齐但仅作为辅助证据')
p('作为客户提出路径的可复核实现，曾在Naive和ActCD4中依次要求：疾病DE原始P<0.05；同状态轨迹P<0.05；轨迹平均斜率与疾病效应同向；真实时间总体检验P<0.05；至少一个D7/D28/D56相对D0对比名义显著且与疾病效应同向。所有阶段成员名单和未通过原因均保留。')
s=R('customer_sequence_v2/04_candidate_intersection/stepwise_gene_counts.tsv');pv=s.pivot(index='stage',columns='state',values='genes');order=['tested_disease','DE_P05','DE_and_PT_P05','DE_and_PT_direction','plus_time_omnibus','plus_time_direction'];labels=['疾病受检','疾病P<0.05','再满足轨迹P<0.05','再满足轨迹同向','再满足时间总体P<0.05','再满足至少一个同向时间对比']
tab(['步骤','Naive基因数','ActCD4基因数'],[[label,int(pv.loc[k,'Naive']),int(pv.loc[k,'ActCD4'])] for k,label in zip(order,labels)],[300,90,90])
p('最终Naive保留AHNAK；ActCD4保留MDFIC、TTC39C、ST6GALNAC3。该结果不等同于历史六基因来源，也不能据此强制删除IKZF2或STAT1等其他角色的候选。')
p('综合路径不把状态内拟时序设为必要条件。细胞轨迹表示状态排序，真实时间表示取样日龄或暴露时长；二者提供不同信息。轨迹以样本级斜率推断，不把细胞当独立重复；路径和根选择保留敏感性，不以最有利轨迹替代主要结果。')
fig('customer_sequence_v2/figures/04_candidate_screening_flow.png','图5  辅助交集的明确规则和逐步数量。该流程用于检验客户提出的路径，不追溯重写历史选基因过程。',240)
page('9 真实时间支持与外部方向必须分开解释')
p('GSE274372包括D0、D7、D28、D56各3文库，按实际时间分类建模，全部12个非空TLike文库纳入。时间总体检验与各时间对比各自使用完整受检基因BH；此数据中的TLike未确认为纯CD4。')
fig('customer_sequence_v2/figures/03B_candidate_real_time.png','图6  四个辅助候选的真实时间变化，使用独立样本伪批量表达；不是将时间点重新拟合为细胞拟时序。',210)
a=R('customer_sequence_v2/04_candidate_intersection/candidate_evidence_cards.tsv')
tab(['候选','疾病状态','DE P / BH','时间总体P / BH','同向时间'],[[z.gene,z.mouse_state,num(z.DE_P)+' / '+num(z.DE_BH),num(z.time_omnibus_P)+' / '+num(z.time_omnibus_BH),z.nominal_same_direction_days.replace('_D0','')] for z in a.itertuples()],[90,65,110,115,100])
p('GSE314943为20周发现终点，GSE183682为D7/D56，后者每条件每时间只有1只小鼠。AHNAK、ST6GALNAC3的外部方向相对一致；TTC39C的D7口径敏感且D56下降仍存在，MDFIC的外部下降也未消失。不能将不同疾病阶段简单合并或把无重复的方向数据写成显著验证。')
p('时间序列缺少各晚期时间同步saline对照，年龄与暴露持续时间难以分开。该边界用于解释时间证据，不否定已观察到的样本表达方向。')
page('10 创新性决定研究问题而非统计阈值')
tab(['候选组','已有基础与当前定位','实验优先级'],[
['历史六基因','差异响应来源可追溯；不等同新的调控驱动','表达及状态验证保留'],['IKZF2','T细胞状态与调节相关；Core组成敏感','状态分层后评价'],['STAT1','干扰素及矽肺相关研究已有；本项目有接收端关联','优先接收端候选'],['ST6GALNAC3','酶学、肿瘤、IPF表达与CD4相关研究已有','特定硅尘ActCD4功能可探索'],['AHNAK','CD4钙流与NFAT机制已有；人类方向存在差异','机制基础明确的并行备选'],['TTC39C','部分其他模型已有功能研究；跨阶段方向不一致','先解决重复性'],['MDFIC','肺纤维化及其他细胞背景已有研究；外部方向矛盾','当前降低优先级']], [95,280,105])
p('10/6的创新性考虑可以保留为候选排序因素，但不能代替差异来源、样本一致性和实验可验证性。本次不以“文献少”自动认定新靶点，也不因基因已被研究就删除项目中真实存在的统计信号。')
p('较有价值的问题是：在明确的硅尘阶段和CD4状态中，巨噬细胞相关信号如何影响迁移、接收端调控及GSDMD相关表型。CXCL16–CXCR6轴本身及STAT1、GSDMD在其他背景均已有研究，[2,7–11]创新需落在细胞情境及干预关系。文献证据和原始检索记录保留在补充材料，不作“首次发现”的断言。')
page('11 STAT1与IKZF2显示不同的接收端证据')
p('去重叠TF活性与程序关联中，Core STAT1–MitoROS：rho=0.9429，P=0.01667、BH=0.3678；Treg：rho=0.8857，P=0.03333、BH=0.3790。Core IKZF2–ApoExec：rho=1，P=0.002778、BH=0.2743。以上均为名义关联，原完整TF家族未通过BH<0.05。')
fig('integrated_v4/figures/TF_programs.png','图7  两候选TF的各状态程序关联。主种子17；星号为名义P<0.05。TF活性与RNA表达分别解释。',325)
p('STAT1在Treg的疾病RNA差异log2FC=2.1506，P=2.5823e-5、BH=0.005589；这项差异通过原BH，不意味着所有STAT1活性关联也通过校正。两个TF与原PyroExec总评分均没有名义显著关联。')
page('12 固定组成分析收窄IKZF2的机制解释')
fig('integrated_v4/figures/composition_sensitivity.png','图8  Core关联的固定亚群组成敏感性。以全部6样本平均亚群比例为统一权重重新组合状态读数。',280)
p('Core IKZF2–ApoExec从rho=1变为-0.7714，P=0.1028、BH=0.5688。原正相关不能作为IKZF2在相同细胞状态内促进凋亡的证据。IKZF2仍保留为身份或组成相关调控线索，不因这项敏感性结果抹去原差异及状态信息。')
p('Core STAT1–MitoROS固定组成后rho=0.8857，P=0.03333、BH=0.4467，方向保留；进一步控制疾病组别后P=0.1667，提示关联仍受到疾病组间差异影响。TF活性由冻结网络推断，不等于p-STAT1蛋白活化。')
p('去重叠评分排除了预定面板、候选及TF自身相关重叠，但表达匹配的对照基因仍可能重叠。新加入的GBP、IRF等分支若未从原regulon中剔除，其关联仅为非独立探索，不能作为机制链的主要支撑。')
page('13 GSDMD分解分析补出了执行分子线索')
fig('STAT1_bridge_v5/figures/STAT1_GSDMD_samples.png','图9  STAT1活性与GSDMD表达的样本级关联。每点为小鼠，原PyroExec总评分不被替换。',225)
tab(['范围','rho','P','BH','疾病RNA差异'],[['Core','0.9429','0.01667','0.06667','GSDMD logFC +0.7477；P 0.03455；BH 0.2790'],['ActCD4','0.9429','0.01667','0.1000','GSDMD logFC +0.8458；P 0.04606；BH 0.6934'],['Treg','1.0000','0.002778','0.05556','GSDMD logFC +0.7793；P 0.1242；BH 0.4697']],[65,65,75,75,200])
p('固定组成后Core STAT1–GSDMD rho=1，P=0.002778、BH=0.02778，在本次定向基因家族内通过BH；控制疾病组别后P=0.05556、BH=0.5556。此校正后显著属于事后补充的固定组成分析，不改称原主分析全面显著。三个算法种子方向一致。')
p('Core CASP1也有名义关联，但固定组成后不再名义显著；GSDME主要种子未达到P<0.05，IRF1没有相应的一致上调与关联。鼠Casp4即历史caspase-11命名对应[13]，不将缺失的CASP11别名重复计为另一基因。')
p('本节支持转录层面的焦亡相关候选，而非已确认焦亡发生。GSDMD蛋白活化、膜定位和裂解性死亡仍需独立测定。[7,8]')
page('14 CXCL16和原配体提供不同的上游依据')
a=R('upstream_human_v6/tables/mouse_upstream_integrated.tsv');a=a[(a.scope=='ActCD4')&(a.model=='Pooled')&a.x.isin(['TNF','IL15','IL18','C3','CXCL16'])]
tab(['来源','配体','来源logFC','来源P / BH','与STAT1 rho','关联P / BH'],[[z.sender,z.x,num(z.ligand_logFC),num(z.ligand_P)+' / '+num(z.ligand_BH),num(z.rho),num(z.P)+' / '+num(z.BH)] for z in a.itertuples()],[65,55,65,105,75,115])
p('TNF、IL15在NicheNet的接收者上调靶基因预测中排序靠前，且有受体覆盖；但来源疾病差异和与STAT1活性的配对关联尚不足，预测排名不能替代缺失边。IL18来源下降及与ActCD4 IKZF2的负向关联应按负向或补偿假设解释，不写成促焦亡正向链。')
p('CXCL16在AMlike中名义上调，但AMlike CXCL16–ActCD4 STAT1的P=1；AllMac的对应关联P=0.002778，而AllMac CXCL16疾病差异P=0.1487。不能跨来源拼接最有利的两段。病例ActCD4的CXCR6检出率28.17%、30.00%、53.78%，支持将该轴作为募集或滞留相关候选。')
p('C3的AllMac差异通过BH，但ActCD4/Treg的C3AR1覆盖仅1/3病例；广义数据库受体不能替代具体活化产物及特异受体。IL12A单亚基差异也不代表完整IL-12，IL27完整亚基和受体跨样本覆盖不稳。')
page('15 矽肺患者支持保留关键方向差异')
p('GSE174725使用复核后的矽肺注释，未混入旧CWP活性。Core STAT1差值-0.9190 CP10k，P=0.10；GSDMD差值-0.5902，P=0.30。本轮31基因家族BH分别0.8667和1。逐个剔除供者后两者差值始终为负；排除既有doublet和混合标记后方向近似不变。')
fig('upstream_human_v6/figures/human_composition_decomposition.png','图10  等供者Core差值的描述性组成分解。状态内部分不是因果效应，缺失状态不补造表达。',245)
p('STAT1总差值中组成部分约-0.0476、状态内部分约-0.8714；GSDMD分别约-0.0862和-0.5040。关键反向结果主要不能由亚群比例解释掉。MemCD4的CXCR6和GSDMD方向向上，P均0.3333；该状态仅2对照与2病例有细胞，不用其替代Core结果。')
p('BALF_AM仅1个病例可评估，故不计算人类巨噬细胞的疾病推断或配体至CD4相关检验。GSE174725目前提供有限、部分相反的患者证据，不能写成已完成同方向机制验证。')
page('16 两位患者补强CXCL16与CXCR6细胞定位')
fig('upstream_human_v6/figures/patient_localization.png','图11  GSE296682两位矽肺患者的作者细胞标签内表达定位。分母为对应细胞数，统计生物学样本量仍为2位患者。',230)
p('肺泡巨噬细胞CXCL16阳性264/323及601/782，CD4 effector T cell的CXCR6阳性5/33及2/27；这支持发送端与接收端定位。作者CD4 effector标签不直接等同于小鼠ActCD4。两患者中各有1个CXCR6阳性效应T细胞同时检出至少2个T细胞标记，且CXCR6阳性细胞未达到预设髓系混合标记规则。')
p('GSDMD在CD4 effector中检出4/33和2/27，但STAT1/GSDMD同细胞共检出为3/33和0/27，前者3个均有髓系混合标记；CXCR6/STAT1/GSDMD三重共检出均为0。这些结果不能构成完整同细胞机制支持，未检出也不等同于蛋白绝对不存在。')
p('两位患者CD4标签中C3AR1未检出；IL27完整受体在效应CD4中为0/33和4/27，供者覆盖有限。该队列仅患者定位，无对照、不作疾病差异检验，也不计算两供者相关P值。')
page('17 基础实验候选按角色和证据选择')
tab(['角色','当前优先候选','选择依据','首轮要回答的问题'],[
['细胞间作用','CXCL16–CXCR6','鼠来源名义差异；患者配体和受体定位','巨噬细胞信号是否改变CD4迁移或滞留'],['接收端调控','STAT1','Treg RNA差异及GSDMD关联','是否存在p-STAT1变化及功能依赖'],['死亡执行读数','GSDMD与相应caspase','逐基因转录关联','CD4是否有GSDMD依赖的裂解性死亡'],['差异响应指标','IL1R2优先；PENK等并行','历史差异可追溯','表达和细胞状态是否可重复'],['状态相关调控','IKZF2','与身份或组成有关','相同状态内是否仍有调控效应'],['创新性并行候选','ST6GALNAC3；AHNAK备选','辅助交集及背景特异研究空间','先验证表达、状态及作用方向'],['暂缓深入机制','MDFIC；TTC39C阶段问题','外部方向或稳健性尚有矛盾','先解决跨队列或阶段重复性']], [75,100,150,155])
p('不建议一次把全部候选作为同等级敲低靶点。优先组合可以是CXCL16/CXCR6的来源与迁移验证，加上STAT1/GSDMD的接收端和执行验证，使用历史六基因中的明确响应指标作配套读数。两段是否能连接由实验决定，不预设阳性结论。')
p('IL1R2的差异证据较强不等于其为上游驱动；STAT1的已有研究较多不等于不能研究特定CD4情境；ST6GALNAC3的探索空间也不意味着优先于可靠的表型和来源确认。')
page('18 实验按身份 表型 细胞间作用和因果次序推进')
sub('第一阶段 先确认细胞身份与功能表型')
p('在独立样本中明确分选或定位CD3阳性CD4阳性亚群，同时检查髓系混入。复核CXCR6、STAT1及GSDMD的RNA和蛋白。线粒体功能与死亡执行分别测定：线粒体ROS、膜电位、细胞存活，以及GSDMD活化、相应caspase活化和膜通透性等互补指标。仅RNA、LDH或炎症因子升高不作为单独的焦亡结论。')
sub('第二阶段 确认巨噬细胞介导而非颗粒直接作用')
p('使用匹配的巨噬细胞条件培养或共培养，与对照巨噬细胞、直接硅尘作用及残留颗粒相关对照比较；记录CD4亚群比例、迁移及亚群内表型。先判断表型是否出现，再检测上清配体蛋白与接收端p-STAT1的时间先后。')
sub('第三阶段 分开验证募集与死亡机制')
p('CXCL16/CXCR6干预优先回答迁移或滞留问题。若只改变迁移，不写成直接促进STAT1或焦亡。接收端STAT1干预需同时观察GSDMD活化和存活，以区分一般激活或毒性影响；GSDMD干预用于判断执行依赖。必要时采用适当救援与正交手段。')
sub('判读分支')
tab(['实际结果','对应结论和下一步'],[['CXCR6干预改变迁移而无死亡变化','保留募集机制；暂停直接焦亡叙事'],['STAT1干预改变转录但无GSDMD活化','定位为转录调节；不推定焦亡执行'],['GSDMD干预降低裂解性死亡','建立执行依赖，再定位上游'],['清理背景或纯化后信号消失','回到细胞来源，不为保留靶点调整阈值'],['不同阶段或亚群方向相反','按阶段或状态建立假设，独立重复后再合并']], [215,265])
p('独立生物学重复和样本量根据预实验效应与变异规划；技术重复不代替供者。状态内拟时序不是上述功能验证的必要条件。')
page('19 讨论集中说明关键矛盾与推断边界')
sub('环境RNA是当前优先解决的测量问题')
p('细胞身份敏感性分析发现冻结CD4矩阵广泛检出肺上皮标记。Saline_2的424个CD4细胞均检出SFTPC与SCGB1A1，多数也有T细胞标记；该样本GSDMD阳性31个细胞，每个仅1 UMI。按至少2个异源标记严格排除后该样本全部排空，无法维持六样本复核。不能把这些细胞全部判为双细胞，也不能宣称已排除背景RNA影响。原矩阵、注释及结果均未覆盖。')
p('需要以原始液滴或合适的已校正对象复核背景贡献，并以纯化CD4蛋白验证。上述风险限定转录机制解释，不静默撤销已经计算出的效应或调整阈值来追求保留。')
sub('病程和物种差异不能用交集消除')
p('20周硅尘发现、D7/D56方向支持及人类矽肺不是同一生物过程的重复测量。GSE274372还缺少同期晚期对照。患者Core的反向结果在供者、质量和组成检查后仍存在，必须保留；不能单独挑选正向MemCD4作为完整人类验证。')
sub('相关与执行之间尚有实验证据缺口')
p('细胞通讯先验、配体活性排序、TF推断活性与RNA关联各有不同假设。尚未证明具体巨噬配体导致STAT1活化并引发CD4焦亡。当前结果支持分段选择实验，而不是已经建立临床可用靶点。线粒体损伤、凋亡、焦亡与自噬不能在叙事中视为同义或固定顺序。')
page('20 原报告的保留 修正和证据迁移')
tab(['10月6日内容对应位置','本版处理','理由'],[['摘要及1.2人类19供者CWP主分析','改为硅尘主发现与有限矽肺患者支持','疾病标签与研究目标对齐'],['2.1至2.3患者图谱和通讯','旧CWP结果保留归档；不作当前患者验证','未确认纯SiO2病例子集'],['2.4小鼠通讯','保留计算框架；补来源、受体与配对证据','区分预测排名和机制链'],['2.5至2.9 IKZF2及STAT1','保留独立候选；补全TF家族及组成敏感性','避免跨读数和组成混杂'],['2.10至2.12人类GSEA与死亡评分','旧CWP结果迁入历史证据说明','不替代矽肺患者结果'],['2.13小鼠功能程序','恢复所有非空样本并按亚群疾病变化呈现','取消最低细胞数门槛'],['旧稿未充分解释的基因来源','新增火山图、六候选统计及辅助交集规则','回答客户选基因问题'],['3.6后续研究及3.7结论','新增角色优先级、实验判读分支','衔接基础研究']], [145,170,165])
p('原发布PDF未修改，保存在baseline/20261006_content_20261007_release.pdf。新稿不将旧CWP的IKZF2、STAT1活性或氧化磷酸化富集改名后移作矽肺证据。保留的主要是研究问题、分析方法、可复用小鼠结果和历史候选记录；需要改变的疾病结论及核心矛盾已在正文明确呈现。')
p('历史六候选的差异来源与后续综合判断得到补全；基础实验可依据不同角色和证据阶段选择；煤工尘肺与矽肺数据的分析边界已经单列。以上三项构成本次修订对客户问题的对应回答。')
page('21 补充统计与可复现材料')
tab(['目录','内容'],[['supplement/results','冻结筛选数量、六基因全范围统计、原火山图输入和复用清单'],['supplement/customer_sequence_v2','亚群程序、样本级轨迹、12文库真实时间及全部交集成员'],['supplement/novelty_direction_v3','创新性文献、模型身份核查、跨阶段方向诊断'],['supplement/integrated_v4','完整TF及RNA关联家族、来源配体、受体和组成敏感性'],['supplement/STAT1_bridge_v5','执行基因分解、患者方向、细胞身份和背景标记检查'],['supplement/upstream_human_v6','来源整合、患者定位、留一供者及组成分解'],['code','各阶段生成与统计脚本副本；大型原始矩阵保留服务器'],['provenance','基准身份、源文件哈希、报告构建和一致性核对记录']], [190,290])
p('补充目录保留阶段原始标签和统计家族。较早README或阶段结论属于当时检查状态，本版正文及相应最新版统计表作为当前解释依据；禁止将不同阶段的BH家族或样本门槛混用。缺失P/BH保留原样。')
p('代码依赖服务器已冻结输入、网络先验及软件环境；本包不是包含全部原始测序数据的一键运行镜像。未下载PRJNA1289989，未启动新的全量聚类或SCENIC训练。旧结果通过复制复用，新结果保存在独立目录。')
p('输出包含正式PDF、与正文一致的可编辑HTML源稿、图件、完整补充统计和代码。PDF采用嵌入中文字体并进行逐页渲染核查；新报告的正文数字按已核实表格生成或核对，所有文件哈希列入交付清单。')
page('22 参考文献与数据来源')
refs=[
('1','Cassel等，2008。NLRP3炎症小体与矽肺。','https://pmc.ncbi.nlm.nih.gov/articles/PMC2449360/'),
('2','You等，2024。硅尘相关CD4组织驻留记忆T细胞。','https://www.nature.com/articles/s42003-024-06662-z'),
('3','Ding等，2023。硅尘暴露巨噬细胞外泌体与T细胞平衡。','https://pubmed.ncbi.nlm.nih.gov/37918332/'),
('4','SCENIC方法，Aibar等，2017。','https://pubmed.ncbi.nlm.nih.gov/28991892/'),
('5','NicheNet方法，Browaeys等，2020。','https://www.nature.com/articles/s41592-019-0667-5'),
('6','单细胞差异分析中生物学重复的重要性，Squair等，2021。','https://www.nature.com/articles/s41467-021-25960-2'),
('7','矽肺GSDMD相关研究。','https://pmc.ncbi.nlm.nih.gov/articles/PMC9069405/'),
('8','矽肺Caspase-1/Gsdmd及Caspase-3/-8/Gsdme研究。','https://pmc.ncbi.nlm.nih.gov/articles/PMC9718385/'),
('9','STAT1蛋白修饰与矽肺，2024。','https://pmc.ncbi.nlm.nih.gov/articles/PMC11107463/'),
('10','IL27激活CD4 STAT1的原始实验研究。','https://pubmed.ncbi.nlm.nih.gov/15356135/'),
('11','CXCR6与肺部T细胞定位，CD8背景。','https://pmc.ncbi.nlm.nih.gov/articles/PMC6888981/'),
('12','小鼠发现队列20周与样本单位依据，原文图2A。','https://www.nature.com/articles/s41467-026-70193-w'),
('13','NCBI小鼠Casp4与历史caspase-11命名。','https://www.ncbi.nlm.nih.gov/gene/12363/')]
for n,t,url in refs:p(f'[{n}] {t} {url}')
p('数据记录：GSE314943、GSE274372、GSE183682、GSE174725、GSE296682均可在NCBI GEO按编号查询。患者定位采用已核实的官方GSE296682对应版本，来源和范围读取记录随补充材料保存。候选创新性具体文献及检索日志见supplement/novelty_direction_v3/literature。')
page('5 历史六候选的逐基因保留依据')
tab(['基因','当前可以复核的保留依据','基础实验定位'],[
['IL1R2','Core与ActCD4均上调且通过各自原BH；Full与NoBG方向一致','优先表达和状态验证指标'],
['PENK','Core和ActCD4上调；Core原BH显著，ActCD4名义显著；效应较大','复核低表达背景与蛋白后作为响应指标'],
['TNFRSF9','Core与ActCD4上调，Core原BH显著；Full与NoBG方向一致','活化相关响应候选，尚非已证实驱动'],
['LAG3','Core与ActCD4上调，Core原BH显著；Full与NoBG方向一致','抑制相关状态读数，不直接等同功能耗竭'],
['TTC39C','两范围名义上调；新增ActCD4轨迹和时间交集保留；外部阶段有反向','历史与新增筛选重合的候选，先处理阶段差异'],
['FOXP1','Core与ActCD4均名义下调；Full与NoBG方向一致；未通过原BH','下调状态对照候选，不能只保留上调基因']], [65,280,135])
p('以上是本次能够复核的逐基因保留理由，不声称恢复了当时所有人工判断及其权重。已有记录不足以唯一解释为什么在全部差异基因中最终只取这六个，因此本版同时开放完整候选池和统计表，而非把这六个包装为唯一最优解。')
sub('新版采用按角色的准入与排序规则')
p('差异响应基因：对应状态疾病效应原始P<0.05进入候选池，随后看效应、样本方向、背景敏感性和实验可测性。调控候选：以TF活性及表型关联为独立证据类型，保留完整TF家族BH，并检查种子和组成敏感性，不要求其必须来自历史六基因。')
p('通讯候选：先检查发送端和兼容受体在病例样本中实际检出，再分别报告来源差异、靶基因预测及同样本关联；只有预测支持的保留为预测候选，不能升级为已证实驱动。创新性用于说明具体研究问题和同等证据下的优先级，不量化为事后挑选的分数。')
p('辅助交集增加证据维度，不是所有角色的共同硬门槛。新报告的选择依据、未通过条件和现阶段实验角色分别列出，使后续研究者能够按研究目标重新选靶。')
pages.insert(6,pages.pop())
# Remove accidental typo without touching numeric evidence.
for pg in pages:
 for b in pg['blocks']:
  if b['type']=='p':b['text']=b['text'].replace('casapse-11','caspase-11')
(O/'report_content.json').write_text(json.dumps(pages,ensure_ascii=False,indent=2),encoding='utf-8')
styles={'body':ParagraphStyle('body',fontName='CN',fontSize=10.4,leading=16,spaceAfter=8,wordWrap='CJK'), 'title':ParagraphStyle('title',fontName='CNB',fontSize=17,leading=23,spaceAfter=14,textColor=colors.black),'sub':ParagraphStyle('sub',fontName='CNB',fontSize=12,leading=17,spaceBefore=6,spaceAfter=7,textColor=colors.black),'cap':ParagraphStyle('cap',fontName='CN',fontSize=8.5,leading=12,spaceAfter=9,wordWrap='CJK'),'cell':ParagraphStyle('cell',fontName='CN',fontSize=8.4,leading=12,wordWrap='CJK')}
def esc(x):return html.escape(str(x)).replace('\n','<br/>')
story=[];web=['<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>矽肺研究完整修订报告</title><style>body{font:16px/1.8 "Microsoft YaHei",serif;max-width:1050px;margin:35px auto;padding:0 24px;color:#111}h1,h2,h3{color:#000}section{margin:40px 0}table{border-collapse:collapse;width:100%;font-size:14px}td,th{border:1px solid #d9d9d9;padding:8px}th{background:#e6edf3}img{width:100%}.caption{font-size:13px}a{color:#17577e}@media print{section{break-before:page}}</style>']
for idx,pg in enumerate(pages):
 if idx:story.append(PageBreak())
 story.append(Paragraph(esc(pg['title']),styles['title']));web.append('<section><h1>'+esc(pg['title'])+'</h1>')
 for b in pg['blocks']:
  if b['type'] in ['p','sub']:
   story.append(Paragraph(esc(b['text']),styles['body' if b['type']=='p' else 'sub']));tag='p' if b['type']=='p' else 'h2';web.append(f'<{tag}>'+esc(b['text'])+f'</{tag}>')
  elif b['type']=='table':
   data=[[Paragraph(esc(z),styles['cell']) for z in row] for row in [b['headers']]+b['rows']];widths=b['widths'] or [480/len(b['headers'])]*len(b['headers']);t=Table(data,colWidths=widths,repeatRows=1,hAlign='LEFT');t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e5edf3')),('GRID',(0,0),(-1,-1),.4,colors.HexColor('#d9d9d9')),('VALIGN',(0,0),(-1,-1),'MIDDLE'),('LEFTPADDING',(0,0),(-1,-1),5),('RIGHTPADDING',(0,0),(-1,-1),5),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5)]));story.extend([t,Spacer(1,10)]);web.append('<table><tr>'+''.join('<th>'+esc(z)+'</th>' for z in b['headers'])+'</tr>'+''.join('<tr>'+''.join('<td>'+esc(z)+'</td>' for z in row)+'</tr>' for row in b['rows'])+'</table>')
  else:
   im=PILImage.open(O/b['src']);scale=min(480/im.width,b['height']/im.height);img=Image(str(O/b['src']),width=im.width*scale,height=im.height*scale);story.append(KeepTogether([img,Spacer(1,4),Paragraph(esc(b['caption']),styles['cap'])]));web.append('<img src="'+b['src']+'"><p class="caption">'+esc(b['caption'])+'</p>')
 web.append('</section>')
web.append('</html>');(O/'矽肺研究完整修订报告_20261015.html').write_text(''.join(web),encoding='utf-8')
class NumberedCanvas(canvas.Canvas):
 def __init__(self,*args,**kw):super().__init__(*args,**kw);self.saved=[]
 def showPage(self):self.saved.append(dict(self.__dict__));self._startPage()
 def save(self):
  count=len(self.saved)
  for state in self.saved:
   self.__dict__.update(state);self.setFont('CN',8);self.setFillColor(colors.HexColor('#555555'));self.drawString(54,26,'矽肺研究完整修订报告  20261015 v7');self.drawRightString(541,26,f'{self._pageNumber} / {count}');super().showPage()
  super().save()
doc=SimpleDocTemplate(str(O/'矽肺研究完整修订报告_20261015.pdf'),pagesize=(595.28,841.89),leftMargin=55,rightMargin=60,topMargin=45,bottomMargin=48,title='矽肺研究中巨噬细胞与CD4 T细胞表型关联分析',author='项目研究组');doc.build(story,canvasmaker=NumberedCanvas)
shutil.copy2(Path(__file__),O/'code/build_report_v7.py');print('BUILT',len(pages),'sections')