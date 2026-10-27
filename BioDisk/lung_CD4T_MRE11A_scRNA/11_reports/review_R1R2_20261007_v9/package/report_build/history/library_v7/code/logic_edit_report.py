from pathlib import Path
from lxml import html, etree
from collections import Counter
import base64, hashlib, json, re

B=Path(__file__).resolve().parent
D=B/'logic_revision_20261007'
source=B/'review_delivery/provenance/review_render_source.html'
tree=html.fromstring(source.read_text(encoding='utf-8'))
nodes=tree.xpath('//body/*')
assert len(nodes)==247
old=json.loads((D/'provenance/original_blocks.json').read_text(encoding='utf-8'))
new={
0:'矽肺研究中巨噬细胞与 CD4⁺T 细胞表型关联分析',
1:'基于煤工尘肺患者与小鼠硅尘模型的整合证据',
2:'数据及文献核对截至2026年10月6日',
5:'方法：人类CD4 RNA分析纳入19名供者（10名对照、9名煤工尘肺患者），以供者为统计单位比较细胞组成和转录程序，并按分析设置来源及深度协变量。巨噬细胞通讯结合配体–受体表达、NicheNet配体活性和同供者关联。候选RNA表达与SCENIC调控子活性分别评价，使用三个随机种子、去重叠基因及供者Bootstrap进行敏感性分析。GSEA、GSVA、UCell和表达匹配评分用于刻画功能程序。小鼠主队列每组3个文库，用于亚群、程序和信号方向的对照。',
11:'尘肺相关肺纤维化伴随持续的髓系炎症和适应性免疫变化。晶体硅尘模型的遗传学研究表明，NLRP3炎症小体参与肺部炎症及纤维化。[1] 巨噬细胞外泌体影响T细胞平衡的研究，为跨细胞信号提供了生物学依据。[5] 本研究结合患者细胞图谱、配体–靶基因预测和供者关联，考察巨噬信号与CD4 T细胞表型的联系。',
12:'IKZF2编码Helios，与调节性T细胞及活化的常规T细胞状态有关。[6] STAT1参与干扰素应答和辅助性T细胞分化，其作用具有细胞类型和组织环境依赖性。[2,8] 本研究将IKZF2与STAT1作为两个独立候选，分别评价RNA表达、调控子活性、程序关联及跨物种方向。',
13:'分析依次呈现免疫细胞全景与CD4表型、巨噬细胞通讯、两个独立候选及功能程序。细胞组成、供者平均表达、调控子AUC、通路富集、配体活性和相关系数按各自的分母与统计假设解释。三个计算种子用于评价算法敏感性，背景筛选用于评价细胞组成和测量背景对结果的影响。',
15:'人类CD4分析纳入19名供者，包括10名对照和9名煤工尘肺（CWP）患者。GSE135893提供2名对照，GSE227136提供8名对照与2名患者，GSE288887提供7名患者。统计分析以供者为单位；细胞组成和程序评分分别采用基础、来源校正及来源与深度校正模型，其他读数的模型见表S1。',
16:'肺细胞全景图谱覆盖20名供者；CD4分析覆盖其中19名供者。全体CD4（Primary）共6523个细胞，包含CD4state；核心CD4（Core）共5705个细胞，为排除该状态818个细胞后的集合。RNA、比例、程序及调控子分析按各自的细胞数与评分覆盖条件确定供者名单。完整细胞的调控子AUC可评分名单为17人，背景筛选后的覆盖见表1。',
17:'小鼠主要疾病对照来自GSE314943，Saline与Silica每组3份scRNA文库，分别分析核心CD4、Naive、ActCD4和Treg。GSE183682在7天与56天每条件各1个文库，用于固定调控子的跨队列方向比较；其结果报告描述性差值。小鼠分析以文库为单位，文库与动物个体的对应资料及其解释范围见讨论3.4。',
21:'CD4细胞按CD4-CCR7、CD4-GZMA、CD4-IL7R、Treg和CD4state注释分析。亚群比例分别以同供者全体CD4或核心CD4细胞总数为分母；程序评分在相应亚群内汇总。供者具有有效分母时保留比例零值，程序评分则要求该亚群具有实际细胞。按“有细胞”和“至少20个细胞”分别记录分析名单。',
22:'巨噬发送端包括肺泡巨噬细胞富集集合MacAM、单核/巨噬混合集合MonoMac及两者的并集MacWide。MacAM由簇0、2、11、13组成；MonoMac由簇3、7组成，C1Q、MARCO、SPP1与FCN1共同表达支持其混合特征。三个集合分别记录表达和供者覆盖，MacWide用于评估扩大来源范围后的结果。T/NK仅用于竞争性来源敏感性分析。',
24:'人类以供者、小鼠以文库为分析单位。[11] 比例及表达匹配程序的基础模型为Y~疾病；来源校正模型增加研究来源，来源与深度模型再增加相应细胞范围的log1p(平均原始UMI)。GSE227136来源内模型仅比较该来源的两组。固定供者名单上调整协变量与筛除背景较高细胞分别实施，并报告各自的样本覆盖。',
25:'各读数的深度项按其统计定义设置。供者均值RNA/AUC疾病模型使用log10(平均深度)，AUC组间参照模型使用log1p(平均深度)；调控子–程序偏相关对深度取秩后进行残差化。UCell质量校正使用标准化log10(中位UMI+1)及标准化线粒体比例。伪批量edgeR以文库大小和TMM构建偏置项。完整公式及各模型的统计输出条件见表S1。',
26:'各分析报告实际纳入人数、效应、P值及对应完整检验家族的BH结果。12程序疾病模型在细胞范围、背景筛选、模型和指标家族内跨程序校正；较窄程序面板的跨范围家族另行标识。调控子–程序关联按全部受检TF与程序及候选表达特征构成家族。图表中的“原BH”指该完整家族已计算的校正值，展示子集不另设家族。“未计算”与“不可估计”分别标注，具体原因见表S1；缺失统计不参与显著性判断。',
28:'按经核验的CD4接收亚群及巨噬来源建立通讯证据。接收者疾病差异以供者伪批量计数拟合基础或来源校正模型，保留完整原始文库偏置和TMM；分别以P<0.05或FDR<0.05，并结合|log₂FC|≥0.5，定义上调和下调靶集。两端供者须各有≥20个相应细胞；主表达门槛为至少两个病例供者中≥10%细胞检出，5%门槛用于表达覆盖敏感性分析。',
29:'接收者背景取受检基因、至少两个合格供者中达到表达门槛的基因与配体–靶先验的交集。背景≥500、靶集≥10、非靶背景≥10且候选配体≥2时计算配体活性；名义与FDR靶集分别保留标签。配体–受体表达及复合物兼容性、配体活性、预测靶与疾病差异基因的重叠、同供者关联作为四类互补读数分别报告。[10]',
30:'人类使用nichenetr 2.2.1.1的人类经典先验（Zenodo 3260758），小鼠使用v2小鼠先验（Zenodo 7074291）。基因映射记录物种、别名与重复项；跨物种按同源基因及方向比较。校正AUPR衡量配体预测对接收端疾病靶集的区分能力并用于排序；先验为无符号网络，配体的作用方向结合独立证据解释。',
31:'同供者关联要求两组各至少两人、两端细胞覆盖合格且残差自由度≥3，调整疾病、来源和两端秩深度。BH在接收范围、背景筛选、差异模型、靶集门槛及方向内，跨配体–发送端组合计算。该关联检验使用可评分的疾病靶集，其可估计性独立于NicheNet配体活性的靶集大小门槛；因此两类检验分别报告。',
32:'配体活性以各巨噬来源通过单对LR门槛的配体并集为候选，在固定接收端疾病靶集和背景上计算一次，再将同一AUPR关联到各来源的表达记录。IL15接收端复合物检查要求IL2RB和IL2RG各自在至少两名合格病例中达到10%检出，两条链允许由不同病例满足。另行统计严格四基因同供者覆盖：两端各≥20细胞，发送端IL15与IL15RA、接收端IL2RB与IL2RG的边际检出率均≥10%。单对LR、接收端链门槛、四基因同供者覆盖及同细胞共检出分别记录。',
36:'去候选基因定义移除IL1R2、PENK、TTC39C、FOXP1、TNFRSF9、LAG3、MRE11七个指定候选及TF自身；去程序重叠定义进一步移除12程序的基因并集，在相同表达排序及种子下重新评分。各定义记录每个种子的靶基因数。调控子–程序偏相关对读数取秩，残差化疾病、来源及接收者深度。Core组成敏感性在此基础上加入内部亚群比例，并省略一个比例以避免共线性；完整、去候选、去重叠定义及背景筛选分别报告。',
37:'独立时间点分析使用主小鼠队列的固定调控子成员，映射至GSE183682可测基因后评分。每条件每时间点1个文库，结果为描述性组差，P与FDR标记为未计算（无组内生物学重复）。',
38:'IL7R背景扣除疾病程序（excess）取来源校正、名义P<0.05且log₂FC≥0.5的88个上调靶基因，计算平均log1p(CP10k)，再减去三组固定匹配参照评分的等权均值；每组参照含88个基因。在固定16名供者（9病例、7对照）中，将excess与IKZF2去重叠AUC进行偏秩相关，调整疾病、来源和接收者秩深度。该项为同一数据内选择疾病程序后的探索性特异性诊断，报告未经多重校正的P，FDR记为未计算。三种子与原始程序、参照、excess及加入参照协变量的24条诊断均予保留。10%/20%调控子成员删除沿用同一excess终点，六项预定高权重删除关联构成探索性BH家族；匹配删除对照报告描述性结果。',
40:'供者重采样分别考察MacAM IL18与IL7R亚群IKZF2、MacAM TNF与核心CD4 STAT1两组关系。效应程序由去共享基因后的PRF1、IFNG、NKG7组成，分别计算配体–TF、TF–程序和配体–程序相关。本节三基因效应终点与方法1.6、结果2.6的88基因背景扣除疾病程序具有不同成员及计分定义。',
41:'主重采样按“来源×疾病”分层，整行有放回抽取供者，每个方案5000次。固定17人名单内的RNA及三个AUC种子共享抽样索引，用于配对比较；19人与17人的方案分别实施。基础关联模型调整疾病、来源及两端秩深度，TNF关系另加入Treg比例作敏感性分析。报告有效抽样的百分位95%区间，采用线性分位数插值。[12]',
42:'秩不足、残差自由度<3、恒定或接近零残差的抽样记为不可估计，并保留失败次数；有效抽样少于1000时不报告区间。未分层重采样用于评价来源分层对可估计性的影响。三个相关方向同时符合的比例定义为重采样代数相容率，其解释范围见讨论3.3。',
44:'按供者及细胞范围汇总原始计数。人类GSEA排序来自来源校正的edgeR准似然模型，使用filterByExpr、过滤后重置文库大小、TMM及稳健离散度估计。全体与核心CD4排序分别含10348和10961个受检基因，排序值为sign(logFC)×√QLF。通路分析使用过滤后文库大小；通讯疾病差异保留完整原始文库大小，两者分别构建偏置项。',
45:'人类基因集取自msigdbr的GO Biological Process和KEGG Legacy，GSEA使用clusterProfiler的fgsea实现，基因集大小10–500，保留全部P值，eps=0，随机种子20260813。[13] 软件版本为msigdbr 26.1.0、clusterProfiler 4.18.4，实际成员保存在分析对象中。BH按细胞范围和数据库分别计算。GSEA用于评价疾病排序中的基因集富集。',
49:'功能程序分别采用表达匹配评分与UCell秩评分。表达匹配评分由签名基因平均log1p(CP10k)减去匹配对照基因均值得到：按供者等权基线表达划分25个分位箱，每个签名基因最多选择20个表达接近的对照。控制池排除候选、全部签名、背景、线粒体编码及核糖体基因。至少3个可测且表达的签名基因、60%成员覆盖及20个控制基因时生成评分。',
50:'12个固定程序包括核编码氧化磷酸化、线粒体整合应激、氧化还原反应、线粒体自噬、炎症小体、焦亡执行相关、凋亡执行相关、T细胞活化、调节性、抑制性、初始状态及增殖。凋亡执行签名为BAX、BAK1、APAF1、CASP9、CASP3、CASP7、CYCS；焦亡执行签名为PYCARD、CASP1、GSDMD、IL1B、IL18。两类签名分别计分，完整成员见表S2。',
51:'人类UCell在单细胞表达排序上计分，再汇总供者均值。[15] 整体分析每范围每供者≥20细胞，亚群分析≥10细胞。9个面板包含线粒体应激、氧化还原、自噬、炎症小体感知、经典及非经典焦亡、Gasdermin替代通路等，成员见表S3。UCell氧化还原签名含FTL，共13基因；表达匹配签名为不含FTL的12基因。UCell的GSDMAlt按其五基因成员定义，与表达匹配的独立凋亡执行签名分别解释。',
52:'小鼠UCell使用原始计数、maxRank=1500、UCell 2.14.0及种子42，对七个GO模块及手工焦亡执行面板（Casp1、Casp4、Gsdmd、Gsdme）计分，并从成员中移除Mre11/Mre11a。每文库平均评分的组差使用全部20种3+3标签分配的双侧精确置换；主要三个细胞范围的模块共同校正，低细胞Treg另设探索性家族。MRE11相关结果作为独立参考分析列于表S14。',
54:'UMAP使用Seurat对象中保存的二维坐标。全肺对象为GSE288887整合集合，含其复用的GSE135893及GSE227136数据，细胞类型按广谱标记注释lineage_prelim。全肺与免疫子集共享坐标；CD4亚群及来源面板使用CD4对象的独立坐标。细胞标识、供者、疾病组、来源与亚群均与分析清单对应。UMAP展示细胞组成及来源分布，疾病比较采用供者模型。',
55:'GSEA曲线读取gseaResult对象，使用来源校正伪批量排序sign(logFC)×√QLF、对应基因集及检验结果。在Core和Primary的GO-BP、KEGG完整检验集合中分别按|NES|降序选取5条，并列时按通路ID排序。该展示规则适用于全部通路名称。enrichplot::gseaplot2绘制加权running ES、各通路命中条带和完整排序统计量，权重指数为1；rel_heights=(1.5,0.4,1)、subplots=1:3、pvalue_table=FALSE、ES_geom=\'line\'、base_size=11。各曲线极值与enrichmentScore、命中数与setSize逐项核对。独立图为10×8英寸PDF及300 dpi PNG，图下列NES、P及BH；BH覆盖相应范围×数据库的全部检验条目，展示的5条沿用完整家族校正值。',
59:old[59]['text'].replace('真实UMAP','UMAP'),
64:old[64]['text'].replace('完整和背景筛选版本','完整细胞与背景筛选结果').replace('由此，人类患者的CD4状态变化具有亚群特异性，调节性和抑制相关转录读数分别构成后续通讯分析的接收端背景。','这些读数呈现亚群特异的转录表型。通讯分析进一步按亚群定义伪批量疾病靶集及表达背景。'),
69:old[69]['text'].replace('说明不同接收亚群对应不同的巨噬信号组合。 这些AUPR是接收端共用预测值，逐来源表达与关联证据见表2及讨论3.2。','显示配体优先级随接收亚群改变。这些AUPR为接收端共用预测值，各巨噬来源的表达及同供者关联见表2，证据分层见讨论3.2。'),
76:'2.5　IKZF2与STAT1的独立调控网络分析',
82:old[82]['text'].replace('该组为同数据诊断中的名义关联，原分析未计算FDR。','该组为同数据特异性诊断中的名义关联，FDR记为未计算。'),
87:old[87]['text'].replace('该历史组间模型未计算FDR，以上为名义差异。','该AUC组间参照模型的FDR记为未计算，以上按名义差异报告（输出条件见表S1）。'),
98:old[98]['text'].replace('GSEA曲线见图7、图8，','全基因排序中按|NES|选取的Top5曲线见图7、图8；上述凋亡条目的数值见表S4，').replace('及既有剔除敏感性结果','及基因剔除敏感性结果'),
103:old[103]['text'].replace('这些名义结果保留了计分方法、细胞门槛和供者数的具体定义。','两种计分方法的名义效应方向相同，实际供者覆盖及BH结果分别报告。'),
111:old[111]['text'].replace('供者平均评分','文库平均评分'),
114:'患者资料显示亚群特异的组成和转录状态。GZMA调节性相关程序与CCR7抑制相关程序在来源及深度模型中主要获得名义支持，提示相应转录状态的变化；完整Treg抑制功能仍需功能读数评估。核心/IL7R五个固定面板的30项疾病比较均未达名义P<0.05，其配体相关结果因此解释为程序共同变化，而非已确立的疾病表型。',
115:'人类病例诊断为煤工尘肺（CWP），疾病关联结论首先适用于该患者人群；其对矽肺研究的意义结合小鼠硅尘模型作方向比较。来源与疾病分布不平衡：GSE135893为2名对照，GSE288887为7名病例，GSE227136为8名对照与2名病例。GSE227136提供唯一的来源内病例–对照交叉，因此来源校正后的疾病效应主要由其中两名病例支持。全体及核心CD4组成在基础模型中呈名义差异，加入来源后未维持名义显著。该设计限制了结果向其他患者来源和纯矽肺人群的推广。',
116:'亚群覆盖决定各分析的有效样本量：至少20细胞时，CCR7为6对照/3病例，GZMA为7/5，IL7R为7/9，Treg为5/5；CD4state只有1名对照达到门槛，其内部疾病模型及通讯DE缺少足够对照支持。背景筛选将可评分AUC名单由17人变为15人。RNA、AUC、比例和程序分别使用各自合格名单；Bootstrap中的RNA17专指与AUC配对的RNA名单。',
119:'IL15的来源证据存在明显差异。MacAM仅一名合格病例达到IL15的10%检出门槛，低于至少两名病例的条件，因此该来源同供者关联标为未计算。MonoMac与MacWide均有两名病例达标。IL7R的单对LR由IL2RG支持，IL2RB达标病例为零；Core的IL2RB/IL2RG边际链门槛通过。严格四基因同供者10%条件在表2六种组合中均无合格病例，IL15RA转呈仍缺少完整链条支持。这些结果描述表达覆盖与复合物兼容程度。',
124:'MacAM的簇2/11分别有98.8%/97.2%的细胞来自单一对照供者，MonoMac具有单核/巨噬混合特征，MacWide则共享前两者的供者。来源归属需结合单供者主导和注释异质性判断，各集合的结果具有相关性。低检出受体及复合物成员记为覆盖不足；5%敏感性用于评估检测门槛影响，主结论仍按10%条件解释。',
128:'GZMA来源覆盖进一步定位了这一表达关联（表3）。达到≥20细胞的模型名单为12人：G227含2名病例与7名对照，G288含3名病例；G135无人达到门槛。该范围的来源内病例–对照信息来自G227，IKZF2检出也集中于该来源。G288缺少来源内对照，其零检出记录主要反映测量覆盖。效应解释以供者模型为依据，细胞数用于描述检出分布。',
132:old[132]['text'].replace('历史伪批量模型中','伪批量疾病模型中').replace('；不能把供者均值模型的校正结果套用于该模型。','；该检验家族与供者均值模型分别定义。').replace('较新的六个细胞范围、三个供者均值模型','六个细胞范围、三个供者均值模型'),
135:old[135]['text'].replace('组成敏感性使用项目已保存的模型：','Core组成敏感性模型').replace('这些结果支持保留程序关联及其组成敏感性，但不能把去基因重叠或一次组成调整解释为排除全部细胞状态贡献，更不能称为IKZF2抗凋亡、抗焦亡或因果保护作用。','这些结果表明程序关联随组成调整而减弱，仍可能包含细胞状态构成的贡献，适合解释为观察性转录关联。'),
136:'表4　Core中IKZF2与程序的组成调整关联',
143:'供者重采样支持部分TF–程序关联，同时显示配体–TF关系对模型设置敏感。方向相容率描述相关符号的组合；现有证据尚不足以确定TNF/IL18经任一候选引起终末表型的机制。IKZF2与STAT1依据各自的表达、亚群和程序结果独立解释。',
149:'GSE314943的3+3代表六份命名scRNA文库，疾病比较以文库为分析单位。公开元数据尚缺少完整的动物个体或混样至文库映射，因此生物学重复的独立程度仍有不确定性。人鼠对照主要用于比较亚群与程序方向。',
151:'凋亡执行的BAX/BAK1–APAF1–CASP9/3/7与炎症小体–CASP1–GSDMD分别代表不同分子程序，替代Gasdermin面板又使用另一组成员。GO调控条目的富集及候选与转录评分的相关，反映基因集层面的共同变化。IL1B/IL18成熟分泌、CASP/GSDM裂解、膜破裂和实际细胞死亡需要相应功能测量；IKZF2的保护作用也有待此类证据支持。',
153:old[153]['text'].replace('不证明恢复轨迹或疾病进展机制','用于不同时间点的方向对照，尚不足以确定纵向疾病轨迹'),
154:old[154]['text'].replace('支持避免把其作为Treg专属身份标记','说明其表达也涉及常规T细胞状态').replace('Szabo等近期人类多组织研究','Szabo等人类多组织研究').replace('不能直接视作成年尘肺肺CD4的验证','与成年尘肺患者肺CD4在年龄和组织背景上存在差异'),
155:old[155]['text'].replace('这一组合与IKZF2不同，不应据同一配体模型将两者捆绑为共同机制。','这一组合为STAT1提供独立于IKZF2的研究依据。').replace('不能替代CD4特异STAT1操作','主要说明髓系及整体干扰素信号的作用范围，CD4特异STAT1效应仍需另行验证'),
156:old[156]['text'].replace('提醒总体CD4或调节性标签不足以决定功能方向；该研究并未验证本数据的人类IL15候选','为按亚群解释功能提供依据；本研究的人类IL15结果则来自配体–靶预测').replace('也不能作为任一候选的特异性因果实验','其作用范围与候选特异性干预有所不同'),
158:'3.6　资源覆盖与后续研究问题',
159:'软件记录包含msigdbr 26.1.0与clusterProfiler 4.18.4，实际基因集成员已保存；人类MSigDB的确切发布号尚缺少完整记录。人类UCell采用默认接口，调用中未显式记录maxRank，复现应以保存对象、成员和软件环境为依据。MRE11作为独立参考分析列于表S14。',
160:'临床严重度、部分细胞来源及低检出信号的解释仍受数据覆盖限制。本研究提供转录及计算通讯证据，候选敲除、配体阻断、蛋白裂解和代谢通量属于后续功能验证问题。研究设计可按接收亚群、巨噬来源、独立候选与功能终点分别展开。完整检验结果及靶集、背景成员作为可复现资料保存。',
162:'本研究在煤工尘肺患者中呈现亚群特异的CD4组成和转录状态，并获得与接收者疾病程序相对应的配体优先级。IKZF2与STAT1分别关联不同的表达、调控子和程序读数；小鼠硅尘模型补充了相应亚群及程序的方向对照。呼吸链富集、IL7R线粒体自噬相关评分及IKZF2与死亡相关转录程序的关联提出了可进一步检验的功能问题，其支持程度和敏感性按具体终点解释。',
164:'补充表汇总各分析的效应、检验结果及样本覆盖，保留正负方向。统计缺失分为“未计算”和“不可估计”，原因由模型输出条件、表达覆盖或设计可估计性说明。完整候选及配体–靶、靶集、背景成员另存可复现资料。',
183:old[183]['text'].replace('×版本×','×调控子定义×').replace('三个调控子版本','三种调控子定义').replace('八个既有TF','八个指定TF').replace('去候选版本移除既有七基因','去候选定义移除指定七基因').replace('背景筛选同时改变供者数，不能把差异只归因于评分。','背景筛选同时改变供者覆盖及评分输入。'),
187:'所有行n=17、残差df=11；完整细胞、去候选及程序重叠调控子，调整疾病、来源和接收者秩深度。小鼠相同关联模型n=6，报告相关方向；P/BH记为未计算，原因是该关联模型的推断输出要求总分析单位数至少10（表S1）。',
209:old[209]['text'].replace('未生成新P/BH，不设置空FDR列。','Bootstrap报告相关系数及区间，P/BH未计算。').replace('原始重采样表','重采样明细表'),
210:'4.10　小鼠评分、调控子转移及MRE11参考分析',
215:'表S14　人类MRE11参考分析，独立于候选TF分析',
217:'调整模型为面板评分~来源+疾病+MRE11检出比例+标准化log10(中位UMI+1)。该表专门报告MRE11读数。小鼠MRE11的128项聚合相关记录另存参考资料，与IKZF2和STAT1候选分析分别解释。',
219:'分析资源包括人类与小鼠基序数据库、NicheNet人类经典和小鼠v2先验、MSigDB来源的GO/KEGG成员，以及小鼠GO/KEGG映射。源文件、哈希、物种、成员和映射记录与各分析输入对应，用于核对资源及评分定义。',
220:'可复现资料包括504项程序模型、人类通讯全部候选及靶集/背景/配体–靶成员、小鼠结果表、12个网络的完整成员、供者重采样诊断和GSEA/GSVA敏感性结果。代码包提供运行顺序和依赖说明；个体数据、大型对象及模型数据库由服务器单独保存。',
}
new[6]=old[6]['text'].replace('小鼠核心CD4的两候选调控子活性均呈上调。','小鼠核心CD4的两候选调控子活性均呈名义上调。')
# Precise terminology changes without changing the plotted figures or numeric results.
new[61]=old[61]['text'].replace('真实UMAP','UMAP').replace('全肺冻结UMAP','全肺UMAP').replace('原广谱标记初步注释','广谱标记初步注释').replace('；My为广谱髓系标签，未将其全部改称巨噬细胞。','；My为广谱髓系标签。').replace('独立冻结CD4嵌入','独立CD4嵌入').replace('GSE288887新样本','GSE288887自身收录样本').replace('A/B来源为原全肺scored对象','A/B来源为全肺scored对象').replace('，未重新整合','')
for i in [99,100,225,228]:
    new[i]=old[i]['text'].replace('原来源校正分析','来源校正分析').replace('个原检验条目','个检验条目').replace('真实加权running ES','加权running ES').replace('图下列原NES','图下列NES')
new[99]=new[99].replace('KEGG疾病名称为基因集名称，不用于判断患者伴随疾病。','KEGG疾病名称在本图中表示相应基因集。')
new[225]=new[225].replace('KEGG疾病名称为基因集名称，不用于判断患者伴随疾病。','KEGG疾病名称在本图中表示相应基因集。')
for i in [100,228]:
    new[i]=new[i].replace('通路名称描述基因集，本图不单独判定该生物过程已发生。','通路名称用于标识基因集，曲线描述其在疾病排序中的分布。')
new[181]='表S5　IKZF2与凋亡/焦亡评分：各调控子定义'
new[138]=old[138]['text'].replace('×调控子版本×','×调控子定义×')
new[172]=old[172]['text'].replace('不称为官方GO通路','与数据库GO通路分别定义')
new[104]=old[104]['text'].replace('Inflam为炎症小体；','Inflam在人类为炎症小体签名，在小鼠为GO炎症小体介导信号通路；').replace('MitoOrg/OxStress/InflamRsp/DDR为','PyroResp为焦亡炎症反应；MitoOrg/OxStress/InflamRsp/DDR为')
new[222]=old[222]['text'].replace('ApoExec为凋亡执行，','ApoExec在此为GO凋亡执行条目，')

changes=[]
for i,text in new.items():
    e=nodes[i]
    if old[i]['type']=='figure': e=e.xpath('./p[@class="Caption"]')[0]
    assert e.tag!='table'
    before=e.text_content()
    for child in list(e):e.remove(child)
    e.text=text
    changes.append(dict(block=i,type=old[i]['type'],before=before,after=text))

# Table edits are limited to methods/descriptive labels; all numerical result cells remain intact.
table_changes={
19:{(0,1):'分析单位及细胞覆盖',(0,2):'分析范围与解释',(1,2):'描述肺细胞构成；CD4分析另取19人名单'},
167:{(2,0):'供者均值RNA/AUC疾病模型',
     (2,2):'读数×调控子定义×种子×背景×模型内，跨范围×两个TF',
     (3,0):'AUC组间参照模型',
     (3,2):'FDR未计算：本项输出条件为每组≥10个分析单位；报告已计算P',
     (4,2):'物种×背景×种子×调控子定义×范围×模型内全部受检TF×程序/候选表达特征；总分析单位n≥10才输出P/BH',
     (9,2):'报告相关及区间，P/BH未计算；失败不补抽，配对差共享索引'}
}
for i,replacements in table_changes.items():
    rows=nodes[i].xpath('.//tr')
    for (ri,ci),text in replacements.items():
        cell=rows[ri].xpath('./th|./td')[ci]
        before=cell.text_content()
        for child in list(cell):cell.remove(child)
        p=etree.SubElement(cell,'p');p.text=text
        changes.append(dict(block=i,cell=[ri,ci],before=before,after=text))

# Remove remaining workflow-history wording where it does not denote a real method operation.
replace_terms={'在组成比较之后，进一步考察':'按亚群考察','GSEA 的既有排序':'GSEA排序','完整12程序和去重叠版本':'完整12程序和去重叠定义'}
for i,e in enumerate(nodes):
    if e.tag=='table' or i>=232:continue
    target=e.xpath('./p[@class="Caption"]')[0] if old[i]['type']=='figure' else e
    text=target.text_content()
    changed=text
    for a,b in replace_terms.items():changed=changed.replace(a,b)
    if changed!=text:
        for child in list(target):target.remove(child)
        target.text=changed
        changes.append(dict(block=i,before=text,after=changed))

tree.xpath('//title')[0].text=new[0]+'——'+new[1]
out=D/'provenance/logic_render_source.html'
out.write_text('<!doctype html>\n'+etree.tostring(tree,encoding='unicode',method='html'),encoding='utf-8')

def image_hashes(t):
    return [hashlib.sha256(base64.b64decode(x.get('src').split(',',1)[1])).hexdigest() for x in t.xpath('//img')]
orig_tree=html.fromstring(source.read_text(encoding='utf-8'))
assert len(tree.xpath('//img'))==14 and image_hashes(tree)==image_hashes(orig_tree)
assert len(tree.xpath('//table'))==18
for i,e in enumerate(nodes):
    if e.tag=='table' and i not in table_changes:
        actual=[[c.text_content() for c in r.xpath('./th|./td')] for r in e.xpath('.//tr')]
        assert actual==old[i]['text'],i
decimal=re.compile(r'(?<![A-Za-z])[-−+]?\d+\.\d+(?:[eE][−+-]?\d+)?')
numeric_changes=[]
for c in changes:
    before=Counter(decimal.findall(c['before']));after=Counter(decimal.findall(c['after']))
    if before!=after:numeric_changes.append(dict(block=c['block'],cell=c.get('cell'),removed=dict(before-after),added=dict(after-before)))
# Decimal differences here can only be explicitly reviewed section references, not results.
(D/'QA/decimal_diff_review.json').write_text(json.dumps(numeric_changes,ensure_ascii=False,indent=2),encoding='utf-8')
(D/'provenance/editorial_changes.json').write_text(json.dumps(changes,ensure_ascii=False,indent=2),encoding='utf-8')
(D/'QA/image_and_table_invariance.json').write_text(json.dumps(dict(figures=14,figure_hashes_unchanged=True,tables=18,unchanged_full_tables=16,table_text_only_edits={str(k):{str(cell):v for cell,v in values.items()} for k,values in table_changes.items()},all_result_cells_unchanged=True),ensure_ascii=False,indent=2),encoding='utf-8')
print('EDITS',len(changes),'DECIMAL_DIFFERENCES',json.dumps(numeric_changes,ensure_ascii=True))
alltext=tree.text_content()
history=re.compile(r'新旧|旧分析|新增|较新|较早|旧限制|旧规则|旧模型|历史|补跑|本轮|本次|当前版本|曾经|恢复|原流程|原研究|既有|审阅修订版|阶段[一二三四五六七八九十0-9]')
hits=[]
for i,e in enumerate(nodes):
    if i>=232:continue
    text=e.text_content()
    found=history.findall(text)
    if found:hits.append(dict(block=i,hits=found,text=text))
(D/'QA/history_term_scan.json').write_text(json.dumps(hits,ensure_ascii=False,indent=2),encoding='utf-8')
print('HISTORY_HITS',json.dumps(hits,ensure_ascii=True))
