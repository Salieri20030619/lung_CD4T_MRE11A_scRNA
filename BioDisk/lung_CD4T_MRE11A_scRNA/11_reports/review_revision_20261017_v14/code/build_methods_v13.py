from pathlib import Path
import json,re,html,shutil,hashlib
from collections import Counter
import pandas as pd
from bs4 import BeautifulSoup,Tag
import fitz
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,PageBreak,Image,Table,TableStyle,KeepTogether
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors
R=Path('C:/codex_test');T=R/'tmp/silicosis_v13';O=R/'output/silicosis_v13';OLD=R/'output/silicosis_v12';SRC=T/'sources_R2';V4=R/'tmp/silicosis_v9/v4'
for folder in ['assets','figures','plot_data','source_tables','supplement']:shutil.copytree(OLD/folder,O/folder,dirs_exist_ok=True)
for folder in ['provenance','code','provenance/method_sources']:(O/folder).mkdir(parents=True,exist_ok=True)
hs=BeautifulSoup((OLD/'report_v12.html').read_text(encoding='utf-8'),'html.parser');cover=hs.find_all('section')[0];cover.find('p').string='完整方法学修订版 v13｜2026年10月16日'
oldsections=hs.find_all('section')[2:];design=oldsections[0];results=oldsections[3:]
assert results[0].find('h1').get_text().startswith('2 ')
# Only the cross-reference to the expanded scoring Methods changes in Results.
for s in results:
 for p in s.find_all('p'):
  if '评分单位和样本汇总方式与1.2节一致' in p.get_text():p.string=p.get_text().replace('评分单位和样本汇总方式与1.2节一致','评分单位和样本汇总方式与1.3节一致')
methods=[];trace=[];used=[]
def add(title,body,entries):
 k=len(methods)+1;ss=BeautifulSoup(f'<section id="method-{k}" data-chapter="no"><h2>{title}</h2>{body}</section>','html.parser').section;methods.append(ss)
 for topic,param,out,files in entries:
  rows=[]
  for origin,path in files:
   root={'archive':SRC,'v4':V4,'v12':OLD}[origin];f=root/path;assert f.is_file(),f
   dest=O/'provenance/method_sources'/origin/path;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(f,dest)
   used.append({'source':str(f),'package_path':dest.relative_to(O).as_posix(),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
   rows.append(dest.relative_to(O).as_posix())
  trace.append({'方法小节':title,'分析环节':topic,'核实参数与实施':param,'对应结果或补充表':out,'方法来源文件':';'.join(rows),'本轮重算':'否'})
add('1.1 数据预处理、质量控制与细胞注释','''
<h3>主发现队列GSE314943</h3>
<p>采用GEO提供的计数矩阵及已冻结的样本映射，保留3个Saline与3个Silica文库。逐样本保留检测基因数≥200、总UMI≥500、线粒体UMI比例≤15%的细胞。使用Scanpy接口逐样本运行Scrublet，预期双细胞率0.06、随机种子42、自动阈值；从主对象移除预测双细胞。跨样本按Ensembl标识取共同特征，排除符号冲突，原始整数计数保存在counts层。</p>
<p>全肺初步聚类后提取T/NK/ILC候选区室，共6303个细胞。采用CP10k归一化和log1p变换，按样本识别2500个高变基因；从嵌入特征中排除预设背景、线粒体、核糖体和血红蛋白基因。数据缩放截断于10，计算40个主成分。Harmony按sample校正PCA空间，随机种子42、最多30轮；随后使用前30个主成分、20近邻及Leiden分辨率0.6构建图谱。</p>
<p>依据T细胞、CD4/CD8、NK、γδT及Treg等成组标记共同确定身份。保留H0、H2、H4分别作为Naive、ActCD4、Treg，合计2311个Core CD4细胞；未按单细胞Cd4或候选基因阳性筛选。巨噬细胞保留AMlike、Spp1Mac、C1qMac、Folr2Mac等原冻结状态，AllMac为包含这些状态的汇总范围。</p>
<h3>分析对象与表达值分工</h3>
<p>Harmony仅用于嵌入和状态识别。差异分析使用原始计数伪批量；表达展示使用原始总文库UMI归一化值。补充小鼠队列复用各自冻结的CD4或TLike对象和真实取样时间标签。GSE174725重新进行广谱身份复核，GSE296682保留作者细胞亚型标签，具体流程见1.7节。</p>
''',[
('主队列QC/双细胞','genes≥200；UMI≥500；MT≤15%；Scrublet 0.06、seed42、自动阈值、移除预测双细胞','全部主队列结果',[('archive','upstream_scripts/GSE314943_stage2.py')]),
('T/NK重聚类与Harmony','2500 HVG；scale max10；PCA40；Harmony sample/seed42/max30；PC30、邻居20、Leiden0.6','Core及Naive/ActCD4/Treg范围',[('archive','upstream_scripts/GSE314943_stage3.py'),('archive','upstream_scripts/GSE314943_stage4.py'),('archive','upstream_scripts/GSE314943_stage5.py')])])
add('1.2 伪批量差异表达与候选基因选择','''
<h3>CD4疾病差异模型</h3>
<p>按独立小鼠将原始RNA计数求和，分别形成Core、Naive、ActCD4和Treg伪批量矩阵。设计矩阵为~group，Saline为参考组，正log2FC表示Silica上调。使用edgeR的filterByExpr过滤低表达基因，过滤后重新计算文库量，TMM归一化；estimateDisp与glmQLFit均采用robust=TRUE，glmQLFTest检验Silica相对Saline的组别系数。BH在每个范围、每个分析分支的全部受检基因内计算。</p>
<p>Full分支按完整基因输入；NoBG分支先移除预设背景与血红蛋白基因，再重新过滤与TMM归一化。Core、Naive、ActCD4均保留全部六个样本；Treg同样保留六个样本，包括13细胞的对照。正文复用冻结Full结果，NoBG作为特征/归一化敏感性分支。</p>
<h3>发送端差异与候选筛选</h3>
<p>巨噬细胞AllMac与AMlike差异复用对应的原始六文库edgeR准似然模型，设计同为组别对比。该发送端模型在基因过滤后保留原始总文库量，再计算TMM；原阶段每文库≥20细胞，两个范围实际均保留3对照+3病例。其BH保留原始完整受检家族。</p>
<p>历史初筛对每张差异表按原始P&lt;0.05及|log2FC|≥0.5，分别导出上调与下调前20个，再结合亚群效应、样本一致性、文献和背景敏感性确定优先级。后续新增候选以对应效应的原始P&lt;0.05作为名义筛选条件，不统一沿用历史倍数门槛。S17按同状态和变化方向比较高效应基因，并将原始筛选记录与本次文献评估分开记录。</p>
''',[
('CD4 edgeR','~group；filterByExpr；keep.lib.sizes=FALSE；TMM；robust dispersion/QLF；Full/NoBG分开BH','S01-S03及候选差异图',[('archive','02_analysis/GSE314943/stage6_20260905_231807_L44XdT/analysis.R')]),
('巨噬发送端DE','~grp；原始library offset；keep.lib.sizes=TRUE；TMM/robust QLF；AMlike与AllMac各6样本','S11发送端logFC/P/BH',[('archive','02_analysis/mechanism_step3_20260928_103057_lIRL38/analyze.R')]),
('历史及补充选靶','原历史P<0.05、|log2FC|≥0.5、上下各20；S17为后续同方向比较和文献补充','S17与3.3节',[('v12','provenance/baseline_v11/screening_history_excerpts.json'),('v12','provenance/baseline_v11/literature_sources.json')])])
add('1.3 功能程序评分、基因集与样本汇总','''
<p>各细胞表达先按原始总UMI归一化为CP10k，再取log1p。预设程序评分 = 程序成员平均log1p(CP10k) − 表达匹配对照基因均值。参考表达为各样本细胞均值的等权均值；按非零参考表达划分25个分位箱，对每个成员选同箱最多20个最接近的对照，同箱不足时从全部合格池选取，程序内去重。</p>
<table><tr><th>程序</th><th>本队列可测成员或基因集</th><th>匹配对照数</th></tr>
<tr><td>PyroExec</td><td>PYCARD、CASP1、GSDMD、IL1B、IL18（5/5）</td><td>100</td></tr>
<tr><td>ApoExec</td><td>BAX、BAK1、APAF1、CASP9、CASP3、CASP7、CYCS（7/7）</td><td>140</td></tr>
<tr><td>Inflam</td><td>NLRP3、AIM2、NLRC4、PYCARD、CASP1（5/5）</td><td>100</td></tr>
<tr><td>OXPHOS / MitoISR / MitoROS / Mitophagy</td><td>项目预设线粒体程序；全部成员与覆盖见S20</td><td>逐程序记录于S20</td></tr>
<tr><td>Autophagy_GO</td><td>org.Mm.eg.db的GOALL GO:0006914；577个注释符号中553个可测、535个有表达</td><td>7204</td></tr></table>
<p>预设程序的对照池排除历史候选IL1R2、PENK、TTC39C、FOXP1、TNFRSF9、LAG3、MRE11，全部预设程序成员，上皮/髓系背景标记，以及线粒体编码、核糖体和免疫球蛋白基因。可评估条件为至少3个可测且至少3个有表达成员、覆盖≥60%、匹配对照≥20个。All与BG1沿用相同成员和对照。</p>
<p>新增广义自噬评分单独按GO成员构建；其对照池另排除GO集合及IKZF2/STAT1，参考表达秩划分25箱，同箱按距离选最多20个对照并去重。每个样本及状态先汇总细胞评分均值，再对独立样本等权计算组均值。原始程序、GO集合、覆盖与对照清单分别保存于S18/S20。</p>
''',[
('预设功能模块','mean(log1p CP10k members) - mean(matched controls)；等样本参考；25箱；每成员最多20对照；≥60%覆盖等资格条件','图1-3、S04/S05/S18',[('v12','source_tables/presentation_inputs/analyze.py'),('v12','source_tables/presentation_inputs/module_definition.tsv'),('v12','source_tables/presentation_inputs/module_coverage.tsv'),('v12','source_tables/presentation_inputs/score_controls.tsv.gz')]),
('GO自噬','GOALL GO:0006914；参考秩25箱；553可测/535表达/7204对照','S05与S20_GO表',[('archive','11_reports/customer_sequence_20261015_v2/code/compute_v2.py'),('archive','11_reports/customer_sequence_20261015_v2/02_mito_autophagy/autophagy_coverage.tsv')])])
add('1.4 SCENIC调控网络、AUCell与TF靶基因先验','''
<h3>建网与调控子活性</h3>
<p>小鼠GSE314943的Core CD4用于SCENIC建网，All与背景负担≤1%的BG1为独立分支。基因须在至少max(10，1%保留细胞数向上取整)个细胞中检出，并同时位于两个cisTarget基因排名数据库中；不按HVG或疾病差异预选建网基因。输入为原始文库归一化的log1p(CP10k)。</p>
<p>每个样本最多无放回抽取250个细胞用于网络训练，种子为17、29、43；训练矩阵方差≤10⁻⁸的基因剔除。用GRNBoost2推断共表达边，再以mm10 RefSeq r80、cisTarget v10聚类基序进行剪枝，使用TSS上下游各10 kb及上游500 bp/下游100 bp两个窗口，配套MGI基序注释及小鼠TF列表。AUCell随后对全部保留细胞评分，抽样上限只用于建网。</p>
<h3>去重叠评分与正文分析</h3>
<p>DropTarget从每个调控子中去除历史七候选及TF自身；DropOverlap进一步去除冻结功能程序成员的并集。在相同表达/排名宇宙及相同种子下重新计算不加权AUCell；空调控子记录为不可评估，其余保留并记录大小。正文采用Mouse/All/DropOverlap、seed17，29和43用于算法敏感性。最终样本汇总重新纳入全部非空状态，包括低细胞数Treg样本，并对每个样本的细胞AUC求均值。</p>
<p>TF疾病效应以样本平均AUC计算；TF-程序关联用1.8节的样本级检验。STAT1与GSDMD等基因的关联，使用同样本的平均CP10k，并沿用预设目标面板中全部可评估基因的校正家族。GSDMD在所用去重叠STAT1调控子外。</p>
<h3>先验网络图</h3>
<p>图10的TF-靶基因边来自冻结的小鼠DoRothEA A/B/C级网络，保留调控符号与置信级别，再与项目程序成员相交；它与本队列SCENIC推断网络分别记录。图中实线的基因-模块连接表示集合归属。</p>
''',[
('SCENIC建网','pySCENIC0.12.1；GRNBoost2；cap250/样本；seed17/29/43；mm10 v10双窗口；全部细胞AUCell','TF活性图谱与S08',[('archive','02_analysis/scenic1_20260929_103513_xRoEnr/pipeline.py'),('archive','02_analysis/scenic1_20260929_103513_xRoEnr/environment_freeze.txt')]),
('DropOverlap/最终样本','候选+TF自身+冻结程序基因排除；noweights=True；非空调控子；最终所有非空样本重汇总','图7-9、S08/S18',[('archive','02_analysis/scenic2_20261001_092655_ZKNWD6/analyze.py'),('v4','code/prepare_mouse_v4.py')]),
('STAT1执行基因关联','预设21基因面板中可测对象；CP10k与样本AUC；按scope/seed/model完整家族校正','S09/S10',[('archive','11_reports/STAT1_pyro_bridge_20261015_v5/code/bridge_v5.py'),('archive','11_reports/STAT1_pyro_bridge_20261015_v5/code/sensitivity_v5.py')]),
('DoRothEA先验','冻结小鼠ABC网络；不以SCENIC边替代先验边','图10与S20',[('archive','02_analysis/mechanism_step6_20260928_163015_VoRQqA/network_source.json'),('v4','code/composition_v4.py')])])
add('1.5 拟时序、真实时间模型与辅助交集','''
<h3>细胞状态排序及样本级斜率</h3>
<p>DPT复用1734个Naive/ActCD4细胞及Harmony前30个主成分。分别构建15和30近邻图；选择同时含每个状态≥50个细胞、每个样本≥20个细胞的最大连通分量。根候选取Naive中Tcf7/Lef1/Ccr7平均logCP10k上四分位细胞，先取距其PCA中心最近者，再取另一样本中最近者，形成K15_R1、K15_R2、K30_R1、K30_R2四种设置。使用15维diffusion map、DPT前10个扩散分量，不作分支拆分。</p>
<p>后续推断分别在Naive和ActCD4内完成：重用上述DPT排序，按状态缩放到0—1；每个样本拟合表达或程序评分随DPT变化的斜率，调整log1p原始文库深度。六个样本斜率的均值通过全部64种符号翻转检验；基因在状态×轨迹设置的完整受检家族内BH校正，程序按每个设置的全部状态-程序组合校正。主设置为K15_R1，其他设置为敏感性。</p>
<h3>真实时间模型</h3>
<p>GSE274372的TLike伪批量按D0、D7、D28、D56分类时间建模，每组3个文库，全部12个非空文库纳入，包括8细胞D0文库。edgeR采用~0+day设计、filterByExpr、保留原始总文库量后的TMM、稳健离散度与准似然拟合。分别进行时间总体检验和D7/D28/D56相对D0对比，每项在全部受检基因内单独BH校正。GSE183682按每时间每条件1只小鼠输出表达方向。</p>
<h3>辅助候选交集</h3>
<p>交集依次要求：疾病差异P&lt;0.05；同状态DPT斜率P&lt;0.05且方向与疾病效应一致；真实时间总体P&lt;0.05；至少一个时间对比P&lt;0.05且同向。逐步数量、成员和未通过条件全部保留。该交集用于辅助支持，历史六候选与TF候选另按各自证据来源记录。</p>
''',[
('DPT生成','PC30；K15/30；Naive根上四分位；两个供者根；diffmap15/DPT10/no branching','3.2节辅助轨迹',[('archive','02_analysis/GSE314943/stage9_20260906_080538_wN3sjO/analysis.py')]),
('最终轨迹推断','按状态/样本调整log1p文库深度；64斜率符号翻转；完整家族BH','S06/S07',[('archive','11_reports/customer_sequence_20261015_v2/code/compute_v2.py')]),
('真实时间','12文库；~0+day；总体3对比矩阵及3个D0对比；独立BH','图6、S07',[('archive','11_reports/customer_sequence_20261015_v2/code/time_v2.R')]),
('交集规则','疾病P+轨迹P/同向+时间总体P+至少1同向时间对比P；均<0.05','S06/S07',[('archive','11_reports/customer_sequence_20261015_v2/code/intersect_v2.py')])])
add('1.6 配体受体覆盖、NicheNet与发送端关联','''
<h3>配体受体资格与复合物覆盖</h3>
<p>配体受体关系来自冻结的小鼠NicheNet v2网络。按样本和细胞状态汇总原始计数、平均CP10k和阳性比例。病例侧至少2个发送端文库检出配体、至少2个接收端文库检出兼容受体时，标为表达支持。10%细胞检出门槛保留为旧规则敏感性标记，当前主排名不以该比例作为准入条件。</p>
<p>AMlike与AllMac分别统计。受体复合物按同细胞各组分同时阳性计算，如IL2RB/IL2RG、IL18R1/IL18RAP；TNF的TNFRSF1A或TNFRSF1B采用任一阳性规则。配体复合物及IL15反式呈递相关组分单独记录。每个样本分别输出阳性细胞数、分母和检出比例。</p>
<h3>NicheNet活动预测</h3>
<p>采用ligand_target_matrix_nsga2r_final_mouse先验矩阵。各CD4范围的背景是edgeR受检基因、先验靶基因及表达输入的交集；上/下调目标分别按原始P&lt;0.05和log2FC正/负定义。至少10个目标、500个背景、10个非目标及2个潜在配体时运行predict_ligand_activities。输出AUROC、AUPR、corrected AUPR及Pearson预测指标；corrected AUPR为AUPR减去目标数/背景数。</p>
<p>活动值由接收端目标集合计算；再按发送端表达资格形成AMlike和AllMac来源排名。完整Core、Naive、ActCD4、Treg的上/下调预测保留于S19。图15从同一先验矩阵提取配体与死亡程序基因的权重，并标记ActCD4目标/背景资格及差异效应；其成员筛选与预测模型保持独立记录。</p>
<h3>同样本发送端—接收端联系</h3>
<p>在同一六只小鼠内，配对巨噬细胞配体平均CP10k与CD4亚群TF活性或程序评分，使用合并及疾病组别调整的秩关联模型。BH在相应来源、接收范围、种子和模型的完整配体-终点家族内计算；受体覆盖、来源差异、预测活动值和样本关联分别保存。</p>
''',[
('配体受体覆盖','≥2病例文库UMI>0；10%仅旧规则敏感性；复合物ALL、TNF受体ANY','图12-13与S11/S19',[('v4','code/prepare_mouse_v4.py')]),
('NicheNet当前预测','原始P<0.05上下分开；目标≥10/背景≥500/非目标≥10/配体≥2；corrected AUPR','图13-15、S19',[('v4','code/ligand_activity_v4.R'),('v4','code/export_prior_v4.R'),('v4','02_sender/methods_and_session.txt')]),
('先验权重与关联','未改变原矩阵权重；同鼠配体CP10k配对TF/程序；完整家族BH','S11/S19与图15-16',[('v12','code/export_inputs.R'),('v4','code/mouse_stats_v4.py')])])
add('1.7 人类队列复核、患者定位与组成分解','''
<h3>GSE174725独立供者分析</h3>
<p>保留2名暴露对照和3名矽肺病例，使用CSV恢复的无歧义共同基因空间及原始计数。原全CSV质量控制为检测特征≥200、计数≥300、线粒体比例&lt;20%。逐供者Scrublet采用预期率0.06和种子42；质量标记保留，主分析与排除既有双细胞/混合谱系标记的敏感性分别进行。</p>
<p>广谱注释使用CP10k-log1p、2500个HVG，排除预定候选、IKZF2/STAT1及线粒体/核糖体特征；scale截断10，PCA30、20近邻、Leiden0.6、种子42。结合T、髓系、巨噬、B和上皮等正/负标记、簇表达和供者分布冻结身份。CD4沿用已复核的MemCD4、NaiveCD4、Treg范围；AM_like与FCN1/FCGR3A单核细胞分开记录。</p>
<p>表达按每个供者内细胞CP10k均值计算，组间对供者等权。至少每组2个可评估供者时，对绝对组均值差进行全部标签置换；正文STAT1、GSDMD等采用v6预定31基因家族，在每个范围与敏感性版本内BH校正。留一供者只重算描述性差值。该队列的人类证据采用表达、组成和细胞定位，不引入煤工尘肺队列的Human SCENIC活性。</p>
<h3>Core组成分解</h3>
<p>对每组每状态，p为供者内该状态比例的组均值；μ为该状态对Core均值的平均贡献除以p，保留等供者Core定义。总差值分为Σ[(p病例−p对照)×(μ病例+μ对照)/2]与Σ[(μ病例−μ对照)×(p病例+p对照)/2]两部分，并验证二者之和等于Core病例-对照差值。空状态贡献按零处理，但不补造其表达值。</p>
<h3>GSE296682患者细胞定位</h3>
<p>从作者BAL_FINAL处理矩阵按元数据提取B2_RZ725、B6_RZ725两位矽肺患者，保留作者Cell_Subtype标签。以存储表达值&gt;0定义检出，在各患者、各标签内报告阳性数/总细胞数；复合受体及多基因共检出要求同一细胞全部组分阳性。以TCR/CD3、CD4/CD8、NK及成熟巨噬标记核对细胞背景，报告作者标签及身份核查结果。</p>
''',[
('GSE174725注释','QC200特征/300计数/MT<20；Scrublet0.06/seed42；HVG2500/PCA30/K20/Leiden0.6；保留旧CD4标签','人类病例对照与S12/S13',[('archive','11_reports/GSE174725_annotation_validation_20261005_v1/code/g174_build_annotation.py'),('archive','11_reports/GSE174725_annotation_validation_20261005_v1/code/g174_freeze_review.py')]),
('人类31基因/分解','供者CP10k均值；组别精确置换；31基因各scope/variantBH；对称贡献分解等式','S12/S13与图18',[('archive','11_reports/upstream_human_support_20261015_v6/code/support_v6.py')]),
('GSE296682定位','作者处理X>0；作者亚型；同细胞共检出；按每患者分母','S14-S16与图17',[('archive','11_reports/upstream_human_support_20261015_v6/code/patient2_v6.py'),('archive','11_reports/GSE296682_patient_localization_20261005_v1/code/g296_localize.py')])])
add('1.8 统计检验、敏感性分析与绘图规则','''
<h3>统一统计口径</h3>
<p>统计单位为独立小鼠或人类供者。预设目标效应的原始P&lt;0.05为名义显著，BH&lt;0.05为相应家族的校正后显著。BH仅对同类、完整且可评估的假设家族计算；DE、程序组间比较、TF-程序、STAT1-单基因、配体-接收端及外部队列分别记录。常量、缺失或不可评估项保留NA。</p>
<p>小鼠3对3的程序及AUC组间比较，以绝对组均值差枚举20种标签分配。死亡程序在Core三项和三个亚群九项中分别BH；线粒体/自噬在对应Core和亚群家族分别BH。合并样本关联对六样本秩值中心化，枚举720种配对；组别调整先在各组内中心化秩值，再枚举36种组内配对。TF-程序按范围×种子×模型的全部调控子-程序组合校正。留一法逐次删除一个独立样本，记录可计算次数和方向一致比例。</p>
<h3>组成与背景敏感性</h3>
<p>固定组成分析先计算每只小鼠Naive、ActCD4、Treg占Core的比例，再对六只小鼠取平均，得到统一权重46.66%、31.00%、22.34%；逐样本对亚群内AUC、表达或程序评分加权，之后重新计算相关系数及相应完整家族BH。图9使用相同样本编号和坐标范围对照实际组成与统一组成。</p>
<p>BG1定义为预设背景面板UMI/原始总UMI≤1%，与Full/NoBG的基因级分支分别记录。历史20%、10%、5%、2%、1%背景负担阈值、扣除背景UMI分母、逐样本剔除及固定组成用于表达稳健性核查。当前主样本汇总使用全部非空状态，旧阶段20细胞汇总或10细胞时间模型不作为当前结果的纳入规则。</p>
<h3>图件与复现记录</h3>
<p>TF活性热图取方差前22个调控子并加入STAT1/IKZF2，按行标准化、平均连接层次聚类；TF-程序图取三亚群名义关联数前18个TF并加入二者，仍展示其全部关联。配体图按已注明的来源内排名选取，显著性不在展示子集重新校正。记录原始P、BH、检验家族、样本数和缺失原因。</p>
<p>冻结环境包括Scanpy 1.11.5、AnnData 0.12.16、R 4.5.2、edgeR 4.8.2、pySCENIC 0.12.1、Arboreto 0.1.6、ctxcore 0.2.0及nichenetr 2.2.1.1；不同分析使用各自冻结环境。S20提供本节方法、参数、结果位置和来源文件的对照及SHA256。</p>
''',[
('最终统计家族与置换','20组别分配；720秩配对；36组内配对；各scope/source/seed/model完整家族BH','全部P/BH及S08-S11',[('v4','code/mouse_stats_v4.py'),('archive','11_reports/customer_sequence_20261015_v2/code/compute_v2.py')]),
('固定组成','平均样本比例；逐样本加权TF/程序/表达；相应完整家族重新计算','图9/S10/S18',[('v4','code/composition_v4.py'),('archive','11_reports/STAT1_pyro_bridge_20261015_v5/code/sensitivity_v5.py')]),
('版本与绘图','冻结版本；图件选择只作展示；保留原检验家族','S20',[('archive','02_analysis/GSE314943/stage2_20260905_190530_BeH7tA/tables/versions.json'),('archive','02_analysis/GSE314943/stage6_20260905_231807_L44XdT/tables/sessionInfo.txt'),('v12','code/draw_v12.py')])])

sections=[design]+methods+results
conclusion=next(s for s in results if s.find(['h1','h2']).get_text().startswith('9 '));tr=hs.new_tag('tr')
for text in ['S20','完整方法、参数、对应结果及来源记录；线粒体/自噬成员与匹配对照']:
 td=hs.new_tag('td');td.string=text;tr.append(td)
conclusion.find('table').append(tr)
nav=BeautifulSoup('<section><h1>阅读路径</h1><p>数据与完整方法 → 功能表型 → 候选 → TF调控 → NicheNet与通讯 → 人类证据 → 实验选择</p></section>','html.parser').section
for ss in sections:
 if ss.get('data-chapter')=='yes':
  p=hs.new_tag('p');a=hs.new_tag('a',href='#'+ss['id']);a.string=ss.find(['h1','h2']).get_text();p.append(a);nav.append(p)
p=hs.new_tag('p');p.string='方法学按1.1—1.8节覆盖预处理、差异分析、程序评分、SCENIC、轨迹/时间、通讯、人类分析及统计敏感性。图1—18和既有结果保持原值；S20提供可追溯参数表。';nav.append(p)
css=hs.find('style').get_text();(O/'report_v13.html').write_text('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>矽肺完整方法学报告 v13</title><style>'+css+'</style>'+str(cover)+str(nav)+''.join(map(str,sections))+'</html>',encoding='utf-8')
fig_no={Path(img['src']).name:int(re.match(r'图(\d+)',img.find_next_sibling('p',class_='caption').get_text()).group(1)) for ss in results for img in ss.find_all('img')}
engine=(R/'tmp/silicosis_v10/build_integrated.py').read_text(encoding='utf-8').split('# Generate a continuous typeset report')[1].split('# Content conservation checks')[0].split('\n',1)[1].replace('v10','v13')
engine=engine.replace('{1:300,2:300,3:240,4:300,5:140,8:200,10:160,14:250,15:250}','{1:230,2:230,3:300,4:240,5:300,6:140,9:295,11:160,14:260,15:290,17:250,18:250}')
exec(engine,globals())
pd.DataFrame(trace).to_csv(O/'supplement/S20_methods_parameters_results_sources.tsv',sep='\t',index=False,encoding='utf-8')
pd.DataFrame(used).drop_duplicates('package_path').to_csv(O/'provenance/S20_method_source_manifest.tsv',sep='\t',index=False,encoding='utf-8')
for src,dest in [('module_definition.tsv','S20_all_fixed_module_genes.tsv'),('module_coverage.tsv','S20_all_fixed_module_coverage.tsv'),('score_controls.tsv.gz','S20_all_fixed_module_controls.tsv.gz')]:shutil.copy2(OLD/'source_tables/presentation_inputs'/src,O/'supplement'/dest)
for stem in ['autophagy_measured_members','autophagy_control_members','autophagy_coverage']:shutil.copy2(SRC/f'11_reports/customer_sequence_20261015_v2/02_mito_autophagy/{stem}.tsv',O/f'supplement/S20_GO_{stem}.tsv')
(O/'provenance/page_map.json').write_text(json.dumps(bookmarks,ensure_ascii=False,indent=2),encoding='utf-8')
d=fitz.open(O/'silicosis_integrated_report_v13.pdf')
for i,page in enumerate(d):
 assert '\ufffd' not in page.get_text();assert all(b[0]>=-1 and b[1]>=0 and b[2]<=843 and b[3]<=596 for b in page.get_text('blocks')),(i,page.get_text())
 page.get_pixmap(matrix=fitz.Matrix(1,1)).save(T/f'page_{i+1:02d}.png')
shutil.copy2(__file__,O/'code/build_methods_v13.py');print('REPORT COMPLETE',len(d));print(json.dumps(bookmarks,ensure_ascii=False))
