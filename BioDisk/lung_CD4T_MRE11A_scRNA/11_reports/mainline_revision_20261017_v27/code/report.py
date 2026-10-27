from pathlib import Path
import json,html,hashlib
import pandas as pd
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.platypus import BaseDocTemplate,PageTemplate,Frame,Paragraph,Spacer,PageBreak,NextPageTemplate,Table,TableStyle,Image,KeepTogether
from PIL import Image as PILImage
O=Path(__file__).parent;H=O.parent;D=H/'Downloads/矽肺_CD4候选研究_主线修订_v27_20261017';D.mkdir(exist_ok=True)
manifest=json.loads((O/'figure_manifest.json').read_text(encoding='utf-8'));assert len(manifest)==24
pdfmetrics.registerFont(TTFont('CN','C:\\Windows\\Fonts\\ARIALUNI.TTF'))
styles={'p':ParagraphStyle('p',fontName='CN',fontSize=10.5,leading=17,spaceAfter=9,wordWrap='CJK'),'h':ParagraphStyle('h',fontName='CN',fontSize=20,leading=28,spaceAfter=16,textColor=colors.HexColor('#143e55'),wordWrap='CJK'),'h2':ParagraphStyle('h2',fontName='CN',fontSize=13,leading=21,spaceAfter=9,textColor=colors.HexColor('#19677b'),wordWrap='CJK'),'small':ParagraphStyle('small',fontName='CN',fontSize=8.5,leading=13,spaceAfter=6,wordWrap='CJK')}

REPLACEMENTS=json.loads((O/'editorial_replacements.json').read_text(encoding='utf-8'))
def clean(t):
 for old,new in REPLACEMENTS:t=t.replace(old,new)
 return t

story=[];textrecords=[];figpages=[]
def p(t,s='p'):
 textrecords.append(clean(str(t)));story.append(Paragraph(html.escape(clean(str(t))).replace('\n','<br/>'),styles[s]))
def page(template='portrait'):story.extend([NextPageTemplate(template),PageBreak()])
def tab(rows,widths):
 t=Table([[Paragraph(html.escape(clean(str(x))),styles['small']) for x in r] for r in rows],colWidths=widths,repeatRows=1,hAlign='LEFT');t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e7f0f3')),('GRID',(0,0),(-1,-1),.3,colors.HexColor('#b6cbd2')),('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6),('LEFTPADDING',(0,0),(-1,-1),7),('RIGHTPADDING',(0,0),(-1,-1),7)]));story.append(t);story.append(Spacer(1,9))

def overview_image(id,maxh=520):
 fn=(O/'figures' if id=='S_C' else H/'lung_cohort_bridge_20261017_v25/figures')/(id+'.png');im=PILImage.open(fn);k=min(505/im.width,maxh/im.height);story.append(Image(str(fn),width=im.width*k,height=im.height*k,hAlign='CENTER'));story.append(Spacer(1,7))
def cohort_pages():
 page();p('队列关系：6个主队列样本与4个补充样本构成联合分析','h2')
 p('“联合10样本”指GSE314943与GSE183682的样本级分层模型，共5个对照和5个矽尘样本。两个队列的细胞图谱分别展示；联合统计在相同CD4状态内汇总计数，并控制队列与时间层。')
 overview_image('S_A')
 p('示意图A。主队列GSE314943用于全肺及核心CD4图谱；GSE183682提供D7和D56的4个补充样本。三个模型层为GSE314943-D140、GSE183682-D7和GSE183682-D56。联合统计使用CD4样本级数据，巨噬来源的表达覆盖证据仍来自主队列。GSE274372的12个时间样本单独分析，不计入该10样本模型。','small')
 page();p('补充队列GSE183682：联合分析中4个样本的CD4图谱','h2')
 p('GSE183682纳入1,413个Core CD4细胞，包括742个Naive、398个ActCD4和273个Treg。下图按已有队列内坐标展示细胞状态，并列明参与联合模型的全部4个样本及细胞数量。')
 overview_image('S_B')
 p('示意图B。GSM5567501/502分别为D7对照/矽尘，GSM5567503/504分别为D56对照/矽尘。柱图为细胞数量描述，不以细胞数代替样本量；每时间点1+1，不作单独推断性检验。两队列UMAP的坐标独立计算，不能直接比较坐标位置或距离。','small')
def early_framework():
 page();p('候选基因筛选策略：三项证据交集与综合评价','h2')
 p('差异分析提供候选变化证据，随后结合目标损伤/死亡功能注释及同一T细胞状态的巨噬来源信号依据。三项证据交集得到55个基因，再与补充候选形成90个基因的评价集合。创新性与实验可行性参与后续分层，用于确定预实验优先顺序。')
 overview_image('S_C')
 p('示意图C。圆形交叠仅用于说明筛选逻辑，面积不代表集合大小；55为按细胞状态和对比匹配后去重的三项交集基因数，不完整列示各互斥区域数量。功能维度来自GO/Reactome注释及直接功能文献，不以本项目GO富集P值替代单基因功能依据。逐步阈值、集合合并和候选去向详见后文筛选章节。','small')

p('矽肺免疫微环境与\nCD4⁺ T细胞候选研究','h');p('单细胞图谱、功能程序、细胞通讯及实验候选整合报告','h2');p('2026年10月17日｜研究报告')
p('研究问题与数据','h2');p('围绕矽肺背景下巨噬细胞如何影响T细胞表型，以及这种变化是否伴随线粒体损伤和焦亡开展研究。分析材料为小鼠公开单细胞数据：主队列GSE314943为3个对照与3个矽尘样本、37,155个全肺细胞，其中CD4细胞2,311个；与GSE183682的4个样本构成分层联合模型，另以GSE274372的12个样本评价时间背景。')
p('矽尘条件下T细胞出现什么变化','h2');p('Naive比例下降、ActCD4及Treg比例上升，组成检验P分别为0.10、0.20、0.10。Treg中IL1R2和SKIL呈疾病相关上调；功能富集显示翻译、免疫响应等程序变化。组成变化方向与表达差异分别构成表型和分子证据。')
p('巨噬信号与候选在哪些环节相联系','h2');p('通讯预测与NicheNet提供巨噬来源信号及T细胞响应候选：IL-1与Treg的IL1R2呈诱饵结合关系；SKIL具备趋化因子等配体的靶先验；TTC39C主要依据凋亡与线粒体应激程序关联。主队列ActCD4、Treg焦亡相关评分P为0.30、0.20，当前未检出显著组间差异。')
tab([['候选与切入环节','主要数据依据','建议研究读数'],['IL1R2 / Treg表型，优先','联合上调logFC=3.290，P=1.66×10^-6，FDR=0.000662','Treg状态及IL-1环境响应'],['TTC39C / 损伤探索','ApoExec rho=0.943，P=0.00347，FDR=0.0556；联合DE P=0.163','存活及凋亡反应'],['SKIL / 状态调控备选','TAct rho=0.965，P=0.00208，FDR=0.0744','Treg状态稳定性']],[115,250,140])
p('首轮实验连接哪些环节','h2');p('先验证巨噬—T细胞模型与候选可测性，再选一个可实施候选，采用四组设计检验其是否改变刺激响应；先确认膜电位、存活及焦亡相关执行读数的可测性。由结果决定候选推进、上游定位与救援，检验表型变化与线粒体、死亡终点之间的联系。','small')
page();p('研究设计与数据范围','h')
tab([['数据集','用途','分析单位与边界'],['GSE314943','主发现队列；全景、CD4状态、巨噬来源信号','3个Saline+3个Silica样本；37,155个全肺细胞；Core CD4共2,311个'],['GSE183682','疾病方向补充及联合DE/关联','D7、D56各1个Saline+1个Silica；单时间点仅描述方向'],['分层联合疾病模型','GSE314943的6个样本＋GSE183682的4个样本','5个对照＋5个矽尘；控制队列-时间层，包含主队列，不能算另一独立复现'],['GSE274372','真实时间背景与程序联系补充','D0/D7/D28/D56各3个样本；TLike共494细胞；每样本至少10个TLike细胞的敏感性分析纳入11个样本']],[105,180,220])
p('统计与筛选原则','h2');p('样本级推断以独立生物学样本/文库为单位，细胞、轨迹区间和技术重复不替代样本量。样本少于10的结果按探索性证据呈现。原始目标效应P<0.05用于名义信号筛选，FDR保留在预定义同类假设家族内，结合效应量、方向及样本一致性评价候选证据。')
p('综合差异表达、功能注释、细胞通讯网络及候选相关证据，对90个基因进行评价。其中52个符合计算筛选条件，37个仍需补充证据，1个未达到名义差异筛选条件。进一步选择9个基因、10个基因-细胞状态组合开展探索性分析，并结合功能文献提出3个优先预实验对象。该选择发生于数据分析后；计算筛选资格不等同于跨队列一致差异或完整机制证据。')
p('图示阅读顺序','h2');p('示意图A/B说明队列关系与补充队列图谱；示意图C在差异分析前介绍交集框架。图1-2免疫全景；图3-4 CD4亚群；图5-6差异表达；图7-8 GO/KEGG；图9-10 GSEA；图11-12细胞通讯；图13-14转录因子网络；图15-16焦亡/凋亡评分；图17-18状态轨迹；图19-20 NicheNet；附录A图21-22候选筛选流程；图23-24候选与验证。')
cohort_pages()
page();p('分析方法','h')
p('图谱、亚群与差异表达','h2');p('对主队列质控后的细胞进行谱系注释，并采用UMAP/Harmony坐标展示，按全肺谱系及Core CD4状态展示。细胞比例以每样本Core CD4为分母，用6样本全部20种组标签分配进行双侧精确检验，3个状态构成一个BH家族。差异表达采用原始计数样本聚合、完整文库分母、表达过滤、TMM归一化及edgeR准似然模型；联合模型包含队列/时间层和疾病效应。各差异对比均报告完整受检基因集合内的P值及FDR。')
p('GO/KEGG与GSEA','h2');p('使用联合模型中ActCD4及Treg的完整差异表达结果进行功能分析。GO BP采用org.Mm.eg.db的GOALL注释，KEGG采用小鼠基因-通路映射。基因符号经统一处理；以可映射受检基因为ORA背景。原始P<0.05的上、下调基因分别作超几何检验，集合限制10-500个背景成员，零命中集合保留在完整校正家族中。')
p('GSEA对全部受检基因按sign(logFC)×sqrt(QL F)排序，同分按基因名固定顺序。fgseaMultilevel使用10-500成员、seed=20261017、nPermSimple=1000、eps=1e-8。其P为排序基因集合富集模型统计，非新增独立样本处理效应检验。每状态/数据库分别BH校正；13条结果达到P估计下限，按P≤1e-8解释。代表图选择规则固定在图注及绘图数据中。')
p('通讯与转录因子网络','h2');p('通讯以数据库配体-受体相互作用及发送端和接收端的组级表达覆盖为依据，按必要受体组分评估表达相容性。图11通过弦图展示巨噬细胞来源与CD4状态之间满足表达覆盖的相互作用数量，图12展示Treg相关配体-受体候选连接。TF分析采用SCENIC网络，移除regulon中与候选基因及程序成员重叠的基因，按独立样本平均AUCell；代表网络按保留边权重每TF取前2个靶标。')
p('GSEA曲线使用enrichplot::gseaplot2绘制。Raw表示完整细胞集；BG1表示background_fraction≤1%的低背景负担细胞子集。该比例为预定义上皮、髓系及免疫球蛋白标志基因UMI占细胞完整文库UMI的比例，不是环境RNA污染率估计。完整标志集合与计算定义见补充材料。','small')
page();p('评分、关联与状态轨迹方法','h')
p('样本程序评分','h2');p('对样本s与状态k，C_gsk为该样本该状态细胞的基因原始计数和，L_sk为完整原始文库UMI之和。E_gsk=10000×C_gsk/L_sk；候选表达G_gsk=ln(1+E_gsk)，模块评分M_msk为保留映射成员上mean_g ln(1+E_gsk)。先样本汇总、后归一化和取对数。零表达保留；缺失映射不填零；至少3成员且映射覆盖≥50%才评价。9候选从模块移除，实际重叠者DDIT4从MitoISR/Glycolysis移除。')
p('Core等权状态标准化敏感性给Naive/ActCD4/Treg各1/3权重：候选先平均各状态E再log，模块和AUCell平均各状态分值。等权状态标准化与基于六样本平均状态比例的加权估计（0.466639/0.309979/0.223382）对应不同估计目标。轨迹分析采用逐细胞匹配对照评分后按样本聚合，与汇总计数评分的数值尺度不同。')
p('样本关联及精确置换','h2');p('联合秩关联在队列/时间层设计下残差化，块内排列6!×2!×2!=2880。疾病调整模型包含疾病项，组内有效排列3!×3!=36，第二队列4个单例单元不能贡献组内交换。联合组间比较保持层内组别数，共80种分配。真实时间关联控制日龄，4组各3样本时有1296种组内排列，去掉1个低细胞样本后432种。每类模型保留自身完整FDR家族。')
p('候选-表型关联的BH家族由同一队列、细胞集、调整模型与检验类型定义。联合疾病调整模型包含10个候选-状态组合×5个表型程序，共50项；损伤程序、TF、配体及其他模型另行校正。完整家族成员随统计表提供。','small')
p('状态拟时序与NicheNet','h2');p('状态轨迹采用主队列Naive/ActCD4的4套DPT排序，组合K15/K30邻居参数及两种根设置；轨迹分析范围不包含Treg。每状态内将DPT缩放，并在每个样本中调整log文库量后拟合表达斜率；以6个样本斜率作64种符号翻转检验。图中4个分位区间只用于展示。UMAP主干将K15邻域图中的细胞按DPT分为12个分位组，汇总组间连接强度并除以两组细胞数乘积的平方根，以最大生成树显示主要邻域连接；箭头沿既有DPT增大方向。')
p('NicheNet整合配体-靶基因先验、受体表达条件、差异目标集及受检基因背景，展示联合疾病上调目标集的配体活动及当前候选路线。活动Pearson及网络权重作为模型排序指标，保留各状态实际目标集和背景。')
# Results lead paragraphs, one figure per page with explicit source and method caption.
results={
'F01':'全肺注释识别出B、T、巨噬/单核、中性粒、树突等免疫群，并同时保留上皮、间质及血管细胞作为组织环境。后续巨噬-CD4分析均在这一细胞来源框架内展开。',
'F02':'T细胞相关CD3D/TRAC、B细胞相关MS4A1/CD79A及髓系相关标志呈现不同分布，支持主要谱系区分。图中广泛出现的组织相关转录背景也说明候选定位需要结合亚群身份及背景敏感性。',
'F03':'Core CD4由933个Naive、801个ActCD4和577个Treg构成。疾病组更多分布于ActCD4/Treg状态区域，样本层面的组成量化见下一图。',
'F04':'疾病组Naive平均比例降低46.85个百分点，ActCD4增加20.85个百分点，Treg增加26.00个百分点；对应精确P为0.10、0.20、0.10，FDR为0.15、0.20、0.15。方向清楚，但当前样本量下未达到名义显著。',
'F05':'Treg中的IL1R2与SKIL在联合疾病模型呈上调信号。TTC39C的联合疾病DE未达到名义显著，因此后续推荐它主要依据凋亡程序联系，而非稳定跨队列差异。',
'F06':'IL1R2/Treg联合logFC=3.290、P=1.66×10^-6、FDR=0.000662；SKIL/Treg联合logFC=1.438、P=2.15×10^-5、FDR=0.00314。FOXP1/Treg为下调状态候选，支持其作为状态重塑候选。',
'F07':'Treg的下调集合突出细胞质翻译/核糖体相关过程，同时上调集合包含白细胞活化；ActCD4的部分免疫相关过程获得富集支持。具体含义按命中成员解释，涉及免疫球蛋白或广泛背景的条目另表标注。',
'F08':'Treg下调差异基因在核糖体相关KEGG集合中富集，并出现上调细胞因子-受体相互作用线索。ActCD4细胞黏附相关集合原P=0.000284、FDR=0.0804，为名义信号。',
'F09':'全基因排序与阈值差异分析共同显示翻译/核糖体相关集合偏向疾病下调端。部分免疫相关集合偏向上调端，GO/KEGG同名或共享成员的结果不计作独立验证。',
'F10':'代表核糖体通路在Treg和ActCD4均为负NES。此处展示完整排名中的成员累积位置，可与ORA方向互相参照；它反映转录集合变化，不直接测量蛋白合成速度。',
'F11':'在所设表达覆盖规则下，巨噬来源至CD4各状态均存在表达相容入口，疾病组满足覆盖的相互作用数量较多。此结果用于确定可检验信号范围，不能把计数变化解释成实际通讯强度增加。',
'F12':'AMlike中的IL1A、IL1B、CCL2和TNF按独立样本呈现不同表达分布。后续上游验证应保持同一发送范围，结合完整受体和近端响应选择具体路线。',
'F13':'去重叠regulon的AUCell呈现CD4状态与样本间差异，为候选相关程序提供调控环境。优先候选RNA、TF靶网络及活性评分是互补证据，不能相互替代。',
'F14':'FOS/JUN、NFAT、NF-κB、FOXP1与TRP53等转录因子网络可连接具体靶成员。图中代表边用于展示网络构成；完整网络成员用于靶标选择及结果解释。',
'F15':'ActCD4与Treg的PyroExec组间效应方向为上升，但主队列精确P分别为0.30与0.20；Inflam对应P均为0.90。当前组间检验不支持将疾病相关CD4状态变化直接等同于已发生焦亡。',
'F16':'ActCD4的ApoExec差值为+0.0404、P=0.10、FDR=0.274；Treg差值+0.0352、P=0.20、FDR=0.433。基因-程序相关与组间评分差异是不同问题，前者的阳性不替代后者。',
'F17':'四套DPT的细胞排序相关为0.904-0.999，支持当前Naive/ActCD4状态排序对所测设置较稳定。其适用范围为转录状态连续性，Treg不在这条轨迹内。',
'F18':'主要排序中，Naive的PyroExec斜率P=0.0625、ApoExec P=0.125；ActCD4对应P=0.156与0.188。当前数据不支持用这些程序曲线证明沿轨迹发生持续焦亡/凋亡增强。',
'F19':'TNF、TGFB1及若干免疫调节配体在联合疾病上调目标集合中获得模型排序线索。配体优先级由目标集、背景及受体条件共同决定，应在具体发送细胞范围中继续定位。',
'F20':'IL1R2和SKIL在联合名义差异及网络条件下具有代表性配体连接。TTC39C未满足本图对应的联合差异筛选条件，未显示上游连接；其主队列结果、配体-靶基因先验及程序关联分别评价。',
'F21':'IL1R2/Treg阳性细胞145个，均共检至少一个CD3D/CD3E/TRAC。\nTTC39C/ActCD4为125/125；SKIL/Treg为319/331。背景≤1%保留的阳性细胞分别为141、119和300。候选表达方向在完整细胞集/低背景负担细胞子集及样本留一描述中一致。',
'F22':'IL1R2与Inhib、TTC39C与ApoExec、SKIL与TAct均有联合名义关联。推荐以各自最直接的细胞状态/功能问题开展预实验，并保留疾病调整及外部方向结果决定下一步投入。'}

results['F09']='Treg的KEGG富集曲线同时呈现疾病上调端与下调端的通路分布。核糖体相关集合偏向下调端；其他通路需结合其具体成员及细胞状态解释。'
results['F10']='ActCD4的多通路曲线显示不同功能集合在完整基因排序中的富集位置。正、负NES分别对应疾病上调端和下调端，图中并列报告每条通路的P值及FDR。'
results['F11']='巨噬细胞与Naive、ActCD4和Treg之间均存在满足所设表达覆盖的配体-受体候选。分组弦图展示这些候选的状态分布，并保留各发送范围的独立定义。'
results['F12']='Treg候选关系包括IL1A/IL1B-IL1R2、CCL2-CCR2、CXCL10-CXCR3、CD200-CD200R1及SEMA3E-PLXND1，可用于安排后续来源、受体与近端响应验证。'


results['F17']='四套DPT排序的相关为0.904-0.999。叠加的邻域主干连接Naive/ActCD4的主要状态区域，帮助阅读DPT方向；线条来自原始邻域连接的描述性抽象，不作为真实分化过程或损伤进展的证明。'
results['F20']='矩阵直接显示IL1R2的诱饵受体关系与SKIL的配体-靶先验连接。TTC39C在本展示对应的联合差异及网络条件下没有合格连接，其程序关联证据仍单独评价。'


results['F23']=results['F21']
results['F24']=results['F22']
results['F21']='主要差异模型共涉及9,091个去重基因；2,765个在至少一个目标对比中达到原始P<0.05。与目标功能集合相交后保留111个，再匹配相同接收状态的巨噬来源信号得到55个。结合补充候选形成90基因评价集合。'
results['F22']='综合评价最终支持52个计算候选，37个保留为证据待补，1个未达到当前名义差异条件。重点9基因来自完整90基因评价集合；TTC39C虽然处于功能证据待补层，仍可凭凋亡程序关联及功能实验先例进入探索性预实验。'
def screening_methods():
 page();p('附录A｜筛选阈值、集合及候选去向','h')
 p('筛选首先整合差异、目标损伤/死亡功能和巨噬来源信号依据三项证据，再进行细胞身份与作用角色评价。创新性及实验可行性用于确定研究方向，不要求所有候选在四项指标上同时通过一个统一阈值。')
 tab([['维度','实际判定条件','参考证据与解释'],['差异表达','有效样本级模型的目标效应原始P<0.05；按基因-状态-对比保留记录','FDR、|log2FC|≥0.5、跨队列方向及敏感性用于评价强度；不要求所有队列同时显著'],['目标损伤/死亡功能','GO/Reactome内源性凋亡、线粒体损伤相关过程或焦亡集合；明确的直接功能文献可扩展','两类机制取并集；不要求基因同时参与凋亡和焦亡。一般凋亡/状态证据不足以证明目标细胞自身损伤'],['巨噬来源信号依据','配体和完整受体组分在同一条件及接收状态内满足覆盖；候选为受体组分或具有对应配体-靶先验','每样本至少10个目标细胞，阳性比例≥5%，且组内至少一半合格样本满足。网络先验不等于功能激活'],['创新性与可行性','逐项比较已有机制、细胞背景、研究空间与可观察实验表型','作为分层依据。已有发表不自动排除，未检索到文献不自动证明创新；结合机制与细胞背景判断研究空间']],[80,215,210])
 p('交集按同一细胞状态匹配；时间证据单列','h2')
 p('主要DE范围包括主队列4个细胞范围、联合队列4个范围及时间队列3个时间对比和1个总体时间检验，共12个家族。名义差异在这些目标对比间取并集，时间效应与疾病效应分别解释；功能和网络匹配完成后再按基因去重，避免跨状态拼接。','small')
 p('初筛与综合评价的网络范围','h2')
 p('初筛网络分支采用校正AUPR>0的合格配体，按活动排序取前30个，每配体取表达背景内前250个靶基因；直接受体分支单独保留。综合评价时扩展检查有正向活动的合格配体，网络排名作为参考，不因超过30/250而统一排除。','small')
 p('补充集合、资格分层与预实验推荐','h2')
 p('42个补充候选围绕T细胞状态、细胞死亡及相关调控问题纳入，与55个交集基因合并后统一评价；来源及成员逐项列入附表。细胞来源不清或主要作用于其他细胞的输出效应基因进入待补层。重点分析选择IL1R2、DDIT4、SKIL、FOXP1、GSN、AHNAK、TTC39C、MDFIC、ST6GALNAC3；预实验顺序进一步依据功能文献、程序关联和可操作性确定。','small')

results.update({'F10': 'ActCD4的抗原加工与呈递集合正向富集（NES=1.903，P=0.000247，FDR=0.0174），核糖体集合负向富集（NES=-1.965，P=5.81×10^-7，FDR=0.000164）。这些方向将活化状态与免疫响应、翻译程序的变化联系起来。', 'F13': 'Treg中NFKB1的样本平均AUCell由对照0.0600升至矽尘0.0779，FOXP1由0.1539降至0.1354。这里呈现调控程序的描述性差异，为活化与状态变化提供背景，不新增显著性判定。', 'F14': '代表网络中JUNB连接GADD45B与NFKBIA，FOS连接JUNB与GADD45B，形成应激与活化相关靶网络线索。各边为SCENIC调控先验，后续可据此选择响应读数。', 'F15': 'ActCD4和Treg的焦亡相关评分均呈上升方向，精确P分别为0.30和0.20；一般炎症评分的组间P均为0.90。本队列未检出显著组间差异，焦亡作为后续实验的独立观察终点。', 'F16': 'ActCD4的凋亡相关程序ApoExec差值为+0.0404、P=0.10、FDR=0.274；Treg差值为+0.0352、P=0.20、FDR=0.433。两状态方向均向上，未达到名义显著；后续结合候选关联检验具体损伤环节。', 'F19': '联合疾病上调目标集的校正AUPR排序中，Treg以ICAM1（0.0552）和TNF（0.0547）居前；ActCD4以WNT4（0.0241）和IL10（0.0237）居前。排序反映各状态响应集合的模型解释度，候选基因的具体连接见下一图。', 'F20': 'IL1R2的IL-1诱饵结合与SKIL的CXCL10、CCL2、SPP1靶先验属于不同证据类型。TTC39C未纳入本网络路线评价；H2-DM相关记录单列为抗原呈递信息。', 'F21': '三项证据交集包含55个基因，与42个研究问题相关补充候选合并，重叠7个，得到90个综合评价基因。', 'F22': '90个基因形成52个计算支持、37个证据待补及1个未达到名义差异条件的分层。9个重点基因来自完整90基因集合，其中5个计算支持、4个证据待补；3个预实验建议结合功能关联、文献及可行性确定。'})
for m in manifest:
 if m['id']=='F05':early_framework()
 if m['id'] in ['F21','F22']:continue
 template='portrait';page(template);p(m['module']+'｜'+m['title'],'h2');p(results[m['id']],'p')
 im=PILImage.open(m['file']);maxw=505;maxh=540;scale=min(maxw/im.width,maxh/im.height);story.append(Image(m['file'],width=im.width*scale,height=im.height*scale,hAlign='CENTER'));story.append(Spacer(1,6));p('图'+str(int(m['id'][1:]))+'。'+m['caption'],'small')
page();p('三候选证据连接｜从巨噬来源到T细胞响应','h')
p('三者对应同一研究主线的不同环节。下表把表达覆盖、网络预测、程序关联和首轮实验问题分开列明；其中AllMac为全体巨噬来源，AMlike为其相关亚群范围，二者不作独立证据相加。')
tab([['环节','IL1R2 / Treg','TTC39C / ActCD4','SKIL / Treg'],['巨噬来源信号','AllMac的IL1A/IL1B与Treg IL1R2在两组满足覆盖','当前联合差异网络未纳入该候选；上游未定位','AllMac的CXCL10、CCL2、SPP1与相应Treg受体在两组满足覆盖'],['受体及网络关系','IL1R2为诱饵受体；结合不等于胞内传导','首轮以巨噬刺激体系评价，不指定配体','CXCR3、CCR2、CD44或ITGAV+ITGB1；SKIL靶排名192、235、155'],['候选/表型依据','疾病上调，Inhib关联及Treg功能文献','ApoExec联合rho=0.943，P=0.00347，FDR=0.0556','TAct联合rho=0.965，P=0.00208，FDR=0.0744'],['损伤依据','PyroExec P=0.988；当前无直接阳性焦亡关联','TLike中MitoISR时间调整rho=0.733，P=0.0185，FDR=0.567','ApoExec疾病调整P=0.917；当前不支持该损伤关联'],['首轮检验','巨噬刺激下的Treg状态及IL-1环境响应是否受IL1R2影响，并伴随线粒体/死亡变化？','候选是否改变巨噬刺激后的存活/凋亡，并伴随膜电位和独立焦亡读数变化？','候选是否改变巨噬刺激后的Treg状态，是否伴随损伤终点变化？']],[75,143,143,144])
p('SKIL的配体-靶排名来自联合Treg上调目标集先验，表达覆盖以主队列发送端和接收端核对；它不是来源特异的处理效应或已验证因果轴。AMlike的配体表达不能替代同一来源的整条实验验证。')
p('H2-DM相关数据库记录另列抗原呈递网络信息，不作为普通表面配体-CD4作用轴推荐。TTC39C的网络NA不否定其凋亡程序关联。完整连接、组别和数据库来源见补充材料。','small')
page();p('候选比较与外部方向','h')
tab([['候选','当前项目证据','关键未决环节'],['IL1R2 / Treg','原六文库logFC=3.115，P=4.89×10^-5，FDR=0.00853；联合DE及Inhib支持；ENO1关联P=0.000694、FDR=0.0285','ENO1疾病调整P=0.111；PyroExec P=0.988；不能用RNA正相关证明酶活抑制'],['TTC39C / ActCD4','ApoExec联合rho=0.943、P=0.00347、FDR=0.0556；疾病调整rho=0.922、P=0.0278、FDR=0.556','联合DE P=0.163、FDR=0.803；D7方向受归一化影响，D56反向；尚无明确上游配体，先检验巨噬刺激下的候选效应'],['SKIL / Treg','TAct联合rho=0.965、P=0.00208、FDR=0.0744；疾病调整rho=0.897、P=0.0278、FDR=0.8681','ApoExec P=0.0747、调整P=0.917；当前数据不支持疾病调整后的凋亡程序关联；AMlike CCL2不能跨来源拼接']],[90,210,205])
p('SKIL疾病调整关联为名义显著，FDR=0.8681来自上述50项候选-表型检验家族。该结果支持探索性状态线索，未达到校正后显著。','small')
p('真实时间队列提供不同问题的补充证据','h2');p('GSE274372的TLike时间分析中，TTC39C与MitoISR在日龄调整后rho=0.733、P=0.0185、FDR=0.567，为探索性信号。它不是ActCD4-ApoExec同一关联的独立复现。时间、疾病、细胞状态和损伤终点分别评价。')
p('其余候选继续按证据分层保留','h2');p('FOXP1优先解释状态重塑，Core-PyroExec关联在等权状态标准化后P=0.494；MDFIC为损伤关联备选，调整后rho=0.950但P=0.0833；DDIT4保留文献驱动假说，项目内PyroExec P=0.751；GSN保留备选。AHNAK/ST6GALNAC3保留氧化代谢或应激转录线索。')
page();p('功能文献为实验提供依据与方向约束','h')
p('IL1R2：以Treg功能为首轮切入点 [1-4]','h2');p('2022年Treg特异Il1r2缺失研究支持其参与肿瘤Treg克隆扩增和效应特征；2024年Tfr研究支持其调节特定免疫环境中的IL-1作用。2025年IL1R2-ENO1-糖酵解/焦亡链来自巨噬细胞及脓毒症。矽肺相关研究也已有其他细胞对象的Il1r2表达报道。研究空间应聚焦矽肺巨噬环境对CD4/Treg及其功能的影响，允许识别保护性反馈。')
p('TTC39C：优先确认干预后的存活与凋亡方向 [5]','h2');p('人A549/NCI-H1299细胞的shRNA实验中，TTC39C降低后Annexin V阳性增加。项目中TTC39C与ApoExec正相关尚不足以证明其促进凋亡；首轮要判断在CD4中是驱动、补偿还是其他背景依赖作用。')
p('SKIL：同时观察SnoN蛋白与SKI [6-7]','h2');p('SKIL编码SnoN，SKI为另一基因。小鼠Treg分化研究中，在Arkadia缺失背景下，双靶向Ski/SnoN可救援，单独靶向不能完成该救援，提示特定条件下的冗余。SKIL RNA上调不等同SnoN蛋白或功能增强。PM2.5/ROS/SnoN与肺纤维化已有研究联系，广义关联本身不作为全新机制。')
p('蛋白检测与基因干预的可行性','h2');p('IL1R2具有小鼠Tfr流式4E2检测先例及厂家小鼠适用性资料。TTC39C有厂家小鼠蛋白WB资源，但原代CD4中的检测和干预效力仍需确认。SnoN原T细胞研究有抗体验证；当前替代试剂需核对物种、SnoN/SnoN2识别和阴性对照，不能把仅人源验证的产品直接视为小鼠适用。检测资源及出处见实验可行性表；试剂在原代目标细胞中的特异性和效力需通过预实验确认。')
page();p('建议验证路径｜连接巨噬环境、T细胞表型与损伤','h')
p('目前可推荐小规模模型及候选可测性验证。IL1R2优先核查Treg与IL-1环境响应，TTC39C关注存活及线粒体状态，SKIL关注T细胞状态调控。优先级表示验证价值，候选功能投入由模型响应、检测条件和实际资源决定。')
tab([['阶段','建议研究问题','推进依据'],['优先验证','巨噬条件是否改变目标T细胞状态与存活？候选RNA/蛋白、线粒体及焦亡相关读数是否可测？','在对应亚群小规模筛查；选择一个检测与干预最可实施的候选'],['有功能结果后','候选干预是否改变巨噬刺激效应？独立干预、回补及上游阻断能否支持作用？','可重复效应、充分干预及可解释基础毒性'],['扩展研究','动物/患者材料、复杂双基因、蛋白相互作用或深度代谢问题','按研究深度、首轮结果、平台、周期及预算另行确定']],[85,245,175])
tab([['候选','建议优先回答的问题','与损伤终点的联系'],['IL1R2 / Treg','刺激是否改变Treg数量/状态及IL-1相关环境响应？抑制试验作为可选功能读数','观察膜电位和死亡；常规抑制试验不作为唯一成败标准'],['TTC39C / 活化CD4','蛋白与干预是否可测？候选是否改变存活/凋亡？','分别观察线粒体及焦亡相关执行层，不预设作用方向'],['SKIL / Treg','SnoN蛋白和状态响应是否改变？成熟Treg稳定性为项目探索问题','状态变化与膜电位/死亡并行观察；不直接等同诱导Treg文献模型']],[100,225,180])
p('模型阶段即安排有限的线粒体与焦亡可测性检查；具体读数取决于目标细胞量和平台。有可测响应后，在选定候选的功能试验中预设主要结局及并行探索终点。线粒体变化与焦亡分别判定，先后关系由时间和干预实验判断。')
p('本节为研究方向与可选验证路径。执行范围确认后，再由实验室确定细胞来源、分选/递送、剂量、时长、独立重复数与费用。四组示例、检测注意点和必要对照见附录B；选用实验所需的对照应随该实验落实。','small')
page();p('讨论：证据边界与研究空间','h')
p('当前证据支持预实验，机制连接由干预检验','h2');p('疾病差异、程序关联、网络先验和干预功能是不同证据层级。IL1R2具备较清楚的疾病差异和Treg功能研究基础；TTC39C提供损伤相关程序切入点；SKIL具备Treg状态和调控背景。三者均有预实验价值，但不对应已量化的成功概率。')
p('样本与模型边界','h2');p('分析材料为小鼠；主队列为3+3个样本；第二队列各时间点1+1，联合分析含原队列；时间队列为TLike范围。精确枚举保证既定排列空间计算完整，不证明可交换性或模型假设成立。MDFIC等强效应而P有限的结果不能全部归因于分辨率。')
p('生物学与技术解释','h2');p('转录评分不能替代线粒体膜电位、呼吸、蛋白裂解或死亡实测。富集结果中的Coronavirus disease、synapse等名称可能由共享核糖体/翻译成员驱动，不表示感染或神经细胞功能；免疫球蛋白/组织背景成员需结合背景敏感性解释。图谱中Like及Unresolved标签表示细胞身份仍存在不确定性。')
p('PROS1-AXL的胶质瘤文献背景不能直接用于小鼠CD4。TYRO3/MERTK的转录表达未满足所设覆盖条件；RNA低检出不是蛋白不存在的证明。其他网络边也仅保留计算支持，不视为已经验证。状态DPT是横断面转录排序，不能代替谱系追踪或真实损伤时间。')
p('创新性落在具体细胞关系和作用环节','h2');p('IL1R2-ENO1-焦亡、SKIL/SnoN-TGF-β-Treg及TTC39C-肿瘤细胞存活均已有研究。新的问题应落在矽肺巨噬环境如何改变CD4/Treg候选及功能，以及候选是否参与损伤耐受或细胞死亡。定向检索未取得完整链直接功能证明，不等同全球首次声明。')
page();p('参考文献','h')
refs=[('[1]','IL1R2 increases regulatory T cell population in the tumor microenvironment by enhancing MHC-II expression on cancer-associated fibroblasts. JITC, 2022.','https://pmc.ncbi.nlm.nih.gov/articles/PMC9438093/'),('[2]','Loss of T follicular regulatory cell-derived IL-1R2 augments germinal center reactions via increased IL-1. JCI Insight, 2024.','https://insight.jci.org/articles/view/174005'),('[3]','Critical Role of IL1R2-ENO1 Interaction in Inhibiting Glycolysis-Mediated Pyroptosis for Protection Against Lethal Sepsis. Advanced Science, 2025.','https://pmc.ncbi.nlm.nih.gov/articles/PMC12533375/'),('[4]','Epithelial NLRP3 drives silica-induced lung injury and fibrosis through IL-18 and pro-fibrotic neutrophil recruitment. Particle and Fibre Toxicology, 2026.','https://pmc.ncbi.nlm.nih.gov/articles/PMC13214385/'),('[5]','Ttc39c is a potential target for the treatment of lung cancer. BMC Pulmonary Medicine, 2022;22:391.','https://pmc.ncbi.nlm.nih.gov/articles/PMC9615393/'),('[6]','Arkadia-SKI/SnoN signaling differentially regulates TGF-β-induced iTreg and Th17 cell differentiation. JEM, 2021.','https://pmc.ncbi.nlm.nih.gov/articles/PMC8421263/'),('[7]','PM2.5 regulates TGF-β1/Smads-mediated pulmonary fibrosis via ROS/SnoN in vitro and in vivo. Int Immunopharmacol, 2025.','https://pubmed.ncbi.nlm.nih.gov/40154177/')]
refs.extend([('[8]','Galluzzi et al. Molecular mechanisms of cell death: recommendations of the Nomenclature Committee on Cell Death 2018. Cell Death & Differentiation, 25:486-541.','https://www.nature.com/articles/s41418-017-0012-4'),('[9]','Thermo Fisher Scientific. TMRM technical guidance and membrane-potential controls.','https://www.thermofisher.com/order/catalog/product/T668/faqs')])
refs.append(('[10]','Reactome. H2-DMb1, lysosomal membrane; antigenic peptide loading. R-MMU-2130673.','https://www.reactome.org/content/detail/R-MMU-2130673'))
refs.extend([('[11]','Acute Silica Exposure Triggers Pulmonary Inflammation Through Macrophage Pyroptosis: An Experimental Simulation. Front Immunol, 2022.','https://doi.org/10.3389/fimmu.2022.874459'),('[12]','Pore-forming activity and structural autoinhibition of the gasdermin family. Nature, 2016.','https://doi.org/10.1038/nature18590'),('[13]','Chemotherapy drugs induce pyroptosis through caspase-3 cleavage of a gasdermin. Nature, 2017.','https://doi.org/10.1038/nature22393')])
for n,title,url in refs:p(n+' '+title+'\n'+url,'small')
p('补充材料','h2');p('完整统计、图值及来源索引随本报告提供，见附录C及配套材料。','small')
screening_methods()
for m in manifest:
 if m['id'] not in ['F21','F22']:continue
 page();p(m['module']+'｜'+m['title'],'h2');p(results[m['id']]);im=PILImage.open(m['file']);k=min(505/im.width,540/im.height);story.append(Image(m['file'],width=im.width*k,height=im.height*k));p('图'+str(int(m['id'][1:]))+'。'+m['caption'],'small')
page();p('附录B｜建议研究问题与可选实验路径','h')
p('本附录供客户与执行实验室确定研究范围。先验证模型与可测性，再选一个候选进入功能干预；独立干预、救援及机制研究按结果推进。动物、患者材料、双基因、蛋白相互作用及深度代谢属于另行确定的扩展研究。')
p('第一阶段：建立巨噬—T细胞作用体系','h2');p('优先采用与发现数据一致的小鼠体系，选择一种巨噬来源和一个目标T细胞亚群。Treg与活化CD4分别对应不同候选，混合总CD4不能替代亚群证据。若以骨髓来源巨噬细胞预筛，关键作用宜回到与假说匹配的肺来源体系确认。')
p('比较对照与矽尘处理巨噬条件对T细胞状态及存活的影响，同时检查线粒体/死亡相关读数的可测性。条件培养基适于可溶性作用预筛；接触或呈递假说需要相应共培养，培养基试验阴性不能排除全部巨噬作用。')
tab([['模型控制','建议记录或对照'],['矽尘材料','结晶型、颗粒特征、批次与内毒素资料；实际刺激范围在模型内优化'],['颗粒及发送端效应','直接颗粒及残留匹配对照，核实去除效果；记录巨噬细胞死亡、数量与培养条件'],['目标T细胞','用活细胞富集进行功能培养，平行样本验证亚群；活化/支持条件一致，计数与死活结合'],['检测阳性对照','用于确认方法有效；不能将额外强效死亡诱导剂的作用当作巨噬机制复现']],[115,390])
p('第二阶段：候选可测性与投入选择','h2');p('IL1R2优先核查表达、刺激响应与读数条件。资源允许时，三候选可分别在对应亚群小规模检查RNA和蛋白；功能干预优先选择一个可测且可实施对象。原代T细胞须确认递送、蛋白改变及基础存活，不直接移用肺癌细胞的条件。技术失败不等同基因无功能。')
p('实施前确认：可获得细胞与数量、分选/递送平台、抗体适用性、各读数细胞需求、周期与预算；据先导变异和目标效应确定正式独立样本量。','small')
page();p('四组设计示例｜检验候选是否改变巨噬刺激响应','h')
p('适用于已选择候选的功能验证阶段。培养体系确定后，以巨噬刺激条件×接收端T细胞候选干预建立因子设计；各组匹配递送、处理、活化与培养支持。')
tab([['组别','巨噬条件','目标T细胞处理','回答问题'],['A','对照巨噬','非靶向对照','基础参照'],['B','矽尘处理巨噬','非靶向对照','巨噬刺激效应'],['C','对照巨噬','候选干预','干预主效应及基础毒性'],['D','矽尘处理巨噬','候选干预','干预是否改变刺激响应']],[40,115,140,210])
p('同时报告各组值、B−A刺激效应、C−A干预本身效应及交互Δ=(D−C)−(B−A)。连续、比例或计数结局使用相应模型，预先说明Δ在原始尺度还是模型连接尺度上计算。交互未达显著而干预主效应明确，可提示一般T细胞作用；不能宣称巨噬刺激特异调控，也不直接判为无功能。')
tab([['必要控制','用途'],['基因与递送','RNA及蛋白效率；非靶向、载体/递送对照及基础毒性；干预先施加于接收端'],['亚群与归属','目标T细胞门内分析；蛋白检测采用足够纯度细胞；固定透化样本不继续用于功能培养'],['混合上清','LDH、IL-1β可来自巨噬细胞或残留条件培养基，不能据此定位T细胞死亡。[11]'],['细胞数量','活细胞数、死亡/丢失数量与比例同时解释，避免存活选择造成表型假象'],['方法有效性','所选死活、凋亡、膜电位及执行层检测配置适用对照，区分技术不可测与生物学阴性']],[115,390])
p('以独立动物/供体作为重复；分装、复孔、重复时间点不增加独立n。记录巨噬与T细胞各自来源，尽量以同一独立来源覆盖四组，平衡批次、随机分配并盲法分析。','small')
page();p('IL1R2路径｜Treg状态与IL-1环境响应','h')
p('原始功能研究涉及肿瘤Treg—IL-1—成纤维细胞环境及Tfr免疫调节；ENO1/焦亡机制来自巨噬细胞。[1–3] 本项目优先检验这些作用能否迁移到矽尘相关巨噬—Treg环境。')
tab([['层次','建议研究问题与读数'],['可测性','先确认目标Treg的Il1r2 RNA与适用蛋白检测；按假说后续区分表面捕获、可溶性或胞内作用'],['状态与数量','巨噬刺激是否改变活Treg数、FOXP3/CD25等状态？结合存活评价，而非仅看阳性比例'],['IL-1环境','IL1R2干预是否改变环境中的IL-1相关响应？实际响应细胞、传导受体和近端读数由表达及功能证据确认'],['功能选择','等量活Treg抑制试验可用作功能读数，记录Treg回收数与应答T细胞增殖；不将其作为IL1R2价值唯一判断'],['损伤联系','根据平台并行观察膜电位、存活与焦亡相关执行读数，检验状态/环境响应是否伴随损伤变化']],[105,400])
p('两种功能设计回答不同问题','h2');p('回收并洗涤Treg后再做抑制试验，检测处理后保留的功能；持续共培养则包含即时环境作用。前者需控制等量活细胞，后者需区分巨噬/处理对Tconv的直接影响。若缺少IL-1及相关响应细胞，常规抑制试验可能未覆盖微环境作用，这是实验设计判断。')
p('干预类型与后续分支','h2');p('检测抗体、功能阻断、阳性细胞清除和基因敲低回答不同问题；流式可用不等于可阻断。表面IL1R2阻断不等同胞内IL1R2—ENO1干预。若出现相应功能及代谢线索，再考虑ENO1酶活、代谢或蛋白相互作用；RNA正相关不证明酶活抑制。')
p('推进依据','h2');p('先确定所选模型具备相应环境与读数，再预设主要终点。可重复候选效应后推进独立干预/回补；若干预加重损伤，可研究保护性反馈。抑制试验单项阴性不能排除其他IL-1环境作用。')
page();p('TTC39C路径｜先确认原代CD4可测与可干预','h')
p('肺癌细胞研究中，降低TTC39C后凋亡增加，为存活/凋亡方向提供先例。[5] 本项目ApoExec关联属于观察性依据，不预设敲低保护T细胞，也不直接指定p53或某个配体为上游。')
tab([['阶段','建议读数与判定'],['原代目标细胞可测性','明确活化CD4标志与培养状态；确认TTC39C蛋白及试剂适用性'],['干预可行性','目标蛋白改变、递送效率与基础存活达到可解释水平；不能移用癌细胞病毒体积或条件'],['候选功能','四组试验可用预设时点存活/凋亡为主要结局；Annexin V与膜完整性联合，结合caspase-3/7或相应裂解指标'],['线粒体联系','适用膜电位与对照；若稳定改变，再选择独立功能读数，如呼吸或ATP，按活细胞量等归一化'],['独立焦亡观察','模型阶段核查适用gasdermin及裂解/膜损伤可测性；凋亡指标不代替焦亡执行证据']],[120,385])
p('测量与时间顺序','h2');p('活细胞门内记录膜电位，同时记录死亡/丢失数量；未处理和去极化阳性对照用于确认测量有效。[9] 先导阶段可用少量早、中、晚时点区分候选、线粒体与死亡变化；正式比较预先固定主要时点，不选择效应最大的时间作唯一结论。')
p('后续选择','h2');p('死亡或线粒体效应按实测方向解释。可重复效应后，选择独立靶向与回补，再评价执行层和作用顺序。若主要终点阴性而并行终点改变，保留探索结果，在后续预设实验中验证；不追认原主要终点已验证。')
page();p('SKIL路径｜区分成熟Treg状态与诱导分化','h')
p('本项目提出成熟Treg状态调控问题，首先确认巨噬刺激下SnoN蛋白与状态是否可测变化。文献中的Arkadia缺失、诱导Treg分化及Ski/SnoN双靶向救援属于另一实验背景，不能称为本成熟Treg单基因方案的直接复现。[6]')
tab([['层次','建议验证与解释'],['可测性','确认SKIL RNA、SnoN蛋白及实际干预效果；核对物种、SnoN/SnoN2识别、条带和阴性对照'],['状态问题','若模型存在状态响应，可预设成熟Treg的FOXP3维持等主要结局；配合活细胞数和存活'],['功能支持','有稳定状态变化后，按问题选择等量活Treg功能试验；FOXP3比例不等同抑制能力'],['损伤联系','根据可获得细胞及检测条件观察膜电位、存活和焦亡相关执行层，判断是否与状态变化伴随'],['补偿分支','SKI可作为观察指标；仅在单基因干预充分且有补偿依据时，考虑双基因设计']],[110,395])
p('独立设计的扩展问题','h2');p('初始CD4诱导Treg分化与成熟Treg稳定性分别设定起点、结局和解释。诱导体系可选择FOXP3阳性产率与活细胞数；成熟体系关注状态维持。TGF-β/SMAD近端读数须依据实际问题配置，不由SKIL RNA上调推断通路增强。')
p('推进或调整','h2');p('若巨噬刺激没有可测状态变化，先优化或重新评估模型，暂不投入复杂单/双基因干预。RNA改变而蛋白未改变时先解决检测与调控层解释；仅细胞数量下降需排除死亡选择。只有双干预有效时，定位为联合调控或补偿，不能全部归于SKIL。')
page();p('检测与推进｜分别判定线粒体变化和焦亡','h')
p('有限初筛与机制验证分层','h2');p('模型阶段核对目标T细胞相关gasdermin蛋白、裂解性死亡读数和检测条件。按表达与结果选择GSDMD/GSDME等执行轴，不默认采用巨噬NLRP3—caspase-1—GSDMD流程。[12–13] 焦亡机制结合目标T细胞执行层、膜损伤/死亡以及后续执行轴干预；Annexin V、LDH、IL-1β或caspase-3单项均不足以独立判定。[8]')
p('线粒体初筛结合膜电位、存活和细胞数量；稳定变化后再选独立功能指标。时间与干预用于判断线粒体变化和焦亡的先后关系。对照验证检测方法有效，不证明巨噬刺激已引发同一机制。')
p('主要终点、模型及独立重复','h2');p('在确定候选与体系后预设主要终点、时点、交互尺度和分析模型，区分平行探索。按连续/比例/计数结局及配对来源选择模型，纳入巨噬条件、候选干预及交互；重复时间点按重复测量处理。单个预设主要对比不需FDR；多终点联合判断按预定义同类家族校正。报告效应、区间、P、FDR及独立n。')
tab([['观察','建议下一步'],['模型/蛋白不可测，或干预不足/毒性大','优化模型和检测；记录技术限制，不判断候选无功能'],['交互不显著但干预主效应清楚','按一般T细胞作用解释；不宣称巨噬刺激特异性'],['模型明确、干预充分，预设终点无支持','保留阴性并评估研究问题或候选；不事后更换主要终点'],['候选效应跨独立来源可重复','按研究深度选独立干预、适当回补与上游阻断；不要求三候选全部同时推进'],['需动物/患者、双基因、互作或深度代谢','结合首轮结果、周期和预算另定范围']],[215,290])
p('实施前由实验室确认细胞供给、分选/递送能力、抗体用途及各读数细胞需求，依据先导方差与目标效应计算样本量。剂量、时长与技术条件在选定模型确定；不将本建议直接视为已定稿操作规程。','small')
page();p('附录C｜图表与补充材料索引','h')
p('配套补充材料按图号提供图值表、完整统计及来源说明。图号沿用各分析模块的固定标识；图21—22为筛选细节，位于附录A。')
tab([['图号','内容','配套路径']]+[['图'+str(int(m['id'][1:])),m['title'],'plot_data/'+m['id']+'.tsv'] for m in manifest[:12]],[45,320,140])
page();p('附录C｜图表索引（续）','h')
tab([['图号','内容','配套路径']]+[['图'+str(int(m['id'][1:])),m['title'],'plot_data/'+m['id']+'.tsv'] for m in manifest[12:]],[45,320,140])
p('统计表、图值来源及方法','h2');p('补充材料README.txt列明各目录用途；file_inventory.json提供文件哈希。tables包含完整关联、程序效应、差异、GO/KEGG/GSEA及候选检测；evidence包含背景指标定义、H2-DM来源核对及三候选连接来源。原始P、FDR、模型与检验家族分别保留。','small')
# Layout uses landscape pages for wide scientific panels and portrait for long text/tall panels.
pdf=D/'矽肺免疫全景_CD4候选与验证_报告_v27.pdf'
def footer(c,d):
 w,h=c._pagesize;c.setStrokeColor(colors.HexColor('#c5d8de'));c.line(42,34,w-42,34);c.setFont('CN',8);c.setFillColor(colors.HexColor('#59707a'));c.drawString(42,22,'矽肺免疫微环境与CD4候选研究');c.drawRightString(w-42,22,str(d.page))
class Doc(BaseDocTemplate):
 def afterFlowable(self,f):
  if isinstance(f,Paragraph) and f.style.name in ['h','h2']:figpages.append({'page':self.page,'heading':f.getPlainText()})
doc=Doc(str(pdf),pagesize=(595.28,841.89),leftMargin=45,rightMargin=45,topMargin=40,bottomMargin=45)
doc.addPageTemplates([PageTemplate(id='portrait',pagesize=(595.28,841.89),frames=[Frame(45,45,505.28,756.89,leftPadding=0,rightPadding=0,topPadding=0,bottomPadding=0)],onPage=footer),PageTemplate(id='landscape',pagesize=(841.89,595.28),frames=[Frame(45,45,751.89,510.28,leftPadding=0,rightPadding=0,topPadding=0,bottomPadding=0)],onPage=footer)])
doc.build(story)
(O/'report_text.txt').write_text('\n\n'.join(textrecords),encoding='utf-8');(O/'report_page_index.json').write_text(json.dumps(figpages,ensure_ascii=False,indent=2),encoding='utf-8');(O/'report_delivery.json').write_text(json.dumps({'pdf':str(pdf),'sha256':hashlib.sha256(pdf.read_bytes()).hexdigest(),'figures':len(manifest),'explanatory_displays':3,'total_visuals':len(manifest)+3},ensure_ascii=False,indent=2),encoding='utf-8');print('REPORT',pdf)