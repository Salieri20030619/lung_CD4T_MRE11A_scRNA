from pathlib import Path
import json,html,hashlib
import pandas as pd
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.platypus import BaseDocTemplate,PageTemplate,Frame,Paragraph,Spacer,PageBreak,NextPageTemplate,Table,TableStyle,Image,KeepTogether
from PIL import Image as PILImage
O=Path(__file__).parent;H=O.parent;D=H/'Downloads/矽肺_CD4候选研究_正式报告_v21_20261017';D.mkdir(exist_ok=True)
manifest=json.loads((O/'figure_manifest.json').read_text(encoding='utf-8'));assert len(manifest)==22
pdfmetrics.registerFont(TTFont('CN',r'C:\Windows\Fonts\ARIALUNI.TTF'))
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
p('矽肺免疫微环境与\nCD4⁺ T细胞候选研究','h');p('单细胞图谱、功能程序、细胞通讯及实验候选整合报告','h2');p('版本 v20｜2026年10月17日｜客户讨论与基础实验设计用')
p('主要结论','h2');p('主小鼠队列的全肺图谱定位了主要免疫群，核心CD4细胞呈Naive、ActCD4与Treg三类状态。疾病样本的Naive占比降低、ActCD4及Treg占比升高；在3+3个独立样本的精确比较下，这些组成变化作为探索性方向呈现。')
p('疾病差异与功能文献共同支持优先验证IL1R2/Treg；TTC39C/ActCD4提供较强的凋亡相关转录程序线索，作为损伤方向探索备选；SKIL/Treg适合作为状态调控备选。三个候选的首轮实验问题不同，采用分方向推进。')
p('当前最明确的是疾病相关状态差异及部分程序关联。主队列焦亡/凋亡评分的组间比较、部分疾病调整和补充队列结果未提供一致支持，尚不能将任一候选写成已确认的“巨噬细胞诱导CD4线粒体损伤或焦亡”机制。')
tab([['建议顺序','优先细胞/终点','核心理由'],['IL1R2 主选','Treg功能、代谢反馈及损伤耐受','疾病DE较明确；已有Treg/Tfr功能干预研究'],['TTC39C 探索备选','ActCD4存活、凋亡及线粒体状态','ApoExec关联较强；干预作用方向需先确认'],['SKIL 状态备选','Treg分化、稳定性及抑制功能','TAct关联与TGF-β/Treg研究；关注SnoN蛋白、SKI冗余']],[105,175,225])
p('报告含11个图示模块、22张代表图。各图给出样本范围、统计单位与图意说明；完整结果保留原始P及各自家族FDR。')
page();p('研究设计与数据范围','h')
tab([['数据集','用途','分析单位与边界'],['GSE314943','主发现队列；全景、CD4状态、巨噬入口','3个Saline+3个Silica样本；37,155个全肺细胞；Core CD4共2,311个'],['GSE183682','疾病方向补充及联合DE/关联','D7、D56各1个Saline+1个Silica；单时间点仅描述方向'],['联合分析','增强疾病效应估计并控制队列/时间层','共10个样本，包含主队列6个；不能算另一独立复现'],['GSE274372','真实时间背景与程序联系补充','D0/D7/D28/D56各3个样本；TLike共494细胞；min10敏感性为11样本']],[105,180,220])
p('统计与筛选原则','h2');p('样本级推断以独立生物学样本/文库为单位，细胞、轨迹区间和技术重复不替代样本量。样本少于10的结果按探索性证据呈现。原始目标效应P<0.05用于名义信号筛选，FDR保留在预定义同类假设家族内，不作为统一入池硬门槛。')
p('候选池由55个原多条件交集基因与42个历史基因合并，重叠7个，得到90个复评基因。当前A1=52、A2=37、A3=1。52表示计算支持入池，并非52个跨队列一致差异基因。重点9基因/10状态组合是发现后的重点复评，三个预实验建议也不是全部候选创新性的无偏排序。')
p('图示阅读顺序','h2');p('图1-2免疫全景；图3-4 CD4亚群；图5-6差异表达；图7-8 GO/KEGG；图9-10 GSEA；图11-12细胞通讯；图13-14转录因子网络；图15-16焦亡/凋亡评分；图17-18状态轨迹；图19-20 NicheNet；图21-22候选与验证。')
page();p('分析方法','h')
p('图谱、亚群与差异表达','h2');p('复用主队列冻结质控、注释和UMAP/Harmony坐标，按全肺谱系及Core CD4状态展示。细胞比例以每样本Core CD4为分母，用6样本全部20种组标签分配进行双侧精确检验，3个状态构成一个BH家族。差异表达采用原始计数样本聚合、完整文库分母、表达过滤、TMM归一化及edgeR准似然模型；联合模型包含队列/时间层和疾病效应。所有DE及原完整家族FDR沿用已复核结果。')
p('GO/KEGG与GSEA','h2');p('为与当前筛选一致，使用冻结联合DE中的ActCD4及Treg完整结果补充功能分析。GO BP来自本次环境org.Mm.eg.db的GOALL注释，KEGG复用项目冻结小鼠映射。基因符号按既有规范统一；以可映射受检基因为ORA背景。原始P<0.05的上、下调基因分别作超几何检验，集合限制10-500个背景成员，零命中集合保留在完整校正家族中。')
p('GSEA对全部受检基因按sign(logFC)×sqrt(QL F)排序，同分按基因名固定顺序。fgseaMultilevel使用10-500成员、seed=20261017、nPermSimple=1000、eps=1e-8。其P为排序基因集合富集模型统计，非新增独立样本处理效应检验。每状态/数据库分别BH校正；13条结果达到P估计下限，按P≤1e-8解释。代表图选择规则固定在图注及绘图数据中。')
p('通讯与转录因子网络','h2');p('通讯以冻结数据库相互作用、发送端和接收端组级覆盖为依据，按必要受体组分评估表达相容性。图11展示合格相互作用数量；图12展示样本级发送配体表达。PROS1相关功能支持已撤回。TF部分复用SCENIC网络及v18去候选/程序重叠后的regulon，按独立样本平均AUCell；代表网络按保留边权重每TF取前2个靶标。')
page();p('评分、关联与状态轨迹方法','h')
p('样本程序评分','h2');p('对样本s与状态k，C_gsk为该样本该状态细胞的基因原始计数和，L_sk为完整原始文库UMI之和。E_gsk=10000×C_gsk/L_sk；候选表达G_gsk=ln(1+E_gsk)，模块评分M_msk为保留映射成员上mean_g ln(1+E_gsk)。先样本汇总、后归一化和取对数。零表达保留；缺失映射不填零；至少3成员且映射覆盖≥50%才评价。9候选从模块移除，实际重叠者DDIT4从MitoISR/Glycolysis移除。')
p('Core等权状态标准化敏感性给Naive/ActCD4/Treg各1/3权重：候选先平均各状态E再log，模块和AUCell平均各状态分值。它与旧版六样本平均比例权重（0.466639/0.309979/0.223382）分开。历史轨迹使用逐细胞匹配对照评分后聚合，不与当前汇总计数评分直接比较数值。')
p('样本关联及精确置换','h2');p('联合秩关联在队列/时间层设计下残差化，块内排列6!×2!×2!=2880。疾病调整模型包含疾病项，组内有效排列3!×3!=36，第二队列4个单例单元不能贡献组内交换。联合组间比较保持层内组别数，共80种分配。真实时间关联控制日龄，4组各3样本时有1296种组内排列，去掉1个低细胞样本后432种。每类模型保留自身完整FDR家族。')
p('状态拟时序与NicheNet','h2');p('状态轨迹复用主队列Naive/ActCD4的4套冻结DPT排序（K15/K30及两种根），不重新选择根或重建Treg分支。每状态内将DPT缩放，并在每个样本中调整log文库量后拟合表达斜率；以6个样本斜率作64种符号翻转检验。图中4个分位区间只用于展示。')
p('NicheNet复用冻结配体-靶先验、表达受体资格、差异目标集和背景，展示联合疾病上调目标集的配体活动及当前候选路线。活动Pearson及网络权重作为模型排序指标，保留各状态实际目标集和背景。')
# Results lead paragraphs, one figure per page with explicit source and method caption.
results={
'F01':'全肺注释识别出B、T、巨噬/单核、中性粒、树突等免疫群，并同时保留上皮、间质及血管细胞作为组织环境。后续巨噬-CD4分析均在这一细胞来源框架内展开。',
'F02':'T细胞相关CD3D/TRAC、B细胞相关MS4A1/CD79A及髓系相关标志呈现不同分布，支持主要谱系区分。图中广泛出现的组织相关转录背景也说明候选定位需要结合亚群身份及背景敏感性。',
'F03':'Core CD4由933个Naive、801个ActCD4和577个Treg构成。疾病组更多分布于ActCD4/Treg状态区域，样本层面的组成量化见下一图。',
'F04':'疾病组Naive平均比例降低46.85个百分点，ActCD4增加20.85个百分点，Treg增加26.00个百分点；对应精确P为0.10、0.20、0.10，FDR为0.15、0.20、0.15。方向清楚，但当前样本量下未达到名义显著。',
'F05':'Treg中的IL1R2与SKIL在联合疾病模型呈上调信号。TTC39C的联合疾病DE未达到名义显著，因此后续推荐它主要依据凋亡程序联系，而非稳定跨队列差异。',
'F06':'IL1R2/Treg联合logFC=3.290、P=1.66×10^-6、FDR=0.000662；SKIL/Treg联合logFC=1.438、P=2.15×10^-5、FDR=0.00314。FOXP1/Treg为下调状态候选，继续保留在证据池。',
'F07':'Treg的下调集合突出细胞质翻译/核糖体相关过程，同时上调集合包含白细胞活化；ActCD4的部分免疫相关过程获得富集支持。具体含义按命中成员解释，涉及免疫球蛋白或广泛背景的条目另表标注。',
'F08':'Treg下调差异基因在核糖体相关KEGG集合中富集，并出现上调细胞因子-受体相互作用线索。ActCD4细胞黏附相关集合原P=0.000284、FDR=0.0804，为名义信号。',
'F09':'全基因排序与阈值差异分析共同显示翻译/核糖体相关集合偏向疾病下调端。部分免疫相关集合偏向上调端，GO/KEGG同名或共享成员的结果不计作独立验证。',
'F10':'代表核糖体通路在Treg和ActCD4均为负NES。此处展示完整排名中的成员累积位置，可与ORA方向互相参照；它反映转录集合变化，不直接测量蛋白合成速度。',
'F11':'在冻结覆盖规则下，巨噬来源至CD4各状态均存在表达相容入口，疾病组满足覆盖的相互作用数量较多。此结果用于确定可检验信号范围，不能把计数变化解释成实际通讯强度增加。',
'F12':'AMlike中的IL1A、IL1B、CCL2和TNF按独立样本呈现不同表达分布。后续上游验证应保持同一发送范围，结合完整受体和近端响应选择具体路线。',
'F13':'去重叠regulon的AUCell呈现CD4状态与样本间差异，为候选相关程序提供调控环境。优先候选RNA、TF靶网络及活性评分是互补证据，不能相互替代。',
'F14':'FOS/JUN、NFAT、NF-κB、FOXP1与TRP53等冻结调控网络可连接具体靶成员。图中代表边用于展示网络构成；完整保留成员用于后续复核和靶标选择。',
'F15':'ActCD4与Treg的PyroExec组间效应方向为上升，但主队列精确P分别为0.30与0.20；Inflam对应P均为0.90。当前组间检验不支持将疾病相关CD4状态变化直接等同于已发生焦亡。',
'F16':'ActCD4的ApoExec差值为+0.0404、P=0.10、FDR=0.274；Treg差值+0.0352、P=0.20、FDR=0.433。基因-程序相关与组间评分差异是不同问题，前者的阳性不替代后者。',
'F17':'四套冻结DPT的细胞排序相关为0.904-0.999，支持当前Naive/ActCD4状态排序对所测设置较稳定。其适用范围为转录状态连续性，Treg不在这条轨迹内。',
'F18':'主要排序中，Naive的PyroExec斜率P=0.0625、ApoExec P=0.125；ActCD4对应P=0.156与0.188。当前数据不支持用这些程序曲线证明沿轨迹发生持续焦亡/凋亡增强。',
'F19':'TNF、TGFB1及若干免疫调节配体在联合疾病上调目标集合中获得模型排序线索。配体优先级由目标集、背景及受体条件共同决定，应在具体发送细胞范围中继续定位。',
'F20':'IL1R2和SKIL在当前联合名义差异及网络条件下保留代表连接。TTC39C在该图的当前联合DE路线中为空缺；其原队列/历史先验与程序关联另行保留，不补造上游链条。',
'F21':'IL1R2/Treg阳性细胞145个，均共检至少一个CD3D/CD3E/TRAC；TTC39C/ActCD4为125/125；SKIL/Treg为319/331。背景≤1%保留的阳性细胞分别为141、119和300。候选表达方向在Raw/BG1及样本留一描述中一致。',
'F22':'IL1R2与Inhib、TTC39C与ApoExec、SKIL与TAct均有联合名义关联。推荐以各自最直接的细胞状态/功能问题开展预实验，并保留疾病调整及外部方向结果决定下一步投入。'}
for m in manifest:
 template='portrait' if m['id'] in ['F15','F16'] else 'landscape';page(template);p(m['module']+'｜'+m['title'],'h2');p(results[m['id']],'p')
 im=PILImage.open(m['file']);maxw=505 if template=='portrait' else 752;maxh=465 if template=='portrait' else 315;scale=min(maxw/im.width,maxh/im.height);story.append(Image(m['file'],width=im.width*scale,height=im.height*scale,hAlign='CENTER'));story.append(Spacer(1,6));p('图'+str(int(m['id'][1:]))+'。'+m['caption'],'small');p('数据索引：'+m['id']+'；图值、来源表及脚本见配套材料索引。','small')
page();p('候选比较与外部方向','h')
tab([['候选','当前项目证据','关键未决环节'],['IL1R2 / Treg','原六文库logFC=3.115，P=4.89×10^-5，FDR=0.00853；联合DE及Inhib支持；ENO1关联P=0.000694、FDR=0.0285','ENO1疾病调整P=0.111；PyroExec P=0.988；不能用RNA正相关证明酶活抑制'],['TTC39C / ActCD4','ApoExec联合rho=0.943、P=0.00347、FDR=0.0556；疾病调整rho=0.922、P=0.0278、FDR=0.556','联合DE P=0.163、FDR=0.803；D7方向受归一化影响，D56反向；上游未补齐'],['SKIL / Treg','TAct联合rho=0.965、P=0.00208、FDR=0.0744；疾病调整rho=0.897、P=0.0278','ApoExec P=0.0747、调整P=0.917；p53/凋亡链未补齐；AMlike CCL2不能跨来源拼接']],[90,210,205])
p('真实时间队列提供不同问题的补充证据','h2');p('GSE274372的TLike时间分析中，TTC39C与MitoISR在日龄调整后rho=0.733、P=0.0185、FDR=0.567，为探索性信号。它不是ActCD4-ApoExec同一关联的独立复现。时间、疾病、细胞状态和损伤终点分别评价。')
p('其余候选继续按证据分层保留','h2');p('FOXP1优先解释状态重塑，Core-PyroExec关联在等权状态标准化后P=0.494；MDFIC为损伤关联备选，调整后rho=0.950但P=0.0833；DDIT4保留文献驱动假说，项目内PyroExec P=0.751；GSN保留备选，PROS1功能入口支持已撤回。AHNAK/ST6GALNAC3保留氧化代谢或应激转录线索。')
page();p('功能文献为实验提供依据与方向约束','h')
p('IL1R2：以Treg功能为首轮切入点 [1-4]','h2');p('2022年Treg特异Il1r2缺失研究支持其参与肿瘤Treg克隆扩增和效应特征；2024年Tfr研究支持其调节特定免疫环境中的IL-1作用。2025年IL1R2-ENO1-糖酵解/焦亡链来自巨噬细胞及脓毒症。矽肺相关研究也已有其他细胞对象的Il1r2表达报道。研究空间应聚焦矽肺巨噬环境对CD4/Treg及其功能的影响，允许识别保护性反馈。')
p('TTC39C：优先确认干预后的存活与凋亡方向 [5]','h2');p('人A549/NCI-H1299细胞的shRNA实验中，TTC39C降低后Annexin V阳性增加。项目中TTC39C与ApoExec正相关不能直接改写成促进凋亡；首轮要判断在CD4中是驱动、补偿还是其他背景依赖作用。')
p('SKIL：同时观察SnoN蛋白与SKI [6-7]','h2');p('SKIL编码SnoN，SKI为另一基因。小鼠Treg分化研究中，在Arkadia缺失背景下，双靶向Ski/SnoN可救援，单独靶向不能完成该救援，提示特定条件下的冗余。SKIL RNA上调不等同SnoN蛋白或功能增强。PM2.5/ROS/SnoN与肺纤维化已有研究联系，广义关联本身不作为全新机制。')
p('检测与干预可行性核查','h2');p('IL1R2具有小鼠Tfr流式4E2检测先例及厂家小鼠适用性资料。TTC39C有厂家小鼠蛋白WB资源，但原代CD4中的检测和干预效力仍需确认。SnoN原T细胞研究有抗体验证；当前替代试剂需核对物种、SnoN/SnoN2识别和阴性对照，不能把仅人源验证的产品直接视为小鼠适用。具体资源与出处保留在实验可行性表，不将“可采购”写成“已在目标模型有效”。')
page();p('基础实验建议与投入判断','h')
p('第一阶段：确认模型、亚群与蛋白','h2');p('比较对照与矽尘处理巨噬条件对目标CD4亚群的影响。IL1R2/SKIL优先Treg，TTC39C优先ActCD4；检测RNA及对应蛋白，同时观察细胞数量、功能及存活。设置颗粒直接作用/残留相关对照，依据假说选择条件培养基、隔离或接触共培养，区分可溶性、接触依赖及颗粒作用。')
p('第二阶段：检验候选是否参与巨噬刺激效应','h2');tab([['巨噬条件','非靶向对照','候选干预'],['对照巨噬','基础状态','干预本身效应/毒性'],['矽尘处理巨噬','刺激效应','干预是否改变刺激响应']],[120,190,195])
p('以“巨噬条件×候选干预”的交互效应为预设主要比较，同时记录基础毒性和干预效率。独立培养/独立动物来源作为重复，技术孔不替代生物学重复。预实验用于估计变异和效应，再确定正式样本量；不依据单次阳性或更换阈值决定投入。')
p('IL1R2优先Treg抑制功能及损伤耐受；TTC39C优先存活、凋亡执行层与线粒体状态；SKIL优先Treg分化/稳定性及抑制功能，同时检测SnoN与SKI。若主张焦亡，需增加相应gasdermin裂解及膜完整性等执行层证据，一般死亡或炎症因子变化不能替代。')
p('第三阶段：对可重复有效候选做救援与上游定位','h2');p('先用独立干预及适当回补/救援确认候选相关性，再从冻结上游范围选择同来源、完整受体及近端响应相容的信号进行阻断或补充。模型无稳定表型时先完善模型；干预充分仍无目标表型变化时，调整候选或终点。')
page();p('讨论：证据边界与研究空间','h')
p('项目支持候选选择，不预先决定湿实验方向','h2');p('疾病差异、程序关联、网络先验和干预功能是不同证据层级。IL1R2具备较清楚的疾病差异和Treg功能研究基础；TTC39C提供损伤相关程序切入点；SKIL具备Treg状态和调控背景。三者均有预实验价值，但不对应已量化的成功概率。')
p('样本与模型边界','h2');p('主队列为3+3个样本；第二队列各时间点1+1，联合分析含原队列；时间队列为TLike范围。精确枚举保证既定排列空间计算完整，不证明可交换性或模型假设成立。MDFIC等强效应而P有限的结果不能全部归因于分辨率。')
p('生物学与技术解释','h2');p('转录评分不能替代线粒体膜电位、呼吸、蛋白裂解或死亡实测。富集结果中的Coronavirus disease、synapse等名称可能由共享核糖体/翻译成员驱动，不表示感染或神经细胞功能；免疫球蛋白/组织背景成员需结合背景敏感性解释。原始细胞图谱的未确定标签保持原状，不因候选需要而重定义。')
p('PROS1-AXL的胶质瘤文献背景不能直接用于小鼠CD4。TYRO3/MERTK在现有覆盖规则下未补回其功能入口；RNA低检出不是蛋白不存在的证明。其他网络边也仅保留计算支持，不视为已经验证。状态DPT是横断面转录排序，不能代替谱系追踪或真实损伤时间。')
p('创新性落在具体细胞关系和作用环节','h2');p('IL1R2-ENO1-焦亡、SKIL/SnoN-TGF-β-Treg及TTC39C-肿瘤细胞存活均已有研究。新的问题应落在矽肺巨噬环境如何改变CD4/Treg候选及功能，以及候选是否参与损伤耐受或细胞死亡。定向检索未取得完整链直接功能证明，不等同全球首次声明。')
page();p('参考研究与可复核材料','h')
refs=[('[1]','IL1R2 increases regulatory T cell population in the tumor microenvironment by enhancing MHC-II expression on cancer-associated fibroblasts. JITC, 2022.','https://pmc.ncbi.nlm.nih.gov/articles/PMC9438093/'),('[2]','Loss of T follicular regulatory cell-derived IL-1R2 augments germinal center reactions via increased IL-1. JCI Insight, 2024.','https://insight.jci.org/articles/view/174005'),('[3]','Critical Role of IL1R2-ENO1 Interaction in Inhibiting Glycolysis-Mediated Pyroptosis for Protection Against Lethal Sepsis. Advanced Science, 2025.','https://pmc.ncbi.nlm.nih.gov/articles/PMC12533375/'),('[4]','Epithelial NLRP3 drives silica-induced lung injury and fibrosis through IL-18 and pro-fibrotic neutrophil recruitment. Particle and Fibre Toxicology, 2026.','https://pmc.ncbi.nlm.nih.gov/articles/PMC13214385/'),('[5]','Ttc39c is a potential target for the treatment of lung cancer. BMC Pulmonary Medicine, 2022;22:391.','https://pmc.ncbi.nlm.nih.gov/articles/PMC9615393/'),('[6]','Arkadia-SKI/SnoN signaling differentially regulates TGF-β-induced iTreg and Th17 cell differentiation. JEM, 2021.','https://pmc.ncbi.nlm.nih.gov/articles/PMC8421263/'),('[7]','PM2.5 regulates TGF-β1/Smads-mediated pulmonary fibrosis via ROS/SnoN in vitro and in vivo. Int Immunopharmacol, 2025.','https://pubmed.ncbi.nlm.nih.gov/40154177/')]
for n,title,url in refs:p(n+' '+title+'\n'+url,'small')
p('数据与图件复核范围','h2');p('v19独立复核报告确认其实际检查范围：2,852个评分值、2,544个相关系数、2,268个穷举置换P、34个设计、221项程序效应及193个原BH家族。276个Monte Carlo P、原始细胞层edgeR重拟合、全套注释/SCENIC等不在该独立复算范围。','small')
p('本次新增通路、比例及代表图由分析端另行核查，不能扩大为全部分析已独立重跑。完整统计、图值、方法、文献/试剂依据及校验清单保存在v20配套材料。原v19与更早版本保留，候选文献优先于旧实验排序表述。报告生成不代表已向外部发送。','small')
# Layout uses landscape pages for wide scientific panels and portrait for long text/tall panels.
pdf=D/'矽肺免疫全景_CD4候选与验证_正式报告_v21.pdf'
def footer(c,d):
 w,h=c._pagesize;c.setStrokeColor(colors.HexColor('#c5d8de'));c.line(42,34,w-42,34);c.setFont('CN',8);c.setFillColor(colors.HexColor('#59707a'));c.drawString(42,22,'矽肺免疫微环境与CD4候选研究');c.drawRightString(w-42,22,str(d.page))
class Doc(BaseDocTemplate):
 def afterFlowable(self,f):
  if isinstance(f,Paragraph) and f.style.name in ['h','h2']:figpages.append({'page':self.page,'heading':f.getPlainText()})
doc=Doc(str(pdf),pagesize=(595.28,841.89),leftMargin=45,rightMargin=45,topMargin=40,bottomMargin=45)
doc.addPageTemplates([PageTemplate(id='portrait',pagesize=(595.28,841.89),frames=[Frame(45,45,505.28,756.89,leftPadding=0,rightPadding=0,topPadding=0,bottomPadding=0)],onPage=footer),PageTemplate(id='landscape',pagesize=(841.89,595.28),frames=[Frame(45,45,751.89,510.28,leftPadding=0,rightPadding=0,topPadding=0,bottomPadding=0)],onPage=footer)])
doc.build(story)
(O/'report_text.txt').write_text('\n\n'.join(textrecords),encoding='utf-8');(O/'report_page_index.json').write_text(json.dumps(figpages,ensure_ascii=False,indent=2),encoding='utf-8');(O/'report_delivery.json').write_text(json.dumps({'pdf':str(pdf),'sha256':hashlib.sha256(pdf.read_bytes()).hexdigest(),'figures':len(manifest)},ensure_ascii=False,indent=2),encoding='utf-8');print('REPORT',pdf)