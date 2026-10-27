from pathlib import Path
import json,html,hashlib
import pandas as pd
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.platypus import BaseDocTemplate,PageTemplate,Frame,Paragraph,Spacer,PageBreak,NextPageTemplate,Table,TableStyle,Image,KeepTogether
from PIL import Image as PILImage
O=Path(__file__).parent;H=O.parent;D=H/'Downloads/矽肺_CD4候选研究_实验方案修订_v26_20261017';D.mkdir(exist_ok=True)
manifest=json.loads((O/'figure_manifest.json').read_text(encoding='utf-8'));assert len(manifest)==24
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

def overview_image(id,maxh=520):
 fn=H/'lung_cohort_bridge_20261017_v25/figures'/(id+'.png');im=PILImage.open(fn);k=min(505/im.width,maxh/im.height);story.append(Image(str(fn),width=im.width*k,height=im.height*k,hAlign='CENTER'));story.append(Spacer(1,7))
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
 page();p('差异分析前的筛选框架：三项证据交集与后续分层','h2')
 p('差异分析提供候选变化证据，随后结合目标损伤/死亡功能注释及同状态巨噬入口。三项证据交集得到55个基因，再与补充候选形成90个基因的评价集合。创新性与实验可行性参与后续分层，不作为第四个统一硬门槛。')
 overview_image('S_C')
 p('示意图C。圆形交叠仅用于说明筛选逻辑，面积不代表集合大小；55为按细胞状态和对比匹配后去重的三项交集基因数，不完整列示各互斥区域数量。功能维度来自GO/Reactome注释及直接功能文献，不以本项目GO富集P值替代单基因功能依据。逐步阈值、集合合并和候选去向详见后文筛选章节。','small')

p('矽肺免疫微环境与\nCD4⁺ T细胞候选研究','h');p('单细胞图谱、功能程序、细胞通讯及实验候选整合报告','h2');p('2026年10月17日｜研究报告')
p('研究目的','h2');p('本研究围绕“矽肺巨噬细胞环境与CD4⁺ T细胞状态及损伤反应有何联系”展开，利用公开单细胞数据筛选可进入基础实验的候选基因，并确定优先验证的细胞亚群和功能终点。')
p('数据与分析流程','h2');p('以GSE314943的3个对照和3个矽尘样本为主队列，构建37,155个全肺细胞的图谱，重点分析其中2,311个CD4细胞。另纳入GSE183682的4个样本开展分层联合分析，并以GSE274372在4个时间点的12个样本提供时间背景；联合模型包含主队列，不作为另一次独立验证。')
p('分析依次包括免疫谱系与CD4亚群组成、样本级差异表达、GO/KEGG及GSEA、巨噬细胞至CD4的配体-受体候选、转录因子网络、焦亡/凋亡及线粒体相关程序评分、状态拟时序和NicheNet配体-靶基因预测。结合功能文献，对90个基因进行综合评价，并重点分析9个基因的10个细胞状态组合。')
p('主要发现与推荐方向','h2');p('疾病组呈Naive比例下降、ActCD4和Treg比例上升的方向；对应组成检验未达到P<0.05。综合筛选得到52个符合计算条件的基因，但各候选的证据侧重不同。结合疾病差异、程序关联和既往功能研究，建议优先验证以下3个方向。')
tab([['候选与定位','直接依据','首轮实验终点'],['IL1R2 / Treg\n优先验证','联合疾病上调：logFC=3.290，P=1.66×10^-6，FDR=0.000662；有Treg功能干预研究','Treg功能及损伤耐受，按结果再评价代谢机制'],['TTC39C / ActCD4\n损伤方向探索备选','ApoExec名义关联：rho=0.943，P=0.00347，FDR=0.0556；联合疾病差异未达名义显著','存活、凋亡与线粒体状态；先确认干预作用方向'],['SKIL / Treg\n状态调控备选','TAct名义关联：rho=0.965，P=0.00208，FDR=0.0744；有Treg分化研究','Treg状态与功能，同时检测SnoN蛋白及SKI补偿']],[105,225,175])
p('现有结果支持候选选择及预实验设计，尚未证明“巨噬信号→候选基因→CD4线粒体损伤或焦亡”的完整因果链。下一步用巨噬刺激与候选干预的组合实验，检验候选是否实际改变T细胞表型。','small')
page();p('研究设计与数据范围','h')
tab([['数据集','用途','分析单位与边界'],['GSE314943','主发现队列；全景、CD4状态、巨噬入口','3个Saline+3个Silica样本；37,155个全肺细胞；Core CD4共2,311个'],['GSE183682','疾病方向补充及联合DE/关联','D7、D56各1个Saline+1个Silica；单时间点仅描述方向'],['分层联合疾病模型','GSE314943的6个样本＋GSE183682的4个样本','5个对照＋5个矽尘；控制队列-时间层，包含主队列，不能算另一独立复现'],['GSE274372','真实时间背景与程序联系补充','D0/D7/D28/D56各3个样本；TLike共494细胞；min10敏感性为11样本']],[105,180,220])
p('统计与筛选原则','h2');p('样本级推断以独立生物学样本/文库为单位，细胞、轨迹区间和技术重复不替代样本量。样本少于10的结果按探索性证据呈现。原始目标效应P<0.05用于名义信号筛选，FDR保留在预定义同类假设家族内，不作为统一入池硬门槛。')
p('候选池由55个原多条件交集基因与42个历史基因合并，重叠7个，得到90个复评基因。当前A1=52、A2=37、A3=1。52表示计算支持入池，并非52个跨队列一致差异基因。重点9基因/10状态组合是发现后的重点复评，三个预实验建议也不是全部候选创新性的无偏排序。')
p('图示阅读顺序','h2');p('示意图A/B说明队列关系与补充队列图谱；示意图C在差异分析前介绍交集框架。图1-2免疫全景；图3-4 CD4亚群；图5-6差异表达；图7-8 GO/KEGG；图9-10 GSEA；图11-12细胞通讯；图13-14转录因子网络；图15-16焦亡/凋亡评分；图17-18状态轨迹；图19-20 NicheNet；图21-22候选筛选流程；图23-24候选与验证。')
cohort_pages()
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
'F21':'IL1R2/Treg阳性细胞145个，均共检至少一个CD3D/CD3E/TRAC。\nTTC39C/ActCD4为125/125；SKIL/Treg为319/331。背景≤1%保留的阳性细胞分别为141、119和300。候选表达方向在Raw/BG1及样本留一描述中一致。',
'F22':'IL1R2与Inhib、TTC39C与ApoExec、SKIL与TAct均有联合名义关联。推荐以各自最直接的细胞状态/功能问题开展预实验，并保留疾病调整及外部方向结果决定下一步投入。'}

results['F09']='Treg的KEGG富集曲线同时呈现疾病上调端与下调端的通路分布。核糖体相关集合偏向下调端；其他通路需结合其具体成员及细胞状态解释。'
results['F10']='ActCD4的多通路曲线显示不同功能集合在完整基因排序中的富集位置。正、负NES分别对应疾病上调端和下调端，图中并列报告每条通路的P值及FDR。'
results['F11']='巨噬细胞与Naive、ActCD4和Treg之间均存在满足所设表达覆盖的配体-受体候选。分组弦图展示这些候选的状态分布，并保留各发送范围的独立定义。'
results['F12']='Treg候选关系包括IL1A/IL1B-IL1R2、CCL2-CCR2、CXCL10-CXCR3、CD200-CD200R1及SEMA3E-PLXND1，可用于安排后续来源、受体与近端响应验证。'


results['F17']='四套DPT排序的相关为0.904-0.999。叠加的邻域主干连接Naive/ActCD4的主要状态区域，帮助阅读DPT方向；线条来自原始邻域连接的描述性抽象，不作为真实分化过程或损伤进展的证明。'
results['F20']='矩阵直接显示IL1R2的诱饵受体关系与SKIL的配体-靶先验连接。TTC39C在本展示对应的联合差异及网络条件下没有合格连接，其程序关联证据仍单独评价。'


results['F23']=results['F21']
results['F24']=results['F22']
results['F21']='主要差异模型共涉及9,091个去重基因；2,765个在至少一个目标对比中达到原始P<0.05。与目标功能集合相交后保留111个，再匹配相同接收状态的巨噬入口得到55个。结合补充候选形成90基因评价集合。'
results['F22']='综合评价最终支持52个计算候选，37个保留为证据待补，1个未达到当前名义差异条件。重点9基因不是52基因的继续硬筛子集；TTC39C虽然处于功能证据待补层，仍可凭凋亡程序关联及功能实验先例进入探索性预实验。'
def screening_methods():
 page();p('候选筛选：证据维度、阈值与集合关系','h')
 p('筛选首先整合差异、目标损伤/死亡功能和巨噬来源入口三项证据，再进行细胞身份与作用角色评价。创新性及实验可行性用于确定研究方向，不要求所有候选在四项指标上同时通过一个统一阈值。')
 tab([['维度','实际判定条件','参考证据与解释'],['差异表达','有效样本级模型的目标效应原始P<0.05；按基因-状态-对比保留记录','FDR、|log2FC|≥0.5、跨队列方向及敏感性用于评价强度；不要求所有队列同时显著'],['目标损伤/死亡功能','GO/Reactome内源性凋亡、线粒体损伤相关过程或焦亡集合；明确的直接功能文献可扩展','两类机制取并集；不要求基因同时参与凋亡和焦亡。一般凋亡/状态证据不足以证明目标细胞自身损伤'],['巨噬来源入口','配体和完整受体组分在同一条件及接收状态内满足覆盖；候选为受体组分或具有对应配体-靶先验','每样本至少10个目标细胞，阳性比例≥5%，且组内至少一半合格样本满足。网络先验不等于功能激活'],['创新性与可行性','逐项比较已有机制、细胞背景、研究空间与可观察实验表型','作为分层依据。已有发表不自动排除，未检索到文献不自动证明创新；不产生二元创新性P值']],[80,215,210])
 p('交集按同一细胞状态匹配；时间证据单列','h2')
 p('主要DE范围包括主队列4个细胞范围、联合队列4个范围及时间队列3个时间对比和1个总体时间检验，共12个家族。名义差异在这些目标对比间取并集，时间效应与疾病效应分别解释；功能和网络匹配完成后再按基因去重，避免跨状态拼接。','small')
 p('初筛与综合评价的网络范围','h2')
 p('初筛网络分支采用校正AUPR>0的合格配体，按活动排序取前30个，每配体取表达背景内前250个靶基因；直接受体分支单独保留。综合评价时扩展检查有正向活动的合格配体，网络排名作为参考，不因超过30/250而统一排除。','small')
 p('补充集合、资格分层与预实验推荐','h2')
 p('42个补充候选围绕T细胞状态、细胞死亡及相关调控问题纳入，与55个交集基因合并后统一评价；来源及成员逐项列入附表。细胞来源不清或主要作用于其他细胞的输出效应基因进入待补层。重点分析选择IL1R2、DDIT4、SKIL、FOXP1、GSN、AHNAK、TTC39C、MDFIC、ST6GALNAC3；预实验顺序进一步依据功能文献、程序关联和可操作性确定。','small')

for m in manifest:
 if m['id']=='F05':early_framework()
 if m['id']=='F21':screening_methods()
 template='portrait';page(template);p(m['module']+'｜'+m['title'],'h2');p(results[m['id']],'p')
 im=PILImage.open(m['file']);maxw=505;maxh=540;scale=min(maxw/im.width,maxh/im.height);story.append(Image(m['file'],width=im.width*scale,height=im.height*scale,hAlign='CENTER'));story.append(Spacer(1,6));p('图'+str(int(m['id'][1:]))+'。'+m['caption'],'small');p('数据索引：'+m['id']+'；图值、来源表及脚本见配套材料索引。','small')
page();p('候选比较与外部方向','h')
tab([['候选','当前项目证据','关键未决环节'],['IL1R2 / Treg','原六文库logFC=3.115，P=4.89×10^-5，FDR=0.00853；联合DE及Inhib支持；ENO1关联P=0.000694、FDR=0.0285','ENO1疾病调整P=0.111；PyroExec P=0.988；不能用RNA正相关证明酶活抑制'],['TTC39C / ActCD4','ApoExec联合rho=0.943、P=0.00347、FDR=0.0556；疾病调整rho=0.922、P=0.0278、FDR=0.556','联合DE P=0.163、FDR=0.803；D7方向受归一化影响，D56反向；上游未补齐'],['SKIL / Treg','TAct联合rho=0.965、P=0.00208、FDR=0.0744；疾病调整rho=0.897、P=0.0278','ApoExec P=0.0747、调整P=0.917；p53/凋亡链未补齐；AMlike CCL2不能跨来源拼接']],[90,210,205])
p('真实时间队列提供不同问题的补充证据','h2');p('GSE274372的TLike时间分析中，TTC39C与MitoISR在日龄调整后rho=0.733、P=0.0185、FDR=0.567，为探索性信号。它不是ActCD4-ApoExec同一关联的独立复现。时间、疾病、细胞状态和损伤终点分别评价。')
p('其余候选继续按证据分层保留','h2');p('FOXP1优先解释状态重塑，Core-PyroExec关联在等权状态标准化后P=0.494；MDFIC为损伤关联备选，调整后rho=0.950但P=0.0833；DDIT4保留文献驱动假说，项目内PyroExec P=0.751；GSN保留备选，PROS1功能入口支持已撤回。AHNAK/ST6GALNAC3保留氧化代谢或应激转录线索。')
page();p('功能文献为实验提供依据与方向约束','h')
p('IL1R2：以Treg功能为首轮切入点 [1-4]','h2');p('2022年Treg特异Il1r2缺失研究支持其参与肿瘤Treg克隆扩增和效应特征；2024年Tfr研究支持其调节特定免疫环境中的IL-1作用。2025年IL1R2-ENO1-糖酵解/焦亡链来自巨噬细胞及脓毒症。矽肺相关研究也已有其他细胞对象的Il1r2表达报道。研究空间应聚焦矽肺巨噬环境对CD4/Treg及其功能的影响，允许识别保护性反馈。')
p('TTC39C：优先确认干预后的存活与凋亡方向 [5]','h2');p('人A549/NCI-H1299细胞的shRNA实验中，TTC39C降低后Annexin V阳性增加。项目中TTC39C与ApoExec正相关不能直接改写成促进凋亡；首轮要判断在CD4中是驱动、补偿还是其他背景依赖作用。')
p('SKIL：同时观察SnoN蛋白与SKI [6-7]','h2');p('SKIL编码SnoN，SKI为另一基因。小鼠Treg分化研究中，在Arkadia缺失背景下，双靶向Ski/SnoN可救援，单独靶向不能完成该救援，提示特定条件下的冗余。SKIL RNA上调不等同SnoN蛋白或功能增强。PM2.5/ROS/SnoN与肺纤维化已有研究联系，广义关联本身不作为全新机制。')
p('检测与干预可行性核查','h2');p('IL1R2具有小鼠Tfr流式4E2检测先例及厂家小鼠适用性资料。TTC39C有厂家小鼠蛋白WB资源，但原代CD4中的检测和干预效力仍需确认。SnoN原T细胞研究有抗体验证；当前替代试剂需核对物种、SnoN/SnoN2识别和阴性对照，不能把仅人源验证的产品直接视为小鼠适用。具体资源与出处保留在实验可行性表，不将“可采购”写成“已在目标模型有效”。')
page();p('基础实验推荐方案｜总体路线与实施前提','h')
p('推荐先完成小规模、分阶段预实验：以IL1R2为优先对象，TTC39C为损伤方向探索备选，SKIL为Treg状态备选。各候选采用相同的核心分组，但预先指定不同主要终点；出现可重复作用后再投入救援及机制验证。')
tab([['阶段','需要回答的问题','主要产出/推进条件'],['1 模型与表达','巨噬细胞条件是否影响目标T细胞，候选是否可测？','目标亚群身份、候选RNA/蛋白、刺激表型及可接受的基础存活'],['2 候选功能','T细胞内候选干预是否改变巨噬刺激响应？','干预效率、四组交互效应、基础毒性与独立来源重复'],['3 特异性与救援','观察到的作用是否与候选本身相关？','独立干预策略及适当回补/救援，排除干预系统的非特异影响'],['4 上游与机制','哪些巨噬信号连接候选和T细胞表型？','同来源配体-受体证据、近端响应及候选/终点的顺序支持']],[85,195,225])
p('细胞体系的选择','h2');p('优先采用与发现数据一致的小鼠体系。IL1R2与SKIL对应Treg，TTC39C对应活化CD4状态。巨噬来源优先与研究假说匹配；若以骨髓来源巨噬细胞作可行性预筛，应在肺来源巨噬细胞体系中确认关键作用。单细胞ActCD4标签不等同于某一个流式标志，需预先确定活化表型及验证标志。')
p('首轮先限定一种巨噬来源、一种培养方式和一个目标T细胞亚群。可溶性信号假说可从条件培养基开始；接触依赖作用另设隔离与直接接触比较。候选干预首先施加于接收端T细胞，避免同时改变发送端和接收端而无法定位作用。')
p('实施前需确定的参数','h2');tab([['待确定项目','确定依据'],['物种、原代/细胞系、巨噬及T细胞来源','客户资源及发现数据的细胞背景；若改用人源，重新核对标志与试剂适用性'],['刺激剂量、时长、培养比例、采样时点','先导实验中的响应范围、细胞存活及候选蛋白变化；同一比较采用一致条件'],['干预方式、检测抗体及主要终点','目标亚群中的递送效率、特异性和可重复读数；参数确定后再进入正式验证']],[160,345])
p('本方案提供实验设计与决策路径。具体操作参数须在选定细胞体系中优化，不将其他细胞类型文献中的剂量或时长直接当作本项目已验证条件。','small')
page();p('共同设计｜四组比较与必要对照','h')
p('每个候选以“巨噬条件×接收端候选干预”为核心。下表中的候选干预可为经验证的敲低或其他基因干预；非靶向对照须匹配递送方式、处理流程和培养条件。')
tab([['组别','巨噬条件','目标T细胞处理','作用'],['A','对照巨噬','非靶向对照','提供基础参照'],['B','矽尘处理巨噬','非靶向对照','测量巨噬刺激效应'],['C','对照巨噬','候选干预','评估干预本身效应及基础毒性'],['D','矽尘处理巨噬','候选干预','检验候选是否参与刺激响应']],[40,120,145,200])
p('主要比较：交互效应Δ=(D-C)-(B-A)。它回答“候选干预是否改变巨噬刺激引起的响应”，不能仅凭D与B存在差异就忽略干预自身的基础作用。各组维持相同的T细胞活化与支持条件。')
p('模型、来源及检测对照','h2');tab([['控制问题','建议设置与读取方式'],['T细胞基础状态','增加无巨噬刺激的T细胞参照；记录巨噬细胞存活、数量与条件培养基制备过程'],['颗粒直接作用或残留','设置颗粒直接暴露及与残留量相匹配的对照；验证去除效果，两组采用一致处理，不能假定简单处理后已无颗粒'],['目标亚群身份','采用适合活细胞的富集方式，另取平行样本确认FOXP3或活化标志及纯度；固定透化后的细胞不继续用于功能培养'],['接收细胞特异性','流式在目标T细胞门内分析；蛋白等检测前分离目标细胞。混合上清的LDH/细胞因子不能直接归因于T细胞'],['干预效率与基础毒性','同时测目标T细胞RNA及蛋白；必要时设置递送空白/载体对照，结合活细胞数评价'],['检测体系有效性','为死活染色、凋亡执行层及膜电位检测设置适用阳性/阴性对照，区分读数失效与生物学无变化']],[125,380])
p('建议先完成模型及干预可测性，再对三候选开展各自的四组试验。所有扩展对照围绕实际假说配置，首轮不必同时展开全部配体阻断、代谢和双基因干预。','small')
page();p('IL1R2主方案｜Treg功能与损伤耐受','h')
p('实验问题：矽尘处理巨噬细胞条件是否改变Treg中的IL1R2，Treg内IL1R2干预是否影响抑制功能及损伤耐受？已有Treg/Tfr遗传干预研究提供功能依据；IL1R2-ENO1-焦亡机制来自巨噬细胞，迁移到Treg需单独验证。[1-3]')
p('核心必做','h2');tab([['项目','建议读数与解释'],['候选定位','在Treg中测Il1r2 RNA及蛋白，分开记录表面、胞内和可溶性IL1R2。表面阻断或外加可溶蛋白不等同于改变胞内IL1R2'],['主要终点','预先指定时点及Treg:Tconv比例下的Treg抑制能力，以应答T细胞增殖为读数；同批设置无Treg的应答细胞参照'],['功能实验控制','刺激后回收Treg，按等量活细胞进行抑制试验；控制培养基/处理残留对Tconv的直接作用，并报告Treg回收数量'],['辅助终点','FOXP3/CD25表型、活细胞绝对数量、存活及凋亡相关读数；区分Treg数量减少与单位活细胞功能改变'],['干预判定','采用共同四组设计；RNA与蛋白降低或改变应在目标Treg中得到确认，且需记录基础毒性']],[110,395])
p('有功能表型后再追加的机制分支','h2');p('若Treg功能改变同时伴随代谢线索，再测ENO1蛋白及酶活、乳酸或葡萄糖利用，并在条件允许时采用代谢通量检测。候选作用与ENO1的关系需通过遗传或适当药理干预交叉验证；RNA正相关不能替代酶活抑制。若提出蛋白相互作用，应设置输入、阴性及免疫沉淀对照，确认Treg内的实际定位与结合。')
p('若后续出现裂解性死亡表型，再检查目标Treg内相关gasdermin及执行过程，不能仅凭上清炎症因子或LDH变化定义Treg焦亡。[8] IL1A/IL1B与IL1R2的结合属于诱饵受体关系；研究下游信号时还需确认实际传导信号的受体和响应细胞。')
p('结果如何决定下一步','h2');p('若干预改变刺激相关抑制功能且独立来源重复方向一致，推进独立干预及回补。若干预主要改变Treg存活，应先分析细胞数量效应，再判断单位活细胞功能。若干预加重损伤，可按保护性反馈方向研究；不将IL1R2上调预设为促损伤。')
page();p('TTC39C探索方案｜存活、凋亡与线粒体状态','h')
p('实验问题：TTC39C是否参与巨噬刺激下活化CD4细胞的存活及凋亡反应？本项目的优势是ApoExec程序关联；已有肺癌细胞研究观察到敲低TTC39C后凋亡增加，因此首轮不预设敲低具有保护作用。[5]')
p('核心必做','h2');tab([['项目','建议读数与解释'],['细胞与候选确认','明确所用活化CD4的标志和培养状态，确认TTC39C蛋白可检测；避免直接用混合总CD4结果替代目标亚群'],['主要终点','预先指定时点的Annexin V与膜完整性染料联合读数，区分早期与晚期死亡相关状态，并同步报告活细胞数量'],['凋亡执行层','至少结合一项caspase-3/7活化或相应裂解蛋白指标。若只有Annexin V变化，保留为相关表型，不单凭此确证凋亡机制'],['线粒体伴随终点','在明确死活门及适当对照下测膜电位，必要时增加线粒体ROS、ATP或呼吸指标；单一染料变化不等同全部线粒体功能损伤'],['干预与特异性','四组试验先确认递送和靶向效率；获得表型后，用第二种独立靶向策略及适当回补检验作用']],[110,395])
p('测量与时间顺序','h2');p('膜电位测量应保留未处理参照和适用的去极化阳性对照，并报告所分析的活细胞范围。[9] 先导实验可设置有限的早、中、晚观察点，区分候选变化、线粒体读数和死亡发生的先后；正式比较应预先确定主要时点，不能挑选效应最大的时点作为唯一证据。')
p('条件性追加','h2');p('若存活/凋亡效应可重复，再判断线粒体变化是否先于明显死亡，并结合执行层干预或回补评价其联系。焦亡检测只在裂解性死亡及相应执行分子支持时追加，不能因项目主题而将焦亡设为唯一终点。[8]')
p('结果如何决定下一步','h2');p('敲低后死亡降低或升高均可形成研究方向，前提是排除基础毒性及干预不足。若仅有线粒体读数改变而主要凋亡终点无效，应如实记录，将新的线粒体功能问题作为后续预设假说，而非追认原主要终点验证成功。')
page();p('SKIL备选方案｜Treg状态与SKI补偿','h')
p('实验问题：巨噬条件是否改变Treg的SnoN蛋白及稳定性，SKIL干预是否影响该状态响应？SKIL编码SnoN，SKI为另一基因。既往小鼠Treg分化研究中的双靶向救援提示特定背景下存在冗余，不能推广为任何条件下单基因干预均无效。[6]')
p('核心必做','h2');tab([['项目','建议读数与解释'],['候选与补偿分子','并行测SKIL RNA、SnoN蛋白及SKI RNA/蛋白；重点检查RNA和蛋白变化是否一致'],['主要终点','以成熟Treg体系为起点，预先指定刺激后FOXP3维持率作为状态终点，同时记录活细胞绝对数量及存活'],['功能确认','有稳定状态变化时，用等量活Treg的抑制试验确认功能影响；FOXP3比例不直接等同抑制能力'],['四组单基因试验','先在对照/矽尘巨噬条件下比较非靶向与SKIL干预，确认SnoN实际改变及递送毒性'],['蛋白检测','核对抗体的物种与应用验证；分清SnoN/SnoN2识别及条带特异性，不能只依据RNA判断蛋白功能']],[110,395])
p('分化与补偿分支按结果追加','h2');p('若主要问题转向Treg分化，可另设初始CD4诱导Treg的独立实验，预先指定FOXP3阳性产率及活细胞数量；不要与成熟Treg稳定性实验合并解释。可依据结果增加TGF-β/SMAD相关近端读数，但不因SKIL升高便推定该通路已经增强。')
p('若SKIL单独干预效率充分而效应有限，且有SKI补偿或相关机制依据，再在两种巨噬条件下比较非靶向、SKIL单干预、SKI单干预及双干预。只有双干预出现表型时，应定位为联合调控或补偿问题，不能把效应全部归于SKIL单基因。')
p('结果如何决定下一步','h2');p('SnoN蛋白可测、干预有效且状态/功能响应可重复时推进救援。若RNA改变而蛋白未改变，先解决检测和调控层面的解释；若只有细胞数量下降，应先排除死亡选择造成的表观FOXP3比例变化。线粒体损伤和细胞死亡作为伴随终点，不替代Treg状态主问题。')
page();p('验证推进｜救援、统计与投入决策','h')
p('特异性与上游定位','h2');p('先用独立靶向策略和适当回补确认候选相关作用，回补设对应载体对照并避免异常高表达造成新表型。随后按实际作用方向选择上游干预：核对巨噬来源、配体、完整受体及近端响应，通过发送端或接收端阻断/补充，检测候选与主要终点，并检验候选回补能否支持作用顺序。')
p('上游选择保持同一来源和接收状态。IL1R2优先区分IL-1捕获与胞内作用；SKIL可结合表达相容的信号及TGF-β相关机制提出候选；TTC39C的上游连接不足，应先取得基因功能表型再定位，不为形成完整链而指定缺乏支持的配体。')
p('独立重复、主要比较与结果报告','h2');p('优先以独立动物/供体来源作为生物学重复单位。同一来源细胞的分装培养、技术复孔及重复测量不增加独立样本量；记录巨噬和T细胞各自来源。各独立来源尽量覆盖同一四组设计，平衡批次、随机分配处理并对检测分析盲法编码。')
p('每候选预先指定一个主要终点及一个主要交互对比Δ。模型包含巨噬条件、候选干预及交互项；配对来源或批次按设计纳入，重复时间点采用相应重复测量分析。单个预设主要对比无需FDR校正；同时对多个候选或终点提出联合验证结论时，预先定义同类检验家族并报告FDR。名义P<0.05与FDR显著分别表述。')
p('先导实验用于估计变异、可检测效应和可接受基础存活；正式样本量根据主要终点与设计计算，不预设固定重复数即可保证效力。报告效应量、区间、原始P、独立来源数及全部预设主要结果；不以技术孔数替代n，不靠更换阈值反复寻找阳性。')
tab([['观察结果','下一步决策'],['模型尚无稳定T细胞表型，或目标亚群/蛋白不清','先完善模型、检测及身份确认；暂不据此判断基因无价值'],['干预效率不足或出现明显基础毒性','优化递送及干预条件，检查非靶向/载体对照；不直接作机制解释'],['刺激表型明确、干预充分，主要终点仍无一致变化','记录阴性，调整候选或预设研究问题；不以事后替换终点追认成功'],['主要效应跨独立来源可重复，基础毒性可解释','推进独立干预、回补，再开展上游及机制验证'],['只有探索性终点改变，主要终点未达到支持','保留探索线索，在新的预设实验中确认，区分探索与验证结论']],[235,270])
p('正式实施前形成一页记录表：细胞来源/批次、四组分配、干预方式、剂量与时长、主要时点和终点、独立n、分析模型及排除规则。表型达到上述推进条件后，再决定是否进入动物水平验证。','small')
page();p('讨论：证据边界与研究空间','h')
p('项目支持候选选择，不预先决定湿实验方向','h2');p('疾病差异、程序关联、网络先验和干预功能是不同证据层级。IL1R2具备较清楚的疾病差异和Treg功能研究基础；TTC39C提供损伤相关程序切入点；SKIL具备Treg状态和调控背景。三者均有预实验价值，但不对应已量化的成功概率。')
p('样本与模型边界','h2');p('主队列为3+3个样本；第二队列各时间点1+1，联合分析含原队列；时间队列为TLike范围。精确枚举保证既定排列空间计算完整，不证明可交换性或模型假设成立。MDFIC等强效应而P有限的结果不能全部归因于分辨率。')
p('生物学与技术解释','h2');p('转录评分不能替代线粒体膜电位、呼吸、蛋白裂解或死亡实测。富集结果中的Coronavirus disease、synapse等名称可能由共享核糖体/翻译成员驱动，不表示感染或神经细胞功能；免疫球蛋白/组织背景成员需结合背景敏感性解释。原始细胞图谱的未确定标签保持原状，不因候选需要而重定义。')
p('PROS1-AXL的胶质瘤文献背景不能直接用于小鼠CD4。TYRO3/MERTK在现有覆盖规则下未补回其功能入口；RNA低检出不是蛋白不存在的证明。其他网络边也仅保留计算支持，不视为已经验证。状态DPT是横断面转录排序，不能代替谱系追踪或真实损伤时间。')
p('创新性落在具体细胞关系和作用环节','h2');p('IL1R2-ENO1-焦亡、SKIL/SnoN-TGF-β-Treg及TTC39C-肿瘤细胞存活均已有研究。新的问题应落在矽肺巨噬环境如何改变CD4/Treg候选及功能，以及候选是否参与损伤耐受或细胞死亡。定向检索未取得完整链直接功能证明，不等同全球首次声明。')
page();p('参考研究与可复核材料','h')
refs=[('[1]','IL1R2 increases regulatory T cell population in the tumor microenvironment by enhancing MHC-II expression on cancer-associated fibroblasts. JITC, 2022.','https://pmc.ncbi.nlm.nih.gov/articles/PMC9438093/'),('[2]','Loss of T follicular regulatory cell-derived IL-1R2 augments germinal center reactions via increased IL-1. JCI Insight, 2024.','https://insight.jci.org/articles/view/174005'),('[3]','Critical Role of IL1R2-ENO1 Interaction in Inhibiting Glycolysis-Mediated Pyroptosis for Protection Against Lethal Sepsis. Advanced Science, 2025.','https://pmc.ncbi.nlm.nih.gov/articles/PMC12533375/'),('[4]','Epithelial NLRP3 drives silica-induced lung injury and fibrosis through IL-18 and pro-fibrotic neutrophil recruitment. Particle and Fibre Toxicology, 2026.','https://pmc.ncbi.nlm.nih.gov/articles/PMC13214385/'),('[5]','Ttc39c is a potential target for the treatment of lung cancer. BMC Pulmonary Medicine, 2022;22:391.','https://pmc.ncbi.nlm.nih.gov/articles/PMC9615393/'),('[6]','Arkadia-SKI/SnoN signaling differentially regulates TGF-β-induced iTreg and Th17 cell differentiation. JEM, 2021.','https://pmc.ncbi.nlm.nih.gov/articles/PMC8421263/'),('[7]','PM2.5 regulates TGF-β1/Smads-mediated pulmonary fibrosis via ROS/SnoN in vitro and in vivo. Int Immunopharmacol, 2025.','https://pubmed.ncbi.nlm.nih.gov/40154177/')]
refs.extend([('[8]','Galluzzi et al. Molecular mechanisms of cell death: recommendations of the Nomenclature Committee on Cell Death 2018. Cell Death & Differentiation, 25:486-541.','https://www.nature.com/articles/s41418-017-0012-4'),('[9]','Thermo Fisher Scientific. TMRM technical guidance and membrane-potential controls.','https://www.thermofisher.com/order/catalog/product/T668/faqs')])
for n,title,url in refs:p(n+' '+title+'\n'+url,'small')
p('数据与图件复核范围','h2');p('v19独立复核报告确认其实际检查范围：2,852个评分值、2,544个相关系数、2,268个穷举置换P、34个设计、221项程序效应及193个原BH家族。276个Monte Carlo P、原始细胞层edgeR重拟合、全套注释/SCENIC等不在该独立复算范围。','small')
p('本次新增通路、比例及代表图由分析端另行核查，不能扩大为全部分析已独立重跑。完整统计、图值、方法、文献/试剂依据及校验清单保存在v20配套材料。原v19与更早版本保留，候选文献优先于旧实验排序表述。报告生成不代表已向外部发送。','small')
# Layout uses landscape pages for wide scientific panels and portrait for long text/tall panels.
pdf=D/'矽肺免疫全景_CD4候选与验证_报告_v26.pdf'
def footer(c,d):
 w,h=c._pagesize;c.setStrokeColor(colors.HexColor('#c5d8de'));c.line(42,34,w-42,34);c.setFont('CN',8);c.setFillColor(colors.HexColor('#59707a'));c.drawString(42,22,'矽肺免疫微环境与CD4候选研究');c.drawRightString(w-42,22,str(d.page))
class Doc(BaseDocTemplate):
 def afterFlowable(self,f):
  if isinstance(f,Paragraph) and f.style.name in ['h','h2']:figpages.append({'page':self.page,'heading':f.getPlainText()})
doc=Doc(str(pdf),pagesize=(595.28,841.89),leftMargin=45,rightMargin=45,topMargin=40,bottomMargin=45)
doc.addPageTemplates([PageTemplate(id='portrait',pagesize=(595.28,841.89),frames=[Frame(45,45,505.28,756.89,leftPadding=0,rightPadding=0,topPadding=0,bottomPadding=0)],onPage=footer),PageTemplate(id='landscape',pagesize=(841.89,595.28),frames=[Frame(45,45,751.89,510.28,leftPadding=0,rightPadding=0,topPadding=0,bottomPadding=0)],onPage=footer)])
doc.build(story)
(O/'report_text.txt').write_text('\n\n'.join(textrecords),encoding='utf-8');(O/'report_page_index.json').write_text(json.dumps(figpages,ensure_ascii=False,indent=2),encoding='utf-8');(O/'report_delivery.json').write_text(json.dumps({'pdf':str(pdf),'sha256':hashlib.sha256(pdf.read_bytes()).hexdigest(),'figures':len(manifest),'explanatory_displays':3,'total_visuals':len(manifest)+3},ensure_ascii=False,indent=2),encoding='utf-8');print('REPORT',pdf)