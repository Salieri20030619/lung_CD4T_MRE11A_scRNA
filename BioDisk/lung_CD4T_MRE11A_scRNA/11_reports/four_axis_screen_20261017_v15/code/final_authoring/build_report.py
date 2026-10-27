from pathlib import Path
import pandas as pd,numpy as np,json,html,hashlib
import fitz
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,PageBreak,Table,TableStyle,Image,KeepTogether
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors
W=Path('C:/Users/31937/lung_four_axes_20261017_v15');T=W/'tables';O=W/'deliverables';F=W/'figures';R=W/'report_render';R.mkdir(exist_ok=True)
rd=lambda p:pd.read_csv(T/p,sep='\t');pdfmetrics.registerFont(TTFont('CJK','C:/Windows/Fonts/ARIALUNI.ttf'))
styles={k:ParagraphStyle(k,fontName='CJK',fontSize=fs,leading=ld,spaceAfter=sp,wordWrap='CJK',textColor=colors.HexColor(cl),keepWithNext=k in ['h1','h2']) for k,fs,ld,sp,cl in [('body',10.5,16.5,9,'#172e3b'),('h1',18,25,14,'#153e53'),('h2',12,19,8,'#26778b'),('table',9,13,0,'#172e3b'),('th',9,13,0,'#ffffff'),('caption',9,13,9,'#49616d')]}
story=[];pagemap=[]
class H(Paragraph):
 def __init__(self,t):super().__init__(html.escape(t),styles['h1']);self.heading=t
class Doc(SimpleDocTemplate):
 def afterFlowable(self,f):
  if hasattr(f,'heading'):
   key='p'+str(len(pagemap));self.canv.bookmarkPage(key);self.canv.addOutlineEntry(f.heading,key,0);pagemap.append({'title':f.heading,'page':self.page})
def section(t):
 if story:story.append(PageBreak())
 story.append(H(t))
def p(t,k='body'):story.append(Paragraph(html.escape(t),styles[k]))
def tab(head,rows,widths=None):
 data=[[Paragraph(html.escape(str(v)),styles['th' if i==0 else 'table']) for v in row] for i,row in enumerate([head]+rows)];tb=Table(data,colWidths=widths or [762/len(head)]*len(head),repeatRows=1,hAlign='LEFT');tb.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#183e51')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#eff4f6')]),('VALIGN',(0,0),(-1,-1),'TOP'),('GRID',(0,0),(-1,-1),.25,colors.HexColor('#ccd9df')),('TOPPADDING',(0,0),(-1,-1),4),('BOTTOMPADDING',(0,0),(-1,-1),4)]));story.extend([tb,Spacer(1,9)])
def fig(name,caption,h=260):
 im=fitz.Pixmap(str(F/name));s=min(762/im.width,h/im.height);story.append(KeepTogether([Image(str(F/name),im.width*s,im.height*s),Spacer(1,5),Paragraph(html.escape(caption),styles['caption'])]))
section('四证据线联合筛选与候选核查')
p('v15｜2026年10月17日｜中文内部分析与独立评审材料','h2')
p('冻结规则的四线交集得到1个对象：DUSP5–Treg；追加诊断提示损伤关联不稳健，尚不能据此确定客户机制主线。','h2')
p('本轮从完整差异检验基因出发，分别评估疾病差异、巨噬信号、损伤相关转录程序及创新性。Naive、ActCD4、Treg合计15,913条基因–状态记录中，1,297条达到原始差异P<0.05，1,161条同时通过预设信号质量检查。所有主筛选使用同一状态、同一发送端范围及六只独立小鼠。')
p('在数据层，DUSP5–Treg是唯一同时通过疾病差异、损伤关联、巨噬关联及组别调整要求的对象。文献支持DUSP5的T细胞功能，但其免疫调控、硅尘响应和其他细胞中的焦亡关系均已有报道；创新空间只能限定在特定矽肺Treg和来源关系。')
p('最终来源核查中，CDH1–KLRG1或CD103接触关系具有生物学依据；另三条先验边不足以作为本项目的直接细胞间信号。六样本去污染诊断保留了CDH1–DUSP5联系，却使DUSP5–PyroExec的合并P由0.01667变为0.175，组别调整P由0.02778变为0.36111。该新增诊断不回写冻结交集，但阻止把筛选阳性升级为已定主线。')
tab(['层级','完成结果'],[['冻结规则四线交集','1个：DUSP5–Treg；原统计与判断保持可追溯'],['宽松计算交集','34个基因–状态组合，包含上述1个；另外33个未通过严格数据层要求'],['当前可建议直接投入的机制主线','0个；主要未决点为损伤程序关联对背景处理敏感'],['任务完成状态','固定计划计算、创新性、候选来源核查和数值验收均已完成；不等于机制选题成功']],[190,572])
section('1 数据与方法：固定对象、质量与疾病差异')
p('GSE314943为发现队列，保留3个Saline与3个Silica独立小鼠。沿用冻结QC、去双细胞和细胞状态定义；不重聚类或按候选表达挑细胞。主筛选为Naive、ActCD4、Treg；Core仅作汇总背景。本轮不使用细胞数替代样本量，不把三个状态视作独立验证队列。')
p('完整差异表采用既有原始计数pseudobulk edgeR：~group设计、filterByExpr、过滤后重新计算文库量、TMM、稳健离散度与准似然F检验。原始目标效应P<0.05为差异准入标准，保留全部logFC、P和FDR；未设新的统一倍数或FDR门槛。原始差异家族不因本轮交集缩小而重新计算。')
p('质量检查使用原始CP10k、BG1（预定背景面板UMI/原始总UMI≤1%）及T身份/BG1。T身份要求CD3D、CD3E、CD3G、TRAC中至少2项检出。在表达较高的一组，至少2/3小鼠具有候选≥2 UMI的T身份/BG1细胞；BG1与T身份/BG1的疾病均值方向保持，全部六样本可评价，留一小鼠描述性均值差方向至少5/6一致。预设上皮、髓系及免疫球蛋白背景面板对象不进入主候选。')
p('使用冻结canon对照处理MRE11A、ATP5与PARK2等历史别名，保留source_canonical及Ensembl标识。每个状态有1个ENSA对象不在冻结的无歧义计数映射中，保留差异记录并标记不可评估；未改写其数据值。')
ro=rd('01_quality/sample_scope_roster.tsv');rows=[]
for sc in ['Naive','ActCD4','Treg']:
 for b in ['All','BG1','T_BG1']:
  z=ro[(ro.scope==sc)&(ro.branch==b)].set_index('sample');rows.append([sc,b,','.join(str(z.loc[f'Saline_{i}','cells']) for i in [1,2,3]),','.join(str(z.loc[f'Silica_{i}','cells']) for i in [1,2,3])])
tab(['状态','分支','对照各样本细胞数','病例各样本细胞数'],rows,[95,110,278,279])
section('1.1 数据与方法：损伤、巨噬及创新性')
p('损伤主程序固定为PyroExec、MitoROS、MitoISR，背景程序为Inflam、ApoExec、OXPHOS、Mitophagy。评分为模块平均log1p(CP10k)减匹配对照均值。每次检验均从模块成员及对照中剔除当前候选，至少保留3个成员、20个对照。全部受检基因均分析，不先取显著基因再缩小损伤检验家族。')
p('样本秩关联采用合并720种配对及疾病组别调整36种组内配对。主筛选要求两模型P<0.05、rho同号、差异方向×关联方向>0、BG1方向一致及合并模型留一方向至少5/6一致。BH分别在状态×分支×模型的主3程序和背景4程序完整家族内计算。')
p('巨噬来源分别为AMlike与AllMac。采用冻结NicheNet完整先验中每配体前200个正权重靶基因；病例至少2/3发送端样本检出配体、2/3接收端样本检出兼容受体。要求同来源配体或同状态受体具有疾病差异P<0.05。每条联系按同六只小鼠检验，主筛选要求合并及组别调整P<0.05、rho同号、BG1方向一致、留一方向至少5/6一致。配体有疾病差异时检查其效应方向与候选方向的相容性。BH使用来源×状态×模型内全部受检先验联系。[M01]')
p('先验权重无激活/抑制方向，也不是P值。形成数据交集后，逐条核查其受体复合物、发送端原始计数、其他细胞来源及是否真正具备跨细胞作用资格。创新性仅对主数据交集进行原始文献核查；A为完整同情境机制已知，B为有基础但特定联系未建立，C为基础不足，D为直接矛盾或来源未解决。未找到文献不等于证明不存在。')
p('本轮唯一数据候选触发额外来源诊断：对37,155个全肺细胞、29,809个输入特征按六样本分别运行decontX，使用冻结Broad及既有CD4/Mac精细标签，seed42、maxIter500。保留原始、校正计数/原始分母、校正计数/校正分母三分支；不替代原筛选结果。该诊断在候选出现后追加，单独记录，未伪称为最初冻结的门槛。[M02]')
section('2 分别分析后，仅DUSP5–Treg进入严格数据交集')
ct=rd('05_intersection/stage_counts.tsv').set_index('scope').loc[['Naive','ActCD4','Treg']]
tab(['状态','受检基因','差异P<0.05','差异+质量','再通过损伤线','再通过巨噬线','同时通过'],[[s,int(r.tested),int(r.DE),int(r.Q_DE),int(r.damage_primary),int(r.mac_primary),int(r.data_primary)] for s,r in ct.iterrows()],[88,90,116,116,116,116,120])
fig('F01_screen_axes.png','图1　两条功能线分别在差异与质量合格对象中计数，再取同状态交集。宽松交集34个，包含严格数据交集DUSP5–Treg；Core不参与主筛选。所有范围仍来自同一批六只小鼠。',270)
p('损伤关联完整输出686,924行，巨噬联系完整输出667,420行，包含背景范围、两模型及不可评估项。Naive和ActCD4均无同时通过两线的对象；Treg有1个。此计数反映操作性筛选结果，不是独立验证次数。')
section('3 DUSP5在Treg下调，信号质量分支保持方向')
q=rd('05_intersection/primary_data_candidates.tsv').iloc[0]
p(f'DUSP5在Treg中log2FC={q.logFC:.4f}，原始P={q.PValue:.6f}，FDR={q.FDR:.6f}，属于名义差异。六个原始样本中均有多UMI阳性细胞；T身份/BG1保留24、9、10个对照及145、142、139个病例细胞，未丢失独立样本。')
fig('F02_DUSP5_raw_evidence.png','图2　C1–C3为三只对照，S1–S3为三只病例。左：相同样本在质量分支中的表达；中：候选剔除后的焦亡相关评分；右：同一AllMac来源CDH1与Treg DUSP5。每点为小鼠，原始完整家族BH与组别调整BH分别标注。',260)
p('DUSP5–PyroExec合并rho=-0.9429、P=0.016667、BH=0.377044；调整组别rho=-0.9585、P=0.027778、BH=0.989070。方向符合“DUSP5降低伴随较高程序评分”的筛选规则，但并不确立保护作用。该对象没有同时获得MitoROS或MitoISR主条件支持，因此不能称为已筛定的线粒体损伤基因。')
section('4 四条数值上游联系中，仅CDH1接触关系具备适用基础')
p('以下均为Treg中的DUSP5，不能用Core的其他上游联系补足本节。原始统计完整保留；文献资格核查用于区分预测联系与可解释的细胞间作用。')
li=rd('05_intersection/DUSP5_Treg_primary_links.tsv');rows=[]
for r in li.itertuples():rows.append([r.sender,r.ligand,f'{r.rho:.4f}',f'{r.P:.6f} / {r.BH:.6f}',f'{r.adj_P:.6f} / {r.adj_BH:.6f}'])
tab(['来源','先验发送分子','合并rho','合并P / 完整BH','调整P / 完整BH'],rows,[110,112,100,220,220])
tab(['分子','资格核查与结果'],[['BMP1','前胶原C蛋白酶，不能凭名称按经典BMP受体配体解释；AMlike中12个原始UMI均为单UMI阳性细胞。保留数值，但不用于直接配体受体机制。[L09]'],['OCLN','与TGFBR的相互作用证据不足以证明外源跨细胞配体作用；保留为先验相关，不建立巨噬→Treg连接。[L10]'],['H2-Q1','APLP2相关证据涉及MHC转运；Treg未检出CD8A/CD8B1同细胞组合。AllMac仅6个UMI且均单UMI，不作为本链的可靠发送信号。[L11]'],['CDH1','CDH1与KLRG1及CD103接触有基础。病例Treg KLRG1检出50/178、64/156、68/161；ITGAE/ITGB7共同检出53/178、41/156、52/161。仅保留接触假说，未证明其驱动DUSP5。[L07–L08]']],[105,657])
p('AllMac CDH1自身疾病差异P=0.106507；本联系进入筛选依赖接收端差异：KLRG1 log2FC=2.1109、P=0.000201、FDR=0.020563；ITGAE P=0.001956。不能把该分支写为“巨噬CDH1显著上调”。')
section('5 DUSP5具有功能基础，创新性仅限特定细胞情境')
tab(['已有证据','对本项目的限定'],[['T细胞与ERK调节；Th17/Treg平衡','T细胞模型及关节炎模型已建立DUSP5相关免疫功能。不能宣称首次发现DUSP5调节T细胞或Treg。[L01–L03]'],['小细胞肺癌中的DUSP5-A20-NLRP3焦亡','已有敲低、蛋白和相互作用证据。该研究中DUSP5促进焦亡相关过程；本项目Treg是反向转录关联，不能把文献当作同方向验证。[L04]'],['硅尘暴露的表达响应','A549表达谱中已涉及DUSP5。疾病关键词缺少完全相同机制不等于基因在硅尘背景全新。[L05]'],['KLRG1+Treg与线粒体改变','衰老Treg中已有线粒体改变与KLRG1状态研究；需区分状态标记、DUSP5效应与巨噬来源作用。[L06]'],['巨噬CDH1及CDH1-KLRG1接触','巨噬细胞可以表达CDH1，受体结合也有基础；已发表巨噬CDH1体内研究亦有阴性结果，不能外推为必需致病源。[L07–L08]']],[250,512])
p('创新性评定为B（具体情境待验证），依据是已有功能基础与尚未建立的特定矽肺来源/状态联系。本轮检索未定位到完整“巨噬CDH1—Treg DUSP5—焦亡”机制的直接研究；这是有限检索下的判断，不是文献不存在的证明。')
p('可提出的研究问题是这一具体联系是否存在，而不是直接承诺其成立。该创新性等级保留，但不能抵消损伤评分的诊断不稳定，也不能凭新颖性把先验不合格的BMP1/OCLN/H2-Q1重新纳入主线。完整检索式、日期和原始论文链接见06_novelty。')
section('6 追加模型诊断保留上游联系，但削弱损伤关联')
fig('F03_model_diagnostic.png','图3　六样本decontX诊断，未改变冻结交集。原始与校正CP10k及评分均由相同细胞计算；右图为六样本合计UMI保留比例，仅描述模型对各转录成分的调整，不作细胞级推断。',250)
tab(['关系','原始合并P','模型校正合并P','原始调整组别P','模型校正调整P'],[['Treg DUSP5–PyroExec','0.016667','0.175000','0.027778','0.361111'],['AllMac CDH1–Treg DUSP5','0.002778','0.002778','0.027778','0.027778']],[238,131,131,131,131])
p('校正计数采用原始分母或校正分母时，上述秩检验结论一致。Treg DUSP5总UMI由677调整为671.48，GSDMD由93调整为91.37；IL1B由67调整为8.65。模型主要提示焦亡评分的组成成分需要警惕，不能把它概括为DUSP5本身是假信号，也不能据模型保留GSDMD即确认其内源性。')
p('DUSP5–PyroExec校正后rho=-0.6571，调整组别rho=-0.6162，方向仍为负，但统计证据减弱。原BH不在候选子集中重算；校正后单一损伤诊断终点不需要FDR，四条上游诊断另成家族，并明确是选择后的敏感性分析。')
# Complete exploratory inventory, including the one strict data hit.
e=rd('05_intersection/exploratory_data_candidates.tsv').sort_values(['scope','gene'])
for j in range(0,len(e),17):
 section(f'7 宽松交集完整名单（{j+1}–{min(j+17,len(e))}/34）')
 p('宽松层包含DUSP5，共34个对象；其余33个未同时通过损伤与巨噬的组别调整等严格条件。该表不是主候选推荐榜，也未对其余33个逐一完成创新性审查；原计划只要求主数据交集进入该步骤。')
 rows=[]
 for r in e.iloc[j:j+17].itertuples():rows.append([r.gene,r.scope,f'{r.logFC:.3f}',f'{r.PValue:.5g}',f'{r.FDR:.4f}',str(r.damage_primary),str(r.mac_primary)])
 tab(['基因','状态','log2FC','差异P','差异FDR','严格损伤数','严格上游数'],rows,[125,90,90,117,110,115,115])
section('8 讨论与最终决策：保留筛选阳性，不升级为已定主线')
p('冻结规则回答的是“哪些对象在本队列中同时满足四条证据线的操作性条件”。答案为DUSP5–Treg。追加诊断回答的是“这个筛选阳性是否足以作为稳定损伤主线”。当前答案仍是否定的：损伤程序关联对模型背景处理敏感，而另三条数值上游边缺少适用的跨细胞语义。两种判断应同时保留，不能把后验诊断伪装成预设门槛，也不能因原交集阳性而忽略新增风险。')
p('本轮主要结果来自六只小鼠，均为探索性推断。全基因与全先验筛选的原始P用于发现，FDR用于呈现完整家族证据；DUSP5关键结果未通过原完整家族FDR，但本次不推荐直接投入并非单由FDR决定。来源资格、表型代理和背景敏感性共同影响研究价值。')
p('PyroExec、MitoROS、MitoISR是转录程序，不是实测焦亡、线粒体膜电位或ROS。质量过滤、去自关联和decontX均不能替代细胞纯化、蛋白活化或功能终点；decontX依赖注释、未使用空滴，诊断变化不能独立证明所有被调整转录均为污染。相同数据上的重复建模也不是独立验证。')
tab(['最终记录','决定'],[['DUSP5–Treg','冻结四线交集阳性；创新性B；损伤诊断不稳健，保留待核查'],['成熟机制主线','本轮未形成；不建议作为已筛定主线直接提交客户投入'],['既有IKZF2、STAT1、CXCL16','未获得保留名额；按相同规则记录在排除账本，不用Core或跨来源结果补交集'],['后续是否继续','本轮停止于预定筛选和有依据的候选核查完成；不继续改门槛或扩展筛选追求阳性']],[190,572])
p('如果后续仍研究DUSP5，应先独立确认Treg内死亡/损伤表型以及DUSP5与该表型的方向，再判断CDH1接触是否影响它。本报告不将这些尚未完成的实验写为当前已有支持。')
section('9 完成范围、验收与统计附件')
tab(['检查项目','实际完成'],[['固定计划','PLAN_FROZEN.json哈希未变；中文计划另附。所有筛选层及NA规则保留'],['完整计算','15,913条主状态基因记录；686,924行损伤检验；667,420行巨噬检验；原P/FDR及全家族保存'],['独立数值复核','48个BH家族通过；24项原始计数损伤精确检验、8项来源精确检验通过；48条All/BG1样本细胞数一致'],['来源诊断','六样本正常退出；原始DUSP5表达及Pyro评分与旧输入回读一致；校正计数非负'],['版本保护','v13/v14保持原目录；v15独立保存，未覆盖原结果'],['未重跑/未验证','上游QC、完整edgeR拟合、原SCENIC建网未重跑；未完成独立队列或实验机制验证']],[200,562])
tab(['目录','内容'],[['01_quality','完整质量门槛、样本表达与细胞数；ENSA不可评估记录'],['02_DE','完整原edgeR结果与原始标识'],['03_damage','全基因主/背景程序检验与同状态汇合证据'],['04_macrophage','完整来源内先验联系、统计、受体及疾病支持'],['05_intersection','完整交集/排除账本、DUSP5样本值、最终分层决策'],['06_novelty','创新性原始论文、检索记录、受体资格、来源身份和诊断统计'],['provenance','计划、环境、输入哈希、日志、退出状态、数值验收']],[180,582])
# Bibliography: source-specific short summaries above, no long quotations.
refs=pd.read_csv(W/'literature_evidence.tsv',sep='\t')
for start in [0,5,10]:
 section('参考文献与方法来源'+(f'（{start//5+1}）' if start else ''))
 for r in refs.iloc[start:start+5].itertuples():
  p(f'[{r.reference_id}] {r.title}','caption');story.append(Paragraph(f'<link href="{html.escape(r.url,quote=True)}" color="#167b83">{html.escape(r.url)}</link>',styles['caption']));p(r.assessment_zh,'caption')
 p('检索日期：2026年10月17日。仅以原始研究或官方方法资料作为所列判断依据。')
def footer(c,d):
 c.saveState();c.setFillColor(colors.HexColor('#153e53'));c.rect(0,583,842,12,fill=1,stroke=0);c.setFont('CJK',8);c.drawString(40,22,'四证据线联合筛选｜v15｜内部分析与再次评审');c.drawRightString(802,22,str(d.page));c.restoreState()
f=O/'four_axis_screen_report_v15.pdf';Doc(str(f),pagesize=(842,595),leftMargin=40,rightMargin=40,topMargin=35,bottomMargin=40,title='四证据线联合筛选与候选核查 v15').build(story,onFirstPage=footer,onLaterPages=footer)
d=fitz.open(f)
for i,page in enumerate(d):
 assert '\ufffd' not in page.get_text();assert all(b[0]>=0 and b[1]>=0 and b[2]<=842 and b[3]<=595 for b in page.get_text('blocks')),(i,page.get_text());page.get_pixmap().save(R/f'page_{i+1:02d}.png')
(W/'report_page_map.json').write_text(json.dumps(pagemap,ensure_ascii=False,indent=2),encoding='utf-8');print('REPORT_PAGES',len(d));print(json.dumps(pagemap,ensure_ascii=False))