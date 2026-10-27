from pathlib import Path
import json,html
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,PageBreak,Table,TableStyle
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors
W=Path('C:/Users/31937/lung_four_axes_20261017_v15');O=W/'deliverables';O.mkdir(exist_ok=True)
pdfmetrics.registerFont(TTFont('CJK','C:/Windows/Fonts/ARIALUNI.ttf'));body=ParagraphStyle('body',fontName='CJK',fontSize=10.5,leading=17,spaceAfter=10,wordWrap='CJK');title=ParagraphStyle('title',parent=body,fontSize=18,leading=25,textColor=colors.HexColor('#173e53'),spaceAfter=16);head=ParagraphStyle('head',parent=body,fontSize=12,leading=19,textColor=colors.HexColor('#28758a'),spaceBefore=8)
story=[]
def p(t,s=body):story.append(Paragraph(html.escape(t),s))
p('四条证据线联合筛选：冻结实施计划',title);p('v15｜2026年10月17日｜先冻结规则，再分析交集')
p('目标与统计单位',head);p('寻找同时具有CD4疾病差异、巨噬上游联系、损伤相关转录程序关联及具体机制创新空间的基因。以GSE314943的3只对照与3只硅尘小鼠为独立样本，分别分析Naive、ActCD4、Treg；Core只作背景汇总，不替代状态内证据。保留全部edgeR受检基因，不给既有候选保留名额。')
p('共同质量门槛',head);p('使用冻结原始计数、BG1背景负担≤1%及T身份支持（CD3D/E/G、TRAC中至少2个检出）。在表达较高的一组，至少2/3小鼠具有候选≥2 UMI的T身份/BG1细胞；六个样本均可评价，BG1及T身份/BG1中的疾病均值方向与原差异方向一致；逐一剔除小鼠的均值方向至少5/6一致。预设上皮/髓系/免疫球蛋白背景面板对象不进入主候选。此为转录信号可信度操作标准，不等于证明蛋白内源性。')
p('证据线一：疾病差异',head);p('按相同CD4状态复用完整原始计数pseudobulk edgeR目标效应，原始P<0.05作为硬筛选标准，报告logFC、P、FDR及完整检验家族；不要求FDR通过，不新增统一倍数门槛。保留原标识及显式基因别名映射，缺失映射记为不可评估。')
p('证据线二：损伤相关程序',head);p('主程序为PyroExec、MitoROS、MitoISR；Inflam、ApoExec、OXPHOS、Mitophagy为背景解释。对每个候选，将自身从模块成员和匹配对照中同时剔除，再计算样本评分；至少保留3个成员、20个对照。全部基因分别计算合并720种配对和组别调整36种组内配对的双侧精确秩关联。')
p('主候选要求合并P<0.05、组别调整P<0.05且rho同号，差异方向×关联方向>0，BG1关联方向一致，合并模型留一方向至少5/6一致。仅通过合并模型者保留为探索性层。BH在每状态、分支、模型及主3程序/背景4程序的完整家族内计算。程序反映转录变化，不能代替实际损伤或死亡测量。')
story.append(PageBreak())
p('证据线三：巨噬信号',title);p('分别分析AMlike、AllMac；不得跨来源拼接。使用完整NicheNet先验每配体前200个正权重靶基因，保留权重和排名。病例侧至少2/3发送端样本检出配体，至少2/3接收端样本检出兼容受体。要求同来源配体或同状态受体具有原始P<0.05疾病差异。')
p('对所有合格先验配体-靶基因联系计算同六只小鼠的关联。主候选要求合并及组别调整P<0.05、rho同号、BG1方向一致、留一方向至少5/6一致；配体有疾病差异时，配体差异方向×关联方向须与候选差异一致。仅受体差异支持时单独标注，不补造配体差异。BH采用来源×状态×模型内全部受检配体-靶基因联系。先验无方向性，不将权重解释为因果或激活证据。')
p('对数据交集中的对象再核查完整受体复合物、发送端条件，并比较全肺其他细胞来源。转录共检出不能替代受体蛋白、分泌或跨呈递。')
p('证据线四：创新性',head);p('数据交集形成后，对全部主数据候选查原始文献，按基因+CD4/T细胞+线粒体/焦亡以及基因+矽肺/巨噬细胞查询。A：相同完整机制已建立；B：有功能基础，但具体疾病/来源/状态联系未验证；C：未找到直接基础，创新性和机制价值不确定；D：直接矛盾或来源身份未解决。只有B级且通过最终受体/来源核查者进入计算层面的主候选。未检索到不等于文献不存在。')
p('交集为空时的处理',head);p('不改门槛、不反向为原候选寻找支持。保留主交集、探索性交集及全部排除原因；若主交集为空，按已冻结的程序支持数、同来源联系数、样本一致性、绝对差异效应和基因名排序，复核最多10个探索性/最近候选的文献与缺口，仍不得升级为主候选。')
p('完成与验收',head);p('完成质量、差异、损伤、巨噬、创新性各自全表；同基因、同状态交集与排除账本；原P/FDR、完整BH家族、NA与样本名单；关键值和BH独立核查；中文完整报告及统计附件交付下载目录。代码、原始数据、中间结果和日志保留服务器/工作目录。交集为空是本轮筛选完成，但不代表客户机制选题已经成功。')
p('计划与运行边界',head);p('复用冻结环境，可按需下载和在项目独立环境安装依赖，不升级共享环境。所有分析为同一发现队列上的探索性筛选，不作为独立验证。本计划的正式机器可读依据为PLAN_FROZEN.json，SHA256：1852715bf9c6a5a38aae4e0dae903833a40b7562df26902cef06989894073ac5。')
def foot(c,d):c.setFont('CJK',8);c.drawString(40,22,'四线联合筛选｜冻结计划v15');c.drawRightString(802,22,str(d.page))
SimpleDocTemplate(str(O/'four_axis_frozen_plan_v15.pdf'),pagesize=(842,595),leftMargin=40,rightMargin=40,topMargin=35,bottomMargin=40).build(story,onFirstPage=foot,onLaterPages=foot)
print('PLAN_PDF_CREATED')