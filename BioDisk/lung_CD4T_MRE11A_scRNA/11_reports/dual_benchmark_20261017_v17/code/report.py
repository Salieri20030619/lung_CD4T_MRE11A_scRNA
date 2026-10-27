from pathlib import Path
import pandas as pd,json,html,hashlib,zipfile,subprocess,datetime
from reportlab.platypus import SimpleDocTemplate,Paragraph,Table,TableStyle,Spacer,PageBreak
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
O=Path(r'C:\Users\31937\lung_dual_benchmark_20261017_v17');B=Path(r'C:\Users\31937\lung_three_cohort_20261017_v16');A=O/'assessment';D=Path(r'C:\Users\31937\Downloads\线粒体焦亡_双基准复评_v17_20261017');D.mkdir(exist_ok=True)
y=pd.read_csv(A/'90基因双基准总表.tsv',sep='\t');x=pd.read_csv(A/'逐状态双基准评估.tsv',sep='\t');pr=pd.read_csv(B/'final_tables/优先验证候选_5基因.tsv',sep='\t');summary=json.loads((A/'summary.json').read_text(encoding='utf-8'))
pdfmetrics.registerFont(TTFont('CN',r'C:\Windows\Fonts\ARIALUNI.TTF'))
styles={k:ParagraphStyle(k,fontName='CN',fontSize=size,leading=lead,wordWrap='CJK',spaceAfter=7,textColor=colors.HexColor('#164c5b') if k=='h' else colors.black) for k,size,lead in [('h',17,25),('b',10,16),('s',8,12),('c',8.3,12)]};story=[]
def p(t,k='b'):return Paragraph(html.escape(str(t)).replace('\n','<br/>'),styles[k])
def add(t,k='b'):story.append(p(t,k))
def page(t):story.append(PageBreak());add(t,'h')
def table(head,rows,widths):
 z=Table([[p(v,'c') for v in head]]+[[p(v,'c') for v in row] for row in rows],colWidths=widths,repeatRows=1);z.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#d7e9ec')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f2f6f7')]),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),5),('RIGHTPADDING',(0,0),(-1,-1),5),('TOPPADDING',(0,0),(-1,-1),4),('BOTTOMPADDING',(0,0),(-1,-1),4)]));story.append(z);story.append(Spacer(1,9))
add('巨噬细胞 - T细胞表型 - 损伤\n双基准二次评估','h');add('v17 | 2026-10-17 | 内部复审版；v16保持原样','s')
add('本轮按上一轮约定，把“候选入池”和“机制主候选”作为两个基准。线粒体损伤/内源性凋亡与焦亡分别记录，不把二者视作同义词，也不从横断面表达数据推定先后因果。')
table(['基准','二次评估结果','含义'],[['A：候选入池','52个计算支持入池；37个条件保留待补；1个本轮未达差异门槛或不可评估','52个具有分样本T细胞表达、名义差异、相容的巨噬入口和特定功能先验；不是52个已证实机制。'],['B：机制主候选','0个达到；90个均明确记录缺口','尚缺本项目目标T细胞自身损伤/死亡终点及同源同状态连续连接；不等于这些基因没有作用。']],[105,170,215])
add('范围为v16的55个交集基因与42个历史回溯基因的并集，共90个。复评344项名义差异对比；补齐4,680条基因-样本-状态表达记录，并检查22,100条不按排名截断的网络先验连接。')
add('主要变化：GZMB、TNFSF10转为“自身损伤/作用其他细胞的角色待澄清”；S100A8转为“来源待确认”。这3个从v16自动交集进入条件保留池，没有删除。新增历史4基因AHNAK、TTC39C、MDFIC、ST6GALNAC3仍保留，但不足以直接升级成目标损伤候选。')
add('现有IL1R2、DDIT4、SKIL、FOXP1仍属于基准A的优先验证对象，GSN为备选；五者均未达到基准B。创新性不作为排除基因的硬门槛。')
page('1  必要证据与参考证据分开执行')
table(['必要环节（A）','实际执行','不允许的替代'],[['目标T细胞表达/身份','冻结细胞身份；逐样本计数、阳性/多UMI、T标记及BG1记录；明显混合来源记待确认','共检出T标记不等于每个RNA都内源；缺证据不强判阴性'],['疾病/表型相关变化','沿用v16独立样本目标效应P<0.05；DE可作为发现阶段的分子表型代理','差异基因不等于功能性T细胞表型已改变'],['巨噬到相应T状态的入口','匹配发送端、接收状态、受体组分和候选；直接结合与下游预测分开','不能用单独配体表达替代候选连接；诱饵结合不等于激活信号'],['目标损伤功能','特定内源性凋亡/线粒体或焦亡成员及直接功能文献；分别记录外推边界','泛代谢、一般凋亡、p53富集本身不自动证明线粒体损伤'],['相容的细胞作用方向','不拼接不同巨噬范围/T状态的显著结果；区分T细胞自身死亡和杀死其他细胞','GZMB导致靶细胞焦亡不等于T细胞自身焦亡']],[108,207,175])
add('只作参考：FDR、|log2FC|≥0.5、跨队列一致性、伪时间、固定组成、组别调整、留一法、背景敏感性、NicheNet排名及创新性。没有把这些条件串成全部必过的显著性链。若参考检查揭示来源错误、状态不相容或单一样本主导等实质风险，必须影响分层。')
add('A1=有计算支持的候选资格；A2=条件保留待补，不是剔除；A3=本轮无合格名义差异或不可评估。三者是证据状态，不是新的统计检验。B要求在本项目内有T细胞自身损伤/死亡表型和连续作用支持；因果成立仍需相应干预证据。','s')
page('2  放宽排名后，缺口从“没进前250”转为具体问题')
add('读取v16所有正corrected-AUPR的合格配体，而非仅前30；对90基因的匹配名义对比检查全部非零靶权重。候选配体仍来自同组表达兼容的AllMac-受体入口。正权重只表示先验覆盖，不是显著调控，更不说明促进或抑制。')
add('共有81个复评基因有匹配网络连接；18个基因至少有一项对比的连接在旧参考排名范围外。排名靠后保留为弱线索，未因排名本身删除。原top30/250参考状态及每条边的配体、权重和排名均保留。')
table(['基因/问题','v17处理'],[['GSDMD','Core原队列最优靶排名405仍保留；原v14低UMI、背景RNA及严格过滤丢失样本的问题未解决，因此是来源待补，不再简单写成没有巨噬连接。'],['ST6GALNAC3','Core可找到较弱网络先验；不因排名581而排除，但目标损伤功能仍待证据。'],['IL1R2','受体结合与候选内源功能耦合分开；新增弱下游先验不能证明IL1R2传导激活信号。'],['GZMB、TNFSF10','从自动功能交集改为作用角色待澄清。需验证它们参与T细胞自身损伤，而不是仅作为对外效应分子。'],['S100A8','保留原始差异和网络信息，但髓系/背景来源未解决，不能直接认定为T细胞内源损伤候选。']],[105,385])
add('上述变化没有改变原始DE、P或FDR。低排名和缺少资源覆盖都不是生物学反证；相反，任何非零网络权重也不是自动通过机制基准的理由。')
page('3  v14新增候选的证据已重新分级')
table(['基因','新增核查','当前定位'],[['AHNAK','小鼠CD4钙流/NFAT/激活有直接研究；另有人癌细胞53BP1-p53/凋亡证据','T细胞表型优先补评。不能把两种细胞的证据拼成T细胞线粒体损伤。'],['TTC39C','A549/H1299敲低增加AnnexinV；能量/p53通路研究','一般凋亡线索；还缺T细胞目标损伤连接。'],['MDFIC','鸡成肌细胞干预及cleaved-caspase9结果','保留跨物种凋亡启动线索；尚缺哺乳T细胞及线粒体终点支持。'],['ST6GALNAC3','本次未核实到特定T细胞线粒体/焦亡直接证据','保留探索。检索到的同名lncRNA不能替代蛋白编码基因证据。']],[95,220,175])
add('v14的S07仍作为参考：AHNAK、ST6GALNAC3在第二小鼠队列两个时间点同向；MDFIC、TTC39C方向与原候选方向相反。伪时间/时间方向可以提高或降低优先级，但不能替代目标损伤功能。')
add('检索依据（直接原始研究，非仅综述）：','s')
for u in ['https://pmc.ncbi.nlm.nih.gov/articles/PMC2350190/','https://pmc.ncbi.nlm.nih.gov/articles/PMC8221568/','https://pmc.ncbi.nlm.nih.gov/articles/PMC9615393/','https://pmc.ncbi.nlm.nih.gov/articles/PMC6678156/','https://www.ncbi.nlm.nih.gov/nuccore/NM_152996']:add(u,'s')
add('一般凋亡与线粒体/焦亡分支的界限对所有候选一致执行。既有GO内源性凋亡成员可支持基准A功能先验，但也不能据此直接宣称已观察到损伤。')
page('4  机制基准未通过的原因是链条证据缺口')
table(['主线环节','目前能支持什么','尚缺什么'],[['巨噬细胞作用','同组配体-受体覆盖、下游模型与部分近端转录代理','来源特异功能作用、配体至候选的实际响应方向'],['T细胞表型改变','样本级候选表达、状态/组成及部分调控程序变化','候选导致特定功能表型改变的直接连接'],['T细胞自身损伤/死亡','内源性凋亡/焦亡相关转录程序与功能文献','同一目标T细胞中足够明确的线粒体或死亡执行终点'],['“然后发生”','部分独立时间队列的方向信息','同一作用链时间顺序；不能由静态通讯和伪时间推定因果']],[105,195,190])
add('v14并非没有表型线索：Core PyroExec分数差+0.0564，ApoExec差+0.00957，MitoROS差+0.0579；ActCD4 MitoROS差+0.0436。完整P、BH和样本范围随附。它们是转录代理，不是直接的膜电位、ROS实测、细胞色素c释放、gasdermin切割或死亡执行测量。')
add('3+3双侧精确置换最小P=0.10，因此不能用这些评分未达P<0.05来否定变化。基准B未通过并非FDR原因，而是目标终点与连续连接尚未建立。这里不规定必须使用某一种实验，只要求证据能识别所声称的目标表型。')
add('GZMB的直接原始研究说明它可在靶细胞内切割GSDME，恰好提示必须区分作用细胞。TNFSF10、IL10的输出作用与自身作用同样单列，当前待澄清不等于证明其绝无T细胞内在作用。','s')
add('参考：https://www.nature.com/articles/s41586-020-2071-9；死亡方式定义：https://pmc.ncbi.nlm.nih.gov/articles/PMC5864239/','s')
page('5  优先验证对象保留，但不越级命名')
rows=[]
for r in pr.itertuples():rows.append([r.gene,r.scope,f'{r.logFC:+.2f}',f'{r.PValue:.3g}',f'{r.FDR:.3g}','A1；B待补'])
table(['基因','重点范围','log2FC','P','FDR','双基准'],rows,[65,85,70,80,80,110])
add('DDIT4：外部直接研究最贴近线粒体损伤与ROS-TXNIP-NLRP3焦亡连接，优先检验该连接能否在T细胞出现。IL1R2：重点为保护性/反馈调节与ENO1，不是诱饵受体直接激活焦亡。SKIL：Treg应激-p53/内源性凋亡调节候选。FOXP1：跨队列下降较强，但巨噬来源不特异。GSN：保留为备选，关注胞内/分泌型、物种及对照异质性。')
add('AHNAK、TTC39C、MDFIC、ST6GALNAC3保留在A2补评池。优先补的是缺失环节，不是追加多个相关性检验直到显著。TNFRSF9、STAT1等有既往研究的基因仍可作比较/参照，不因“发表过”就排除。')
add('执行与验收：复用v16完整DE家族；补齐90基因分样本表达；取消配体/靶排名硬排除；完成344项状态对比和90项基因结论；逐项保存必要环节、参考证据、缺口和原始P/FDR。v14模块、近端检验、复合物及组成敏感性原表另存供审查。')
add('本报告是二次证据评估，不新增动物或实验结果，不把22个样本合并成单一疾病设计，不重标旧FDR，也不把两队列/多状态当作独立验证。原始矩阵与拟合对象仍在v16服务器目录。完整追溯见附件protocol_v17.json、assessment与tables。','s')
# 30 rows per appendix page, retain three-way statuses.
for k in range(0,90,30):
 page(f'附表  90基因双基准结论（{k+1}-{k+30}）')
 rows=[]
 for r in y.iloc[k:k+30].itertuples():rows.append([r.gene,r.A_status.split('_',1)[0],str(r.A_supported_scopes) if pd.notna(r.A_supported_scopes) else '待补/NA','未达到'])
 table(['基因','A状态','A支持范围（逐状态条件见附件）','B机制主候选'],rows,[95,55,245,95])
 add('A1：计算支持入池；A2：条件保留待补；A3：本轮无合格名义差异或不可评估。TLike时间证据不等同纯CD4疾病验证。完整缺口与参考检查见逐状态表，不以本简表代替。','s')
pdf=D/'v17_双基准二次评估报告.pdf'
def foot(c,d):c.setFont('CN',8);c.setFillColor(colors.HexColor('#526b75'));c.drawString(48,25,'v17 | 双基准二次评估 | 内部复审 | 2026-10-17');c.drawRightString(547,25,str(d.page))
SimpleDocTemplate(str(pdf),pagesize=A4,leftMargin=48,rightMargin=48,topMargin=42,bottomMargin=42).build(story,onFirstPage=foot,onLaterPages=foot)
(O/'delivery.json').write_text(json.dumps({'pdf':str(pdf),'directory':str(D)},ensure_ascii=False),encoding='utf-8');print(pdf)