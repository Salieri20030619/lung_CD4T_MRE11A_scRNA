from pathlib import Path
import json,html,os
import numpy as np,pandas as pd
from PIL import Image as PILImage
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import Paragraph,Table,TableStyle
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4,landscape
P=Path(r'C:\Users\31937\lung_death_priority_20261018_v3');A=P/'analysis';T=A/'tables';B=Path(r'C:\Users\31937\lung_death_priority_20261018_v1');BT=B/'analysis/tables'
O=Path(r'C:\Users\31937\Downloads\矽肺_凋亡焦亡优先_5加5证据补全_v3_20261018');O.mkdir(exist_ok=True)
if os.environ.get('REPORT_DRAFT')=='1':O=P/'draft_report';O.mkdir(exist_ok=True)
pdf=O/'矽肺_CD4凋亡焦亡优先_5加5生信证据补全_v3.pdf'
pdfmetrics.registerFont(TTFont('CN',r'C:\Windows\Fonts\ARIALUNI.TTF'))
styles={n:ParagraphStyle(n,fontName='CN',fontSize=s,leading=l,textColor=colors.HexColor('#213c48'),wordWrap='CJK') for n,s,l in [('body',10.5,17),('title',18,26),('sub',12.5,20),('small',8.5,13),('cell',8.5,13)]}
c=canvas.Canvas(str(pdf),pagesize=A4);c.setTitle('矽肺CD4凋亡焦亡优先：完整5+5生信证据补全 v3');c.setAuthor('研究分析')
page=0;pages=[];figure_entries=[];W,H=A4;y=0
def para(s,style='body'):
 global y
 p=Paragraph(s,styles[style]);_,h=p.wrap(W-96,H)
 if y-h<48:new(current_title+'（续）')
 p.drawOn(c,48,y-h);y-=h+11
def new(title,wide=False):
 global W,H,y,page,current_title
 current_title=title
 if page:c.showPage()
 W,H=A4;c.setPageSize((W,H));page+=1;y=H-48
 c.setStrokeColor(colors.HexColor('#287c8e'));c.setLineWidth(2);c.line(48,H-30,W-48,H-30)
 c.setFont('CN',8);c.setFillColor(colors.HexColor('#6b8089'));c.drawString(48,24,'凋亡/焦亡优先 · 独立路线 v3 · 2026-10-18');c.drawRightString(W-48,24,str(page))
 pages.append({'page':page,'title':title,'orientation':'portrait'});para(title,'title')
def table(headers,rows,widths=None):
 global y
 if widths is None:widths=[(W-96)/len(headers)]*len(headers)
 else:widths=[w*(W-96)/sum(widths) for w in widths]
 data=[[Paragraph(html.escape(str(v)),styles['cell']) for v in rr] for rr in [headers]+rows]
 tb=Table(data,colWidths=widths,repeatRows=1)
 tb.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#dcecef')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f2f6f7')]),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),5),('RIGHTPADDING',(0,0),(-1,-1),5),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#8aa9b2'))]))
 _,h=tb.wrap(W-96,H)
 if y-h<48:
  if len(rows)>1:
   split=max(1,len(rows)//2);table(headers,rows[:split],widths);new(current_title+'（续）');table(headers,rows[split:],widths);return
  new(current_title+'（续）')
 tb.drawOn(c,48,y-h);y-=h+13
def fmt(v):
 if pd.isna(v):return 'NA'
 return f'{v:.3g}' if abs(float(v))<.001 else f'{float(v):.4f}'
def rd(n):return pd.read_csv(T/n,sep='\t')
def br(n):return pd.read_csv(BT/n,sep='\t')
manifest={x['id']:x for f in ['manifest_core.json','manifest_comm.json','manifest_joint_network.json','manifest_ORA.json'] for x in json.loads((A/'figures'/f).read_text(encoding='utf-8'))}
old={x['id']:x for x in json.loads((B/'figure_manifest.json').read_text(encoding='utf-8'))}
def fig(fid,title,insight,caption,baseline=False):
 global y
 f=(old if baseline else manifest)[fid];path=Path(f['file']) if baseline else A/'figures'/Path(f['file']).name
 im=PILImage.open(path);iw,ih=im.size;wide=iw/ih>1.45
 new('图'+str(len(figure_entries)+1)+'｜'+title,wide);para(insight)
 figure_entries.append({'figure':len(figure_entries)+1,'page':page,'title':title,'id':fid,'source':('v1冻结筛选表及图'+fid if baseline else f.get('source',''))})
 if fid.startswith('N08_'):
  sc=fid.split('_',1)[1];dd=rd('joint_GSEA_all.tsv');dd=dd[(dd.scope==sc)&(dd.database=='KEGG')].sort_values(['pval','pathway']).head(5)
  caption+=' 图例：'+'；'.join(html.escape(r.pathway+' '+r.description.replace(' - Mus musculus (house mouse)','')) for r in dd.itertuples())+'。'

 cp=Paragraph(caption,styles['small']);_,ch=cp.wrap(W-96,H)
 scale=min((W-96)/iw,(y-60-ch-12)/ih);wi,he=iw*scale,ih*scale
 c.drawImage(str(path),48+(W-96-wi)/2,y-he,width=wi,height=he);y-=he+10;para(caption,'small')
def assoc(g,sc,mod,var='Raw',model='JointStratum'):
 return assocdf[(assocdf.cohort=='Joint')&(assocdf.variant==var)&(assocdf.model==model)&(assocdf.gene==g)&(assocdf.scope==sc)&(assocdf.program==mod)].iloc[0]
def atext(r):return f'ρ={fmt(r.rho)}，P={fmt(r.P)}，FDR={fmt(r.FDR)}'
assocdf=rd('all_associations.tsv');summary=rd('candidate_evidence_summary.tsv');base=br('primary_death_candidates.tsv');inventory=rd('completion_inventory.tsv');inventory=inventory[~inventory.evidence.astype(str).str.startswith('TF_')].reset_index(drop=True)
jqa=json.loads((A/'joint_network/QA_complete.json').read_text(encoding='utf-8'))
assert jqa['status']=='success'
new('完整5+5分析补齐后，三个主候选呈现不同证据方向')
para('矽肺CD4凋亡/焦亡优先研究 · 独立路线 v3','sub')
para('本报告以v1的新路线筛选结果为基准，套用既有分析框架，补齐两个疾病队列的完整5个对照与5个矽尘样本。主候选仍为HIF1A、IL2RA、ARL6IP1；次级时间分支单独评价。')
para('IL2RA在Treg内有稳定的疾病上调、代谢程序关联及联合新网络的NFKB2名义关联，但其证据强度对签名构建和供者剔除敏感；HIF1A在Naive内有焦亡评分的名义关联；ARL6IP1在Naive内有线粒体ROS评分的名义关联，后两者对评分成员构成更敏感。三者在矽尘肺中的定位仍为探索性候选。新增外部真实敲除分析支持IL2RA维持外周Treg存活；多来源人肺Treg的代谢正相关获FDR支持，但同来源疾病上调未复现。')
para('补齐后的主要边界同样明确：主交集没有焦亡核心执行成员；联合新网络的IL2RA-NFKB2关联在疾病调整后未保留；巨噬IL2未达到配体检出门槛。当前结果支持进一步分解候选功能，不支持“巨噬刺激导致CD4焦亡”的机制结论。新增公开真实敲除数据用于区分候选的死亡调控方向，详见外部干预章节。')
para(f'交付状态：完整5+5分析及All/BG1三种子重新建网已完成。原52项核查、成员敏感性20个家族以及联合新网络{jqa["checks"]}项核查通过；主候选关键关联及新网络关联均完成独立数值复算。报告保留完整阳性、主要阴性及敏感性结果。')
if (P/'supplement_QA.json').exists():
 v3qa=json.loads((P/'supplement_QA.json').read_text(encoding='utf-8'))
 para(f'v3新增数值核查：{v3qa["checks"]}项通过，包括42个外部程序效应的独立实现复算；报告全部使用竖版。','small')
para('性质：内部生信证据补全报告。原v1及既往路线保留；本报告不自动构成客户发送或实验因果验证。','small')
new('两个疾病队列组成5+5，独立时间队列回答时间变化')
table(['队列','独立样本','分析角色'],[['GSE314943','3对照＋3矽尘，D140','主队列；CD4与巨噬发送端'],['GSE183682','D7和D56各1＋1，共4只','联合补充；本轮新增全肺巨噬注释'],['疾病联合','5对照＋5矽尘，共10只','按队列/时间分层的疾病效应与关联'],['GSE274372','D0/D7/D28/D56各3只，共12只','独立真实时间；TLike，单独建模']], [90,145,264])
para('疾病联合CD4细胞共3,724个。统计样本量为10只独立小鼠，细胞数仅用于描述覆盖及细胞级评分聚合。D7/D56的1＋1单独只报告方向，不生成独立推断性P值。')
para('新路线保持Naive/ActCD4与Treg两条独立状态轨迹；不将Treg状态排序解释为Naive分化至Treg。时间队列TLike也不等同于疾病队列的纯Naive、ActCD4或Treg。')
para('主分析保留原始P＜0.05的名义筛选规则，所有多重检验保留对应同类假设家族的FDR。FDR未通过不等同于没有方向或效应；只有FDR实际通过时才称为校正后显著。')
for start,end,title in [(0,8,'CD4与当前TF分析覆盖完整联合样本'),(8,len(inventory),'发送端、联合建网与补充敏感性均已完成')]:
 new(title)
 rows=[]
 for r in inventory.iloc[start:end].itertuples():rows.append([r.id,r.analysis,r.sample_scope,r.status])
 table(['编号','分析','样本范围','完成情况'],rows,[32,190,105,172])
 para('“复用”表示原分析已包含完整5+5且本轮未改动统计值；“新增”表示针对新候选或原未覆盖的联合范围补算。完整表名与证据位置见补充材料索引。','small')
 if start==0:para('TF结果全部来自当前完整10样本的All/BG1三种子联合重建网络。')
fig('F09','v1主筛选仍得到3个凋亡相关候选','状态轨迹P＜0.05 ∩ D28对D7 P＜0.05 ∩ 凋亡/焦亡GO的主交集保持不变。','继承v1筛选基准图。焦亡分支及核心执行成员交集为空；本轮功能与通讯分析用于补证据，不反向修改筛选门槛。',True)
fig('F03','联合状态轨迹保留两个独立分析范围','Naive/ActCD4的连续状态排序与Treg内部梯度分别建模。','继承v1联合10样本轨迹。Harmony仅按数据集处理几何；疾病效应统计使用独立小鼠及队列/时间分层。',True)
fig('F04','ActCD4凋亡程序沿状态排序上升','v1的样本级斜率分析支持ActCD4凋亡程序趋势；焦亡程序未获得同样支持。',old['F04']['caption']+' 这些为v1筛选程序；v2候选关联另行使用去除7个候选的程序评分。',True)
fig('F05','Treg自噬和线粒体自噬程序随状态排序上升','v1轨迹中的自噬相关程序为新候选提供功能背景。',old['F05']['caption'],True)
fig('F10','状态轨迹与真实时间的变化方向并不完全一致','HIF1A与IL2RA的状态/时间方向较一致，ARL6IP1保留方向相反这一关键结果。',old['F10']['caption'],True)
new('主候选保留完整筛选效应、P值和FDR')
rows=[]
for r in base.itertuples():rows.append([r.gene+' / '+r.scope,fmt(r.mean_slope),fmt(r.P)+' / '+fmt(r.FDR),fmt(r.time_logFC),fmt(r.time_P)+' / '+fmt(r.time_FDR)])
table(['候选/状态','轨迹斜率','轨迹P / FDR','D28-D7 log2FC','时间P / FDR'],rows,[95,72,116,93,123])
para('三个候选均为名义显著的交集成员，未被重新标为FDR显著。HIF1A为内源性凋亡注释；IL2RA、ARL6IP1为凋亡调控/相关注释；三者均不是本筛选中的焦亡核心执行基因。')
para('状态轨迹的疾病特异斜率检验：HIF1A P=0.875，IL2RA P=0.100，ARL6IP1 P=0.475；当前数据不支持其状态趋势具有明确的疾病特异性。')
fig('F11','IL2RA在四种轨迹设置下保留名义趋势','HIF1A和ARL6IP1对Naive轨迹根节点更敏感；该结果贯穿后续候选定位。',old['F11']['caption'],True)
fig('N01','IL2RA在Treg内呈现名义疾病上调','完整5+5的亚群内差异支持IL2RA上调；HIF1A和ARL6IP1在原提名Naive亚群中未达到疾病差异P＜0.05。','点为独立小鼠。圆点/方点区分队列；差异来自原完整基因家族的edgeR队列/时间分层模型，非图中均值的单独检验。')
fig('N15','IL2RA疾病上调方向不由单个样本决定','IL2RA的10次逐样本剔除拟合均保留上调；HIF1A为9/10同向，ARL6IP1为8/10同向。','区间表示逐一剔除后的效应范围，不是置信区间。每次拟合n=9，FDR保留该次完整基因家族；重复剔除不是新增独立验证。')
fig('N02','扩展程序的联合疾病效应呈现亚群差异','本轮统一比较14类程序，并在候选关联用评分中移除全部7个新路线候选。','颜色为Silica减Saline的分层效应；*表示原始P＜0.05，不代表FDR通过。各程序只在自身尺度内解释，不能跨程序按绝对分数排序。')
fig('N03','主候选分别关联状态、代谢与损伤程序','HIF1A-焦亡评分、ARL6IP1-线粒体ROS、IL2RA-氧化磷酸化/糖酵解得到名义关联；完整主要阴性结果同图保留。','相关为队列/时间调整后的秩残差相关，n=10，2,880种分层置换；*为原始P＜0.05。完整P/FDR见all_associations.tsv。评分已去除全部7个候选。')
fig('N04','功能关联对疾病调整和背景过滤的稳定性不同','HIF1A-焦亡评分及ARL6IP1-ROS在BG1下保留方向；IL2RA的代谢关联在进一步调整疾病组别后减弱。','BG1为既有背景分数≤1%的接收细胞范围。条件调整使用36种块内置换；P值离散性与效应变化需同时解释。敏感性分析的显著性不能替代主分析。')
exec((P/'module_member_report_section.py').read_text(encoding='utf-8'))
for gene,scope,mods,title,explain in [
 ('HIF1A','Naive',['TReg','PyroExec','ApoExec','MitoROS'],'HIF1A连接Naive状态与焦亡评分，但未显示Naive疾病差异','HIF1A-PyroExec在BG1保留名义关联，条件调整后P未过0.05。PyroExec逐成员移除仅2/5保留名义关联，不能将RNA关联解释为HIF1A驱动焦亡。'),
 ('IL2RA','Treg',['OXPHOS','Glycolysis','ApoExec','PyroExec'],'IL2RA的Treg疾病表达与代谢关联最为明确','IL2RA-ApoExec主分析P=0.0531，未达到P＜0.05；BG1下增强，但疾病条件调整后仍不足。IL2RA-PyroExec主分析不支持关联；糖酵解在117次逐成员移除后均保留名义关联。'),
 ('ARL6IP1','Naive',['MitoROS','ApoExec','PyroExec','Autophagy'],'ARL6IP1提供Naive线粒体ROS关联线索','ARL6IP1的Naive疾病差异不支持，真实时间方向与状态斜率相反。ROS关联在BG1及其条件调整中方向保留，逐成员移除为7/12保留名义关联，适合作为功能分解线索。')]:
 new(title);r=summary[(summary.gene==gene)&(summary.scope==scope)].iloc[0]
 para(f'疾病效应：log2FC={fmt(r.DE_logFC)}，P={fmt(r.DE_P)}，FDR={fmt(r.DE_FDR)}；逐样本剔除{int(r.LOO_same_direction)}/10同向。')
 rows=[]
 for mod in mods:
  v=assoc(gene,scope,mod);rows.append([mod,fmt(v.rho),fmt(v.P),fmt(v.FDR),str(int(round(v.loo_same_sign*v.loo_valid)))+'/'+str(int(v.loo_valid))])
 table(['程序','ρ','原始P','FDR','剔除同向'],rows,[125,75,100,100,99]);para(explain)
 para('主分析n=10；关联为探索性结果，使用预设检验家族。疾病差异与程序关联是不同假设，不互相替代。','small')
fig('ORA01','Treg上下调基因均有FDR支持的GO富集',
    'Treg上调/下调GO分别有117/11项FDR＜0.05；Naive上调23项，ActCD4上调6项。另两亚群下调分支保留名义条目。',
    '输入为相应5+5亚群疾病差异中原始P＜0.05的上/下调基因；每亚群/方向选择前2个名义条目并列展示。色值为-log10 FDR，黑边表示FDR＜0.05，点大小为命中基因数；Core结果保留在完整表中。')
fig('ORA02','Treg的KEGG上下调分支均保留校正后富集',
    'Treg上调/下调分别有1/2项KEGG条目通过FDR；Naive与ActCD4的KEGG ORA保留名义信号。不同条目的功能含义按名称及成员解释。',
    '同一预定义数据库、亚群和方向内进行完整BH校正。图中只绘制原始P＜0.05的条目；空白表示该状态/方向未达到名义阈值或不在可评估集合中。完整候选、背景和成员见joint_ORA_all.tsv。')
fig('N06','联合富集显示死亡与代谢信号的亚群差异','整体CD4富集、亚群内富集与固定执行程序评分提供不同层面的证据，应并列呈现。','展示预先围绕死亡、线粒体及代谢问题关注的通路；全部GO_BP/KEGG结果保留。点大小对应原始P，颜色为NES，黑边为FDR＜0.05。')
for sc in ['Naive','ActCD4','Treg']:
 fig('N08_'+sc,sc+'的KEGG富集保留完整排序曲线','按原始P选取该亚群KEGG前5项用于展示，正文不将未过FDR的曲线称为校正后显著。','上：真实加权running ES；中：命中位置；下：全基因sign(logFC)×sqrt(F)排序。曲线ES已与fgsea输出核对。所有通路在各自完整家族内校正。')
fig('N07','Core焦亡炎症富集在组成调整后不再获得支持','调整ActCD4与Treg比例后，焦亡炎症及部分糖酵解/淋巴细胞凋亡条目信号减弱；KEGG凋亡仍保留名义信号。','Core原模型：按队列/时间分层。敏感性模型再纳入两项亚群比例，n=10、残差自由度4、条件数约103。KEGG凋亡调整后P=0.003025、FDR=0.073858；不据此断言所有变化均由组成引起。')
exec((P/'joint_network_report_section.py').read_text(encoding='utf-8'))
fig('N09','第二队列全肺图谱补齐了发送端识别','四个样本经质控和双细胞筛查后保留41,604个细胞，17个簇；两个巨噬富集簇用于联合发送端。','全肺UMAP仅用于注释与覆盖展示。MixedAmbient、Mixed及DC/MonoAPC单独标记，不并入主巨噬范围；细胞级marker排序用于描述，不作为独立样本推断。')
fig('N10','两个巨噬富集簇与单核/DC范围具有可辨别的标志差异','簇2具有AM样标志，簇4具有巨噬及单核来源特征；簇5和9分别保留为单核及DC/单核抗原呈递范围。','色值为四个样本的平均log表达，点大小为阳性细胞比例的样本均值。四样本巨噬富集范围为1,073、2,093、1,452、2,633个细胞。')
fig('N11','完整10个样本均具备巨噬发送端覆盖','主AllMac范围与相同标志规则的MacStrict范围分别计算，未拼接两范围的显著结果。','MacStrict要求单细胞至少表达2个C1QA/C1QB/C1QC/ADGRE1/MERTK/MSR1/MARCO标志；两队列使用同一规则。发送端独立样本数始终为10。')
new('新增发送端质控支持分析，同时保留自动双细胞筛查边界')
qc=rd('macrophage_annotation_QC.tsv')
table(['样本','巨噬细胞','分数＞0.25','≥2 T标志','≥2 B标志'],[[r.sample,int(r.cells),int(r.doublet_score_gt025),int(r.Tmarkers_ge2),int(r.Bmarkers_ge2)] for r in qc.itertuples()],[105,90,105,100,99])
para('第二队列两个矽尘样本的Scrublet自动阈值较高，预测双细胞数为0，不能解释为真实双细胞不存在。巨噬范围内分数＞0.25及混合淋巴标志细胞占比较低；共同标志敏感性与完整来源竞争结果用于保留这一风险。')
para('发送端全基因联合edgeR完成11,559个基因的疾病效应检验，保留完整FDR家族。两队列巨噬亚型构成不同，因此发送端效应同时反映该范围内的表达及组成差异。')
fig('N12','联合通讯表达代理获得名义疾病变化','AllMac有205个可评估“受体状态-配体受体”组合，其中85个原始P＜0.05；MacStrict为209个、89个。两范围均无FDR＜0.05组合。','图示按原始P排序的代表组合；空白表示该状态不满足表达资格。评分是配体及完整受体亚基表达的几何均值代理，不是CellChat概率、分泌量或受体激活测量。')
fig('N13','联合NicheNet提供响应配体排序，而非因果检验','完整发送端和接收端分别对疾病差异及状态轨迹响应建模，两个发送范围合计28个分析组合。','图为AllMac疾病上调/下调响应的前6配体。corrected AUPR是先验网络对响应基因集合的预测指标，不是生物学P值。完整候选权重、排名和响应归属均保留。')
new('上游关联以负向名义信号为主，主候选机制连接仍未闭合')
u=rd('joint_upstream_associations.tsv');z=u[(u.variant=='Raw')&(u.model=='JointStratum')&(u.sender=='AllMac')&u.kind.eq('ligand_candidate_primary')&(u.P<.05)].sort_values('P')
table(['候选/状态','来源配体','ρ','P','FDR'],[[r.gene+' / '+r.scope,r.ligand,fmt(r.rho),fmt(r.P),fmt(r.FDR)] for r in z.itertuples()],[125,110,85,90,89])
para('这些连接均来自相同10只小鼠、相同接收亚群；方向为负，且均未过FDR。不能将细胞接触类配体或负向关联直接解释为促焦亡刺激。')
para('HIF1A与IL2RA分别有24和27个AllMac先验支持配体进入关联评价；ARL6IP1在预设NicheNet top250先验范围中无可评估路线。先验缺失表示当前资源未覆盖，不表示生物学上不存在上游调控。')
fig('N14','IL2RA受体表达并不支持巨噬IL2来源链','Treg中IL2RA/IL2RB/IL2RG的表达资格为10/10样本；巨噬IL2在10个样本均未达到预设10%细胞检出门槛。','竞争来源按队列分别作描述性样本均值。低于检出门槛不等于绝对零表达；不同Broad范围不用于推断分泌贡献。')
new('次级分支的CD44-Treg证据较强，但保持次级定位')
sec=summary[summary.branch=='secondary']
table(['候选/状态','疾病log2FC','P / FDR','功能/联合新TF/上游名义项数'],[[r.gene+' / '+r.scope,fmt(r.DE_logFC),fmt(r.DE_P)+' / '+fmt(r.DE_FDR),f'{int(r.program_nominal)} / {int(r.jointTF_nominal)} / {int(r.upstream_nominal)}'] for r in sec.itertuples()],[130,90,135,144])
para('CD44-Treg疾病上调：log2FC=1.0318、P=0.001427、FDR=0.038544，10次逐样本剔除均同向；并有功能及上游名义关联。其进入依据仍是次级时间总效应分支，不替代主D28对D7筛选。')
para('FOXP1-ActCD4疾病下降为名义信号；PRKCQ与TRPS1保留状态/程序线索及阴性环节。表中计数仅描述已测试结果，不构成“显著项越多越优”的综合评分。')
new('补齐后可明确三条功能线索，而非单一焦亡机制')
table(['主候选','当前最清楚的证据','仍缺的关键环节'],[['IL2RA / Treg','疾病上调稳定；糖酵解逐成员移除117/117保留；NFKB2有名义关联','主分析焦亡关联阴性；巨噬IL2不满足来源门槛；联合NFKB2名义关联在疾病调整后未保留'],['HIF1A / Naive','TReg评分关联FDR通过；PyroExec名义关联，但逐成员移除仅2/5保留','Naive疾病差异不支持；换根敏感；联合网络主分析同亚群TF连接未达到名义显著'],['ARL6IP1 / Naive','MitoROS名义关联；BG1条件调整保留方向，逐成员移除7/12保留名义关联','疾病差异不支持；时间与轨迹方向相反；当前NicheNet先验未覆盖']],[105,185,209])
para('本轮不改变主候选集合。IL2RA最适合围绕Treg状态与代谢开展进一步生信分解；HIF1A适合核查焦亡评分由哪些成员驱动；ARL6IP1适合核查ROS相关性与背景敏感性。')
para('次级CD44-Treg可作为独立备选支线，并明确其时间筛选依据。不得把候选、程序RNA评分和通讯预测合并描述为已验证的死亡执行。')
exec((P/'supplement_v3_report.py').read_text(encoding='utf-8'))
new('讨论：证据强度取决于样本、亚群与测量层级')
para('独立生物学单位少于10的外部转录分析均按探索性证据解释，包括外周Treg的4＋4小鼠、人HIF1A的4名配对供者及小鼠HIF1A的小规模培养比较；方向、效应量和完整P/FDR同时保留。原论文的蛋白与死亡表型实验作为另一个来源层级列示。')
para('跨物种程序使用同名符号投射，实际成员覆盖不同，不直接比较绝对评分。外部IL2RA的脾与胸腺实验在删除设计、组织和批次上不同，分开分析所得的方向差异不能单独归因于组织。外部IL2RA真实敲除提供了Treg存活调控的实验背景，HIF1A真实干预强化了代谢调控方向；二者不等于在矽尘肺内完成候选-指定通路-死亡执行的验证。人肺多来源程序关联与同来源病例差异回答不同问题：前者支持IL2RA代谢关联，后者未复现疾病上调；关键成员与BAL覆盖的未复现或不足亦完整保留。')
para('两个疾病队列在时间、技术和病程上不同。分层模型估计跨层共同疾病效应；第二队列每个时间点只有1＋1，其单独方向不能提供独立推断性验证。完整n=10仍不足以支持精细交互项或稳健因果推断。')
para('候选来自当前数据的轨迹与时间交集，随后同数据的关联是探索性补证据，不是独立确认。名义P、FDR、效应方向和样本一致性均须并列解释。')
para('Core组成调整只有4个残差自由度且条件数约103；调整后信号变化提示组成敏感性，不能证明全部差异仅由组成导致。BG1是背景负荷筛选，不是去污染；严格巨噬范围也不能排除所有双细胞。')
para('联合新建网络、NicheNet和配体受体表达代理提供转录层支持。它们不直接测量TF蛋白活化、线粒体损伤、自噬通量、焦亡蛋白切割或分泌量。当前没有足够依据将这些候选写成可用临床靶点。')
para('主要阴性结果保留：焦亡核心执行主交集为空；三个主候选疾病特异轨迹未显著；IL2RA-NFKB2关联在疾病调整后未保留；巨噬IL2来源门槛未通过。')
new('方法：输入、质控与功能评分')
para('疾病接收端使用已冻结GSE314943和GSE183682的Core原始计数及Naive/ActCD4/Treg注释；时间队列使用GSE274372的494个TLike细胞。输入路径和SHA256保存在manifest。v1轨迹与时间筛选原值保持不变。')
para('第二队列全肺原始矩阵按样本进行基因数≥200、UMI≥300及线粒体比例＜20%过滤，再以Scrublet expected_doublet_rate=0.06、种子20261018筛查。归一化log1p后选2,500 HVG、30维PCA，按样本Harmony、15邻居、Leiden resolution=0.6及UMAP；标志与样本覆盖确定注释。')
para('v2样本程序评分为完整文库CP10k的log1p后成员均值；使用三队列共同测得成员，全部7个新路线候选从所有程序移除。表型-损伤关联还移除两程序共享成员。保留Raw、BG1及Core固定1/3状态权重范围。')
para('独立时间队列保持12样本主分析和≥10细胞样本敏感性；程序/候选时间关联单独报告时间合并与day调整模型。v1用于筛选的程序成员与v2防自相关程序不混写。')
para('基于候选结果追加五项程序关联的逐成员移除敏感性，Raw/BG1及两个模型共660项检验、20个家族；移除成员后的评分仍为剩余成员均值，主筛选与原主分析结果保持不变。')
new('方法：统计单位、分层模型与多重检验')
para('疾病差异采用edgeR robust quasi-likelihood，完整样本伪bulk计数、完整文库分母、TMM及队列/时间分层设计；目标系数为Silica对Saline。原完整基因家族FDR保留。逐样本剔除重新拟合40个模型；Core组成敏感性增加ActCD4和Treg比例。')
para('程序、TF及LR表达代理的疾病效应采用分层线性效应和80种保持层内组别数的精确分配。候选关联为秩残差相关：主模型调整队列/时间，2,880种层内排列；条件调整敏感性为36种层内且组内排列。表中P_model_t仅为另列的近似模型值，不替代主置换P。')
para('FDR按预定义同类家族分别计算：基因差异、GO_BP/KEGG、程序效应、候选-表型/损伤/辅助程序、表型-损伤、TF-候选、TF-程序、LR效应、配体-候选及配体-TF分别分层，主/次候选、范围、变体和模型均保留。')
para('原始P＜0.05为名义筛选硬标准；FDR＜0.05才称校正后显著。完整家族、计划/可评估检验数、样本列表和缺失/不可评估状态均保留。P=0.05不计入P＜0.05。')
new('方法：富集、TF与通讯')
para('GO_BP使用org.Mm.eg.db的GOALL；KEGG使用项目冻结REST映射。ORA按原始P＜0.05的上调/下调基因分别计算超几何检验，背景为被测试且数据库可注释基因。GSEA使用sign(logFC)×sqrt(F)全基因排序和fgseaMultilevel，集合15-500基因。数值边界/精度条目定向提高采样并重算完整家族BH。')
para('TF使用当前联合重建网络的motif支持稳定regulon；去除候选、全部程序成员及TF自身，至少5个目标后进行AUCell，auc_threshold=0.05、seed=17、noweights=True。细胞活性聚合为独立样本/状态均值。')
para('联合新建网络：All/BG1分别使用全部10样本，每样本训练上限250细胞，种子17/29/43。共同基因映射至两套mm10 cisTarget v10窗口，检出阈值max(10,1%细胞)，全库CP10k log1p；GRNBoost2、motif剪枝后取至少2/3种子重复边，去重叠并以全部保留细胞重新评分。')
para('通讯使用冻结CellChat蛋白配体-受体及完整复合物亚基定义。每端亚基最小检出比例≥10%的样本至少3/10，且两队列均至少1个；合格组合计算log1p(sqrt(配体几何均值CP10k×受体几何均值CP10k))。')
para('NicheNet使用冻结mouse ligand-target矩阵；分别以联合疾病差异或状态轨迹P＜0.05的上/下基因为响应，corrected AUPR＞0且排名≤30的配体提供候选先验，候选权重＞0且target rank≤250进入上游关联。两个发送范围独立保留。')
new('复核与复现：数值、输入及主要限制均可追溯')
q=json.loads((A/'logs/numerical_QA_complete.json').read_text())
para(f'v2基础分析的数值与完整性核查共{q["checks"]}项通过，包含{q["independent_correlations"]}项主候选关键关联的独立P/ρ复算、完整BH家族复算、样本数核对、候选/TF去重叠、原始输入SHA256及发送端覆盖。')
para('完成标记与进程退出状态共同判断分析成功。通讯全肺来源导出采用等价的稀疏矩阵投影加速；先前输入就绪/布尔读取调试日志及被替换的慢导出记录保留于工作目录。')
para(f'联合重新建网另通过{jqa["checks"]}项核查，{jqa["independent_primary_TF_tests"]}项主候选TF关联独立复算；评分成员敏感性20个家族复算通过。')
if (P/'supplement_QA.json').exists():
 v3qa=json.loads((P/'supplement_QA.json').read_text(encoding='utf-8'))
 nv=len(pd.read_csv(P/'tables/independent_network_LOO_QA.tsv',sep='\t')) if (P/'tables/independent_network_LOO_QA.tsv').exists() else 0
 para(f'v3追加{v3qa["checks"]}项数值与完整性核查，外部程序独立重拟合42项；逐小鼠网络的IL2RA-NFKB2关键关联另完成{nv}项独立ρ与精确P复算。独立实现复算属于技术验收，不增加生物学样本量。')
if (P/'tables/GSEA_precision_audit.tsv').exists():
 pr=pd.read_csv(P/'tables/GSEA_precision_audit.tsv',sep='\t')
 para(f'新增GSEA统一提高精度{len(pr)}项；相对初次估计，{int(((pr.P_before>=.05)&(pr.P_after<.05)).sum())}项进入、{int(((pr.P_before<.05)&(pr.P_after>=.05)).sum())}项退出名义阈值。全部前后值与原完整家族FDR可追溯。')
para('复核范围为本次输入映射、统计结果、图表及PDF交付完整性；未进行独立外部AI评审，未验证湿实验因果、死亡蛋白执行或临床外推。')
table(['资源/方法','可追溯来源'],[['队列与注释','GSE314943、GSE183682、GSE274372；inputs/manifest.json'],['筛选基准','独立路线v1及其冻结协议、主/次候选完整表'],['基因集与LR先验','org.Mm.eg.db/GOALL、冻结KEGG、CellChatDB.mouse'],['TF与NicheNet','联合17/29/43的All/BG1网络；冻结NicheNet矩阵'],['分析实现','code/、protocol.json、R session与各完成标记']],[140,359])
new('补充材料索引与交付位置')
rows=[['样本/候选','frozen_candidate_pairs.tsv；candidate_expression.tsv；sample_scores.tsv'],['差异/稳定性','candidate_existing_DE.tsv；candidate_DE_leave_one_donor.tsv；joint_sender_DE_all.tsv'],['功能/富集','all_associations.tsv；program_effects.tsv；joint_ORA_all.tsv；joint_GSEA_all.tsv；module_member_sensitivity*.tsv'],['TF','joint_network_TF_*.tsv（联合新网络）；joint_network_candidate_edges.tsv'],['通讯/来源','joint_LR_eligibility.tsv；joint_LR_disease_effects.tsv；joint_NicheNet_*；whole_lung_source_expression.tsv.gz'],['上游/整合','joint_upstream_associations.tsv；joint_network_upstream_TF_associations.tsv；candidate_evidence_summary.tsv；completion_inventory.tsv'],['复核','QA_checks.tsv；joint_network_QA_checks.tsv；module_member_QA.tsv；inputs/manifest.json及additional_provenance.json']]
table(['类别','主要统计表'],rows,[95,404])
para('完整竖版报告保存到Windows下载目录；本轮新增材料归档于v3，v1/v2保持不变。','small')
para('服务器：lung_CD4T_MRE11A_scRNA/11_reports/death_priority_20261018_v2。工作目录：C:/Users/31937/lung_death_priority_20261018_v3。','small')
for start in range(0,len(figure_entries),12):
 new('图表索引：图号、页码及数值来源'+('（续）' if start else ''))
 rows=[[v['figure'],v['page'],v['title'],v['source']] for v in figure_entries[start:start+12]]
 table(['图','页','内容','来源'],rows,[28,28,205,238])
 para('图号在本报告中连续编号。v1基准图仅承接冻结结果；新增图对应v2基准及v3追加统计表。完整数值表位于复核归档analysis/tables。','small')
(P/'figure_index.json').write_text(json.dumps(figure_entries,ensure_ascii=False,indent=2),encoding='utf-8')
c.save();(P/'report_pages.json').write_text(json.dumps(pages,ensure_ascii=False,indent=2),encoding='utf-8');(P/'report_delivery.json').write_text(json.dumps({'pdf':str(pdf),'pages':page,'version':'v3'},ensure_ascii=False,indent=2),encoding='utf-8')
print('PDF',page,str(pdf))
