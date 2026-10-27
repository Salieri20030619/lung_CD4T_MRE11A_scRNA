from pathlib import Path
import json,hashlib,shutil,textwrap
import numpy as np,pandas as pd
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A3,landscape
from reportlab.lib import colors
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import Paragraph,Table,TableStyle
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.utils import ImageReader
from PIL import Image
P=Path('technical_delivery');D=P/'data';S=P/'source_tables';F=P/'figures'
(P/'output/pdf').mkdir(parents=True,exist_ok=True)
def rd(n,where=D):return pd.read_csv(where/n,sep='\t')
def criterion(p,q=np.nan):
 if pd.isna(p):return 'NO_VALID_P'
 if p>=.05:return 'NOT_NOMINAL'
 return 'FDR_PASS' if pd.notna(q) and q<.05 else 'NOMINAL_ONLY'
records=[]
def addrows(df,kind,source,pcol,ecol,qcol='FDR',cohort=None):
 for _,r in df.iterrows():
  p=r.get(pcol,np.nan);q=r.get(qcol,np.nan)
  records.append(dict(evidence_type=kind,cohort=cohort or r.get('species','Human'),scope=r.get('scope',r.get('subtype','')),tf=r.get('tf',r.get('tf_gene','')),feature=r.get('feature',r.get('ligand','')),filter=r.get('filter','All'),seed=r.get('seed',''),variant=r.get('variant',''),sender=r.get('sender',''),sender_filter=r.get('sender_filter',''),sample_set=r.get('sample_set',''),model=r.get('model',''),n=r.get('n',np.nan),n_control=r.get('n_control',np.nan),n_case=r.get('n_case',np.nan),estimate=r.get(ecol,np.nan),p=p,FDR=q,df=r.get('df',np.nan),reason=r.get('reason',''),acceptance=criterion(p,q),source=source))
addrows(rd('existing_RNA_DE_focused.tsv'),'RNA_DE','data/existing_RNA_DE_focused.tsv','PValue','logFC')
eff=rd('S2_TF_group_effects.tsv');eff=eff[eff.variant=='DropOverlap']
addrows(eff,'REGULON_GROUP_EFFECT','data/S2_TF_group_effects.tsv','p','coef')
mi=rd('S3_mouse_TF_identity_adjustment.tsv',S);mi['tf']=mi.tf.str.upper()
addrows(mi,'MOUSE_IDENTITY_ADJUSTMENT','source_tables/S3_mouse_TF_identity_adjustment.tsv','p','coef',cohort='Mouse')
assoc=rd('S2_TF_associations.tsv');assoc=assoc[assoc.variant=='DropOverlap']
addrows(assoc,'TF_MODULE_ASSOCIATION','data/S2_TF_associations.tsv','p','rho')
human=rd('prior_sender_TF_focused_with_acceptance.tsv');human['tf']=human.tf_gene
addrows(human,'PATIENT_SENDER_TF','data/prior_sender_TF_focused_with_acceptance.tsv','p','rho',cohort='Human')
comp=rd('Human_CD4_subtype_composition_tests_corrected.tsv',S)
comp['model']='subject_Wilcoxon';comp['tf']='';comp['scope']=comp.subtype
addrows(comp,'CD4_COMPOSITION','source_tables/Human_CD4_subtype_composition_tests_corrected.tsv','p_value','delta_median','fdr_bh',cohort='Human')
accept=pd.DataFrame(records)
accept.to_csv(D/'nominal_acceptance_audit.tsv',sep='\t',index=False)
accept.groupby(['evidence_type','acceptance']).size().rename('rows').reset_index().to_csv(D/'acceptance_inventory.tsv',sep='\t',index=False)
checklist=[
['免疫全景','已有原图；本轮整理收录','人肺初步lineage UMAP/marker图；鼠whole_Broad；髓系原工作注释已追回','existing_figures/human_preliminary_atlas.png; mouse_whole_atlas.png; F08','描述性图谱不要求人为造P；初步大类与已用工作标签须标清；不是疾病新亚型证明'],
['CD4亚群变化','已有样本级组成分析','人CD4 UMAP/受试者组成；鼠CD4 UMAP/组成','existing_figures/human_CD4_composition.png; source_tables/Human_CD4_subtype_composition_tests_corrected.tsv','人CCR7、CD4state、Treg名义P<0.05；五项FDR均未通过；原Source混杂仍须注明'],
['NicheNet','Core主分析及敏感性已完成','Core FDR集合可用；ActCD4 FDR集合不足按原规则跳过；探索集合及受体审计保留','existing_figures/NicheNet_Core_Up.png; NicheNet_Core_Down.png; source_tables/M3_geneset_status.tsv','排名没有P，不作正式显著通信；ActCD4探索结果不能替代失败的主集合'],
['IKZF2/STAT1 TF网络','SCENIC1-4及患者后处理已完成','RNA、活性、预测边、模块、患者sender-TF全部分开整理','figures/F04*.pdf; F05*.pdf; F06*.pdf; F07*.pdf; F09*.pdf','RNA与活性不是同义；各既定模型名义P/FDR分别审计；三seed不作生物重复'],
['焦亡/凋亡/线粒体评分','12冻结模块已完成','MitoISR/MitoROS/Mitophagy/OXPHOS/PyroExec/ApoExec及表型模块','existing_figures/M1_module_scores.png; M2_human_modules.png; figures/F07*.pdf','RNA评分不证明线粒体损伤、蛋白裂解或死亡；本轮不新增评分、阈值或检验']]
pd.DataFrame(checklist,columns=['交付项','状态','已有内容','主要文件','限制与验收']).to_csv(P/'五项交付清单.tsv',sep='\t',index=False)
diagnostics=[
['RNA vs 活性','IKZF2/鼠ActCD4','原RNA All P=0.7781，BG1 P=0.8523；regulon活性在多种原模型名义通过','可估计但该范围RNA差异未获支持；不等同基因无效','不能用Core RNA显著或活性显著替代ActCD4 RNA结论；保留两指标'],
['RNA/队列','STAT1/人Core','原源校正RNA P=0.4071；鼠ActCD4 RNA P=0.0037；人Core活性方向与鼠相反','队列/组织/范围异质性及可估计的不支持结果','不把患者RNA写成上调，不以鼠替代患者验证'],
['RNA/范围','IKZF2/人Core','原RNA P=0.0169、FDR=0.8388；Primary P=0.0061、FDR=0.5083','名义RNA差异支持，未校正显著；与特定亚群活性指标不同','可列名义层级结果，同时保留范围和FDR'],
['患者主线/混杂','IL18–IKZF2 MacAM→CD4-IL7R','All/Matched SourceDep P=0.00726–0.01169；BothEpi P=0.05256–0.07621；df5→3','精度下降与发送深度/上皮混杂敏感并存，不能只归因于功效','完整展示两模型；不声称强调整后仍名义通过'],
['患者样本量','BG1固定主线','n=9，原规则不给相关P；BothEpi df=1、LOO仅少量可拟合','推断条件不足/模型信息量不足','近1 rho不是稳健证据；不通过改筛选追P'],
['多重检验','患者sender–TF全家族','12960行原模型，8208行有FDR；FDR<0.05为0，最小约0.102747','全家族校正未通过','达到名义者只标N；不重算缩小后的两TF家族FDR'],
['RNA检出','两TF','已有每文库/受试者CP10k与检出率；RNA检验表中部分基因/亚群缺失','缺失检验与零表达不同；检测率、样本量和效应分开诊断','缺失记NA；不凭P>0.05或未入表宣判无生物学作用'],
['髓系注释/供者','r4簇2、11及混合簇','原工作标签已使用；2/11约98.8%/97.2%由单个对照供者主导；NoDom已做','供者主导与混合谱系风险有实际证据','保留原标签与NoDom；不命名疾病新亚型，不重新判群'],
['方法实现','SCENIC4依赖与旧标签包','缺matplotlib已修复；旧reviewed_cells.gz确有截断，本轮按同源hash/原映射恢复31120行','有证据的技术错误，已修复','模型表未改；不能把此文件错误当作统计方向相反的解释'],
['分组/来源','现有人肺队列','原标签为CWP/Control，来自G288/G227/G135；已有来源/深度调整','来源与疾病组成非完全平衡','保留原标签，不改写新诊断或制造独立验证队列'],
['NicheNet主集','鼠ActCD4','FDR集合0–2基因<原10基因要求，原流程跳过','既定方法输入条件不足，不是运行错误','保留SKIP和探索性Nominal集合，不降低阈值补显著'],
['外部验证','SCENIC4/GSE183682','每时间点每组1文库；第56天方向一致，第7天不完全一致','仅辅助方向，不能有效组间推断','不作为正式显著验证，不把seed当重复']]
pd.DataFrame(diagnostics,columns=['诊断域','对象','可核验事实','分类','处理原则']).to_csv(P/'未达标原因诊断.tsv',sep='\t',index=False)
caps=json.loads((P/'figure_manifest.json').read_text(encoding='utf-8'))
for info in caps:
 if info['figure']=='F05_regulon_activity_effects' and '不能称显著反向' not in info['caption_zh']:
  info['caption_zh'] += ' 人Core两TF在Source/SourceDep所有展示规格均未达名义P<0.05；STAT1负方向仅为描述，不能称显著反向。'
# concise Chinese handoff note
readme="""IKZF2 / STAT1 技术交付整理 v1
日期：2026-10-02。目标：为巨噬细胞影响肺CD4 T表型、候选基因表达/活性与线粒体及死亡表型的后续基础实验选择整理证据；不宣称因果链已证实。

本轮完成：
1. 五项交付清单、统一队列/亚群/单位/来源及图注。
2. 8张新技术图（PNG+矢量PDF）：原RNA检验、RNA表达/检出、regulon效应、网络展示、鼠/人模块关联、髓系工作标签审阅、固定患者IL18–IKZF2主线敏感性。
3. 原人髓系工作注释与12960行患者sender–TF模型已从Library追回。先前服务器目录未含该聊天工作区后处理包，不代表未完成；本轮不重做模型。
4. reviewed_cells.tsv.gz截断技术修复：原细胞表SHA256完全匹配，恢复31120行，原样本计数逐项一致；其他22文件与模型表字节未变。
5. 按新验收要求建立逐行名义P/FDR审计与未达标原因诊断；不重拟合、不改阈值、不筛样本追P。

验收标记：
NOMINAL_ONLY：已有合法检验P<0.05，但FDR未通过/原规则未计算。FDR_PASS：原完整家族FDR<0.05。NOT_NOMINAL：可估计但P>=0.05。NO_VALID_P：原规则未提供有效P，不能当正式阳性。正式名义层级结果仍是观察性结果，不等于因果证明。三seed/敏感性规格不是独立样本，也不可用其数量计生物学证据强度。

重要结论：
- 小鼠ActCD4：IKZF2 RNA未达名义；STAT1 RNA达名义但FDR未通过。两者regulon活性与RNA必须分开。
- 人Core：IKZF2 RNA达名义，STAT1 RNA未达名义；人STAT1 regulon方向与小鼠不同。
- 固定IL18–IKZF2主线SourceDep达名义，BothEpi未达名义；需同时保留混杂敏感和df下降解释，不能只说功效不足。全患者sender–TF模型无FDR<0.05。
- 人CD4组成：CCR7、CD4state、Treg达名义；全部FDR未通过。图谱描述和负结果完整保留。
- NicheNet Core主集完成，ActCD4主FDR集合不足按原规则跳过；排名/网络不具独立显著性P。
- SCENIC4每时点每组仅1文库，只是辅助方向，不是正式显著验证。

阅读入口：
output/pdf/IKZF2_STAT1_technical_delivery_zh.pdf：中文技术图册。
五项交付清单.tsv、未达标原因诊断.tsv：交付与诊断。
图注与缩写说明.txt：每张图的单位、来源、适用范围和限制。
data/nominal_acceptance_audit.tsv：所有保留模型的逐行名义/FDR标记。
data/myeloid_working_annotation_audit.tsv：原工作标签及marker/供者/父簇核对。
prior_review_repaired/scenic3_review：修复后的原后处理包内容；模型未改。
provenance/：来源文件hash、方法与Library身份记录。
figures/：新技术图；existing_figures/：原图原样收录。
code/：本轮提取、修复、绘图、排版脚本。

工作注释状态：
原reviewed_annotation.tsv已用于已有分析，但用户未逐标签最终确认。本轮只复用和审阅，不重新定义标签、不据此启动任何新模型。MacAM9901、MonoMac5135、MacWide15036；其他16084细胞保留审阅但排除于原巨噬发送范围。簇2/11为供者主导，不是新疾病亚型。

后续需决策：是否将已交付工作标签作为客户版最终图谱标签；固定主线是否作为“名义支持但对更强调整敏感”的实验候选，不能改写为稳健因果轴。若新增方法或改变科学范围，先说明，不静默启动。
"""
(P/'交付说明_中文.txt').write_text(readme,encoding='utf-8')
with (P/'交付说明_中文.txt').open('a',encoding='utf-8') as fh:
 fh.write('\n补充限定：人Core两TF的Source/SourceDep活性效应均未达名义P<0.05；STAT1负方向不能称显著反向。患者sender–TF保留含ETV7比较项的原完整家族FDR，不为IKZF2/STAT1缩小检验家族。审计行数包含配置/seed，不代表独立发现数。\n')
# PDF, landscape A3 for legible scientific panels
pdfmetrics.registerFont(TTFont('CN',r'C:\Windows\Fonts\msyh.ttc',subfontIndex=0))
W,H=landscape(A3);margin=48
out=P/'output/pdf/IKZF2_STAT1_technical_delivery_zh.pdf'
c=canvas.Canvas(str(out),pagesize=(W,H));c.setTitle('IKZF2 / STAT1 五项技术交付与证据审计');c.setAuthor('Project technical synthesis')
page=0
styles={n:ParagraphStyle(n,fontName='CN',fontSize=s,leading=l,wordWrap='CJK',textColor=colors.HexColor('#233442')) for n,s,l in [('body',13,20),('small',10.5,16),('table',11.5,17),('head',18,25)]}
def para(text,x,y,width,style='body'):
 p=Paragraph(text.replace('&','&amp;').replace('\n','<br/>'),styles[style]);w,h=p.wrap(width,H)
 p.drawOn(c,x,y-h);return y-h
def frame(title,tag='技术整理 / 原结果复用'):
 global page
 page+=1;c.setFillColor(colors.HexColor('#173b4d'));c.rect(0,H-92,W,92,fill=1,stroke=0)
 c.setFont('CN',25);c.setFillColor(colors.white);c.drawString(margin,H-48,title)
 c.setFont('CN',11);c.setFillColor(colors.HexColor('#cfe0e7'));c.drawString(margin,H-75,tag)
 c.setFillColor(colors.HexColor('#6c7c84'));c.setFont('CN',9);c.drawString(margin,24,'2026-10-02 | 保留原分组、阈值、seed、排除规则；无新模型拟合');c.drawRightString(W-margin,24,f'{page}')
 return H-118
def imagebox(path,x,y,w,h):
 iw,ih=Image.open(path).size;scale=min(w/iw,h/ih);ww=iw*scale;hh=ih*scale
 c.drawImage(ImageReader(str(path)),x+(w-ww)/2,y+(h-hh)/2,width=ww,height=hh,mask='auto')
def table(rows,widths,y,fontsize=11.5):
 st=ParagraphStyle('cell',parent=styles['table'],fontSize=fontsize,leading=fontsize*1.45)
 cells=[[Paragraph(str(v).replace('&','&amp;').replace('\n','<br/>'),st) for v in row] for row in rows]
 t=Table(cells,colWidths=widths,repeatRows=1,hAlign='LEFT')
 t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e8f0f3')),('GRID',(0,0),(-1,-1),.5,colors.HexColor('#d4dfe4')),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),8),('RIGHTPADDING',(0,0),(-1,-1),8),('TOPPADDING',(0,0),(-1,-1),8),('BOTTOMPADDING',(0,0),(-1,-1),8)]))
 _,h=t.wrap(W-2*margin,H);assert y-h>40,('Table overflows',page,y,h);t.drawOn(c,margin,y-h);return y-h
def newpage():c.showPage()
y=frame('IKZF2 / STAT1：五项交付与证据等级','技术整理 v1 | 配套图册、完整审计表及旧标签文件技术修复')
y=para('目标：巨噬细胞信号与肺CD4 T细胞候选基因表达/活性、线粒体及死亡表型的证据整理，支持后续基础实验选择。观察关联尚不能确立整条因果顺序。',margin,y,W-2*margin)-16
rows=[['交付项','本轮整理状态','必须保留的限定']]+[[r[0],r[1]+'；'+r[2],r[4]] for r in checklist]
y=table(rows,[135,450,W-2*margin-585],y,12)-16
y=para('验收：已有有效名义P<0.05才标名义层级阳性；FDR另列。图谱、网络排名、无P及负结果照实展示，不以趋势包装正式阳性。不按最小P选模型，不将seed算作独立生物学重复。',margin,y,W-2*margin)-10
para('关键修正：人髓系工作注释和患者sender–TF模型已在原聊天工作区交付，本轮已追回，未重做。旧标签gzip截断已同源修复为31,120行；原模型表及其他22文件未变。',margin,y,W-2*margin);newpage()
y=frame('01-02 | 全景与CD4亚群：原图原样收录')
items=[('human_preliminary_atlas.png','人肺：原初步lineage全景，含非免疫肺细胞'),('mouse_whole_atlas.png','小鼠：GSE314943原Broad全景'),('human_CD4_composition.png','人肺：受试者CD4组成（原校正身份）'),('mouse_CD4_composition.png','小鼠：原Core CD4范围的组成')]
bw=(W-2*margin-28)/2;bh=245
for idx,(name,txt) in enumerate(items):
 col=idx%2;row=idx//2;x=margin+col*(bw+28);bottom=H-135-(row+1)*(bh+30)
 imagebox(P/'existing_figures'/name,x,bottom,bw,bh)
 para(txt,x,bottom+bh+21,bw,'small')
para('人组成检验单位是受试者，原Wilcoxon+BH：CD4-CCR7 P=0.0455、CD4state P=0.0260、Treg P=0.0370，三项名义通过但FDR均约0.0758；CD4-GZMA、CD4-IL7R未达名义。来源混杂仍须保留。鼠组成以保守Core为分母，原六文库置换检验最小非零P=0.1，仅描述，不外推时间轨迹。',margin,96,W-2*margin,'small');newpage()
y=frame('03 | NicheNet：Core主分析与ActCD4探索分开')
imagebox(P/'existing_figures/NicheNet_Core_Up.png',margin,335,bw,340);imagebox(P/'existing_figures/NicheNet_Core_Down.png',margin+bw+28,335,bw,340)
rows=[['接收范围','主FDR集合','状态与解读'],['Core / All、BG1','Up/Down集合均达到原条件','已有主分析排名；AUPR是排序指标，不是通信P或因果概率'],['ActCD4 / All、BG1','Up为1–2基因，Down为0；低于原10基因要求','原流程SKIP；名义集合的探索排名保留，但不替代失败主集合'],['受体/发送支持','原表达支持LR、复合受体共检出及5%/10%敏感性均已做','转录可行性不等于蛋白复合物、空间接触或受体→TF机制']]
table(rows,[180,330,W-2*margin-510],306,12)
para('来源：mechanism_step3（Core图与geneset_status）、Step4受体审计、Step5候选移除敏感性。原模型未重跑；旧候选基因表不机械改名为IKZF2/STAT1通路。',margin,82,W-2*margin,'small');newpage()
y=frame('05 | 焦亡、凋亡、线粒体与T细胞表型：冻结评分')
imagebox(P/'existing_figures/M1_module_scores.png',margin,265,bw,420);imagebox(P/'existing_figures/M2_human_modules.png',margin+bw+28,265,bw,420)
para('原12模块完整保留：OXPHOS、MitoISR、MitoROS、Mitophagy、PyroExec、ApoExec、TAct、TReg、Inhib、Inflam、Naive、Cycle。模块为冻结转录基因面板及匹配表达对照评分，不是新通路激活测定。原独立文库/受试者单位、覆盖率条件、至少20细胞资格与All/BG1规则均未更改。',margin,220,W-2*margin)
para('MitoROS表示抗氧化应答转录信号；PyroExec/ApoExec不证明蛋白裂解或真实死亡。不同队列评分尺度不合并。后续TF–模块图展示既有相关，不能由相关推断“TF改变导致线粒体损伤/焦亡”。',margin,143,W-2*margin)
para('来源：source_tables/M1_module_definition.tsv、M1_sample_scores.tsv、M1_group_effects.tsv、M2_human_source_models.tsv及原方法说明。',margin,70,W-2*margin,'small');newpage()
order=['F04_RNA_existing_DE','F04b_RNA_level_detection','F05_regulon_activity_effects','F06_regulon_network_display','F07_mouse_module_associations','F07_human_module_associations','F09_patient_fixed_axis_robustness','F08_myeloid_working_label_audit']
for name in order:
 info=next(x for x in caps if x['figure']==name)
 frame(info['title'],name+' | 原结果可视化，不重拟合')
 imagebox(F/(name+'.png'),margin,160,W-2*margin,H-285)
 para(info['caption_zh'],margin,143,W-2*margin,'small')
 para('来源：'+info['source'],margin,66,W-2*margin,'small')
 newpage()
y=frame('验收与诊断：哪些达名义，哪些不能下正式阳性结论')
rows=[['证据','已有结果','诊断/处理']]
rows += [
 ['RNA：鼠ActCD4','IKZF2 P=0.778/0.852；STAT1 P=0.0037/0.0059（All/BG1）；FDR均未通过','IKZF2该范围RNA差异未获支持；不能以Core或regulon替代。STAT1可标名义RNA差异。'],
 ['RNA：人Core','IKZF2 P=0.0169、FDR=0.8388；STAT1 P=0.4071','前者为名义层级，后者未名义；不据不显著宣判基因无效。人鼠/亚群/指标不同。'],
 ['固定患者IL18–IKZF2','SourceDep P=0.00726–0.01169；BothEpi P=0.05256–0.07621；df5→3','基础调整达名义，更强调整未达；精度与混杂敏感同时保留，不能说仅功效不足。'],
 ['BG1与全家族','固定主线BG1 n=9，无相关P；全12960行模型无FDR<0.05','NO_VALID_P不能当正式阳性；保留原完整家族FDR，不缩家族追显著。'],
 ['技术问题','SCENIC4依赖已修复；旧标签gzip按同源hash恢复31120行，模型未变','有证据的方法/文件错误已处理；不能用该错误解释所有不一致生物学结果。'],
 ['外部与实验选择','SCENIC4每时间点每组1文库，仅方向辅助；整条因果链未证明','名义支持、负结果及敏感性一起用于实验选择；不新增因果/中介/轨迹方法。']]
rows[2]=['RNA/活性：人Core','RNA：IKZF2 P=0.0169，FDR=0.8388；STAT1 P=0.4071；两TF的Source/SourceDep活性均未达名义','IKZF2 RNA为名义层级；STAT1活性负方向不能称显著反向。不据不显著宣判基因无效。']
y=table(rows,[175,450,W-2*margin-625],y,12)-16
para('当前需要定稿的仅是科学表述与工作标签状态：是否将原已用工作注释作为客户版最终图谱标签；是否把固定主线按“名义支持、但对更强调整敏感”作为实验候选。技术整理及文件修复无需重跑原模型。',margin,y,W-2*margin)
para('完整审计：data/nominal_acceptance_audit.tsv；未达标原因诊断.tsv；五项交付清单.tsv。修复审计位于prior_review_repaired/scenic3_review/technical_repair_audit.json。',margin,75,W-2*margin,'small')
c.save()
summary={'pdf':str(out),'pages':page,'new_figures':len(caps),'acceptance_rows':len(accept),'new_models_fit':0,'source_hash_count':len(rd('source_manifest.tsv',P/'provenance')),'reviewed_cells':31120,'reviewed_sender_model_rows':12960}
(P/'technical_delivery_status.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
print(accept.groupby(['evidence_type','acceptance']).size().to_string())
print('HUMAN_ACTIVITY_CORE_SOURCE')
print(eff[(eff.species=='Human')&(eff.scope=='Core')&(eff.model=='Source')][['filter','seed','tf','coef','p','FDR']].to_string(index=False))
