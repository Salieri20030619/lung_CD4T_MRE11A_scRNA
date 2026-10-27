from pathlib import Path
import json,html,hashlib,collections
import pandas as pd
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,PageBreak,Table,TableStyle,KeepTogether
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors
import fitz
W=Path('C:/Users/31937/lung_review_20261017_v14');O=W/'revision_v14';T=W/'render_response';T.mkdir(exist_ok=True)
items=[
{'id':'R01','severity':'高','type':'必须修复（解释）','status':'已解决解释；功能连接仍待验证','location':'v13第24–28页；v14第27–32页，第5–5.3节','evidence':'AMlike CXCL16差异P=0.002370、BH=0.1065，而对应ActCD4 STAT1关联P=1；AllMac关联P=0.002778、BH=0.4811，但差异P=0.149。ActCD4上调目标的CXCL16来源内NicheNet排名为199/280。','impact':'不同发送端的证据不能串成同一已建立机制，NicheNet不能支持首选驱动措辞。','action':'将“接收端调控明确后”改为“确定接收端候选后”；补入来源断点与预测排名，保留CXCL16–CXCR6、STAT1–GSDMD两条独立验证方向。','files':'S11；S19_NicheNet_all_source_rankings.tsv；S23–S24','acceptance':'核对来源、接收范围和排名；正文无“多方法一致首选驱动”或完整因果链成立的表述。'},
{'id':'R02','severity':'高','type':'必须修复（核心证据边界）','status':'已补核心结果；信号来源未解决','location':'v13第23、30、35页；v14第24–25、35、40页','evidence':'从原始计数复核Saline_2：424个Core细胞全部检出SFTPC与SCGB1A1；GSDMD阳性31个且全部1 UMI。严格排除后该样本剩0细胞。患者3个双阳性效应CD4伴髓系混合标记。','impact':'影响GSDMD是否来自纯CD4内源性表达，分母或背景基因敏感性无法替代来源判断。','action':'风险前移至4.4节，新增原始UMI核对、身份和严格过滤名单，明确六样本敏感性无法评价；未将严格过滤写成通过，未擅自做无空滴依据的背景扣除。','files':'S22_GSDMD_UMI_source_check.tsv；S22_strict_identity_sensitivity_roster.tsv；S16','acceptance':'确认31/424及每细胞1 UMI、0个保留细胞可回查；结论不确认CD4内源性GSDMD。仍需空滴/原始捕获资料或纯化CD4蛋白/原位证据。'},
{'id':'R03','severity':'高','type':'必须修复（统计解释）','status':'已解决说明','location':'v13第11页；v14第11、13、40页','evidence':'3对3枚举20种标签分配；双侧绝对均值差使完全对调分配成对，Pmin=2/20=0.10。已用独立枚举核对。','impact':'该检验设计不能达到P<0.05，未显著不等于无方向或无变化。','action':'方法准确写明分辨率；结果与讨论强调效应、方向与样本分布。保持原检验和阈值，不改单侧、不为显著性更换方法。明确不适用于edgeR或720种相关检验。','files':'provenance/v14/local_numeric_validation.json；S04–S05；S21','acceptance':'数学下限、适用范围和正文解释一致；原差异基因证据未被该下限否定。'},
{'id':'R04','severity':'高','type':'必须修复（模型与证据口径）','status':'已解决；原始组成调整模型新增澄清','location':'v13第23页；v14第24–26页','evidence':'四模型均n=6：原始组成合并P=0.016667/BH=0.066667；原始组成调整组别P=0.277778/BH=0.555556；固定组成合并P=0.002778/BH=0.027778；固定组成再调整组别P=0.055556/BH=0.555556。四项从样本输入复算通过。','impact':'评审所列0.05556必须限定为固定组成后再调整组别；不能泛化为所有组别调整结果或稳定校正后显著。','action':'新增四模型并列表，保留特定补充家族的BH通过与共同组别解释；明确种子和重叠亚群不构成独立重复；IKZF2仅作组成/状态线索。','files':'S21_key_exact_P_recalculation.tsv；S21_complete_STAT1_gene.tsv；S21_complete_STAT1_fixed.tsv；S22模型输入','acceptance':'模型、P、BH、样本及家族一一对应；固定组成和疾病组别调整不再混用。'},
{'id':'R05','severity':'中至高','type':'定向补充（机制连接）','status':'部分解决：可计算证据已补，因果/蛋白连接仍缺','location':'v14第12、33页，第1.9与5.4节','evidence':'冻结15项配体—近端TF组合×2来源×4范围×2模型=240项，200项可评估，40项NA；7项名义显著，0项BH通过。IL15发送/接收组分有检出，但seed17无STAT5A/B regulon；IL18–NFKB1为负向名义关联；AllMac TNF–ActCD4 NFKB2为正向名义关联。','impact':'补充支持通路分层及缺口定位，不能证明巨噬信号驱动STAT1/IKZF2或焦亡。','action':'按IL15、TNF、IL18、IL27、IFNB1、IFNG逐轴呈现完整受体、发送端复合物、近端活动与所有阴性/NA。主要种子未评估的STAT5不临时换种子挑结果。','files':'S23全部计划配对与样本输入；S24原始共检出；S22发送端表达','acceptance':'核对完整复合物与样本分母；受体缺检不等同绝对不存在。蛋白跨呈递、磷酸化、IL18→IFNG时序及到IKZF2/STAT1功能联系仍需实验。'},
{'id':'R06','severity':'中','type':'必须澄清（创新性与来源）','status':'已收窄创新性并补来源；来源特异作用未验证','location':'v14第34、37页，第5.5与7节','evidence':'2024年矽肺CD4驻留研究已涉及CXCR6；2019年肺CD8研究涉及上皮CXCL16定位。全肺冻结标签病例等权CP10k：CXCL16 Mac=2.5848、DC=3.0685、AT1=2.4172；IFNG NKLike=3.7019、T=0.4536、Mac=0.0038。','impact':'CXCR6肺部驻留不是本项目新发现；巨噬表达不能证明特异发送源；T/NK IFNG需保留为竞争解释。','action':'复核原始文献链接，新增所有Broad标签逐样本描述；创新问题限定为阶段、巨噬亚群、CD4状态及来源可区分的功能贡献。','files':'S26_whole_lung_competing_sources.tsv；参考文献[2,11]','acceptance':'报告不把CD8文献直接当CD4验证，不声称巨噬特异性；新来源分析使用冻结Broad标签，未重新验证注释且不作组织分泌量推断。'},
{'id':'R07','severity':'中','type':'补充材料可复核性','status':'约定复核范围内已完成','location':'v14第12、41–42页；S21–S27及本说明','evidence':'106个保存的完整检验家族BH复算均通过；Core STAT1–GSDMD四项精确P复算通过。182个基准文件内容哈希与服务器一致。补充表包含六样本名单、完整检验池、基因集及匹配对照。','impact':'评审不再只凭PDF表述判断统计实现；核心数值和家族可以实际复核。','action':'保留S01–S20，新增S21–S27完整表，不仅提供显著子集；新增样本输入、NA原因、重叠成员及哈希。代码与完整证据在服务器修订目录另存。','files':'S21_BH_recalculation_audit.tsv；S21_input_manifest.tsv；S22_sample_scope_roster.tsv；S25–S27','acceptance':'复算家族必须使用完整池；本轮未重跑QC、注释、edgeR拟合、SCENIC网络及NicheNet训练，不把BH复算通过称为全流程独立审计通过。'},
{'id':'R08','severity':'中','type':'有依据的新增问题（评分构造）','status':'已补敏感性；探索性定位保留','location':'v14第12、26页，第1.9与4.5节','evidence':'DropOverlap模块成员交集为0，但seed17 STAT1与MitoROS匹配对照重叠13个、IKZF2与ApoExec对照重叠15个。重现336个原评分后去除对应regulon/TF对照，STAT1–MitoROS方向保留；IKZF2–ApoExec固定组成rho=0.0286、P=1。','impact':'模块成员去重叠不自动消除匹配对照重叠；IKZF2直接凋亡解释进一步减弱。','action':'新增全部三种子重叠成员、GO独立说明、固定对照剔除敏感性；新BH仅在七模块敏感性家族内，原结果不得用新小家族重标BH。未改写原统计与图件。','files':'S25_regulon_module_control_overlap.tsv；S25_GO_overlap_not_in_TF_test_family.tsv；S27全表','acceptance':'核对原评分重现误差、成员及新家族；GSDMD单基因终点不使用该匹配对照评分，不能把此敏感性当作解决GSDMD来源问题。'}]
(O/'review_response.json').write_text(json.dumps(items,ensure_ascii=False,indent=2),encoding='utf-8');pd.DataFrame(items).to_csv(O/'review_response.tsv',sep='\t',index=False,encoding='utf-8')
pdfmetrics.registerFont(TTFont('CJK','C:/Windows/Fonts/ARIALUNI.ttf'))
style=ParagraphStyle('body',fontName='CJK',fontSize=10.5,leading=17,spaceAfter=9,wordWrap='CJK');h=ParagraphStyle('head',parent=style,fontSize=18,leading=25,textColor=colors.HexColor('#163b50'),spaceAfter=15);sub=ParagraphStyle('sub',parent=style,fontSize=12,leading=19,textColor=colors.HexColor('#28678a'))
story=[]
def p(t,sty=style):story.append(Paragraph(html.escape(t),sty))
p('独立评审逐项回应与补充分析说明',h);p('v14｜2026年10月17日｜供再次独立评审')
p('基准：silicosis_integrated_report_v13.pdf（37页）。修订：silicosis_integrated_report_v14.pdf（43页）。本次依据反馈1017.txt，保留原18幅图和S01–S20的数值，在独立目录补充分析与解释。')
p('结论仍为候选筛选、证据整合和实验优先级：未建立完整巨噬细胞—CD4损伤/焦亡机制。统计口径、风险呈现和可复核材料已补足；细胞来源、蛋白活化及因果连接仍需验证。')
p('最重要的新增结果',sub)
for s in ['106个完整BH家族与4种关键精确P复算通过。','240项计划近端配对中200项可评估，7项名义显著，0项BH通过；全部阴性及40项NA保留。','匹配对照存在regulon重叠；去重叠敏感性保留STAT1–MitoROS方向，但IKZF2–ApoExec固定组成后接近零。','全肺CXCL16有巨噬、DC和上皮等多来源；T/NK IFNG是竞争解释。','0.05556对应固定组成后再调整疾病组别；原始组成直接调整组别P为0.27778。']:p('• '+s)
p('资料使用',sub);p('主报告正文已整合关键新增结果；原S01–S20在supplement目录，新增S21–S27在supplement_v14目录。附表为UTF-8无BOM，使用制表符解析并保留NA。完整脚本、来源与日志保留于服务器review_revision_20261017_v14目录。')
p('本说明的“已解决”只针对列明的报告/计算问题，不代表独立评审已通过。下页起逐项给出严重程度、定位、证据、影响、修订和复核要求。')
for z in items:
 story.append(PageBreak());p(z['id']+'｜'+z['type'],h);p('严重程度：'+z['severity']+'　　状态：'+z['status'],sub)
 for label,k in [('文件与位置','location'),('证据或依据','evidence'),('对结果的影响','impact'),('已完成修订/补充','action'),('对应材料','files'),('复核验收与未决事项','acceptance')]:p(label,sub);p(z[k])
def footer(c,d):
 c.setFont('CJK',8);c.setFillColor(colors.HexColor('#42596a'));c.drawString(40,22,'v14独立评审回应｜原成果保留｜待再次评审');c.drawRightString(802,22,str(d.page))
f=O/'review_response_v14.pdf';SimpleDocTemplate(str(f),pagesize=(842,595),leftMargin=40,rightMargin=40,topMargin=35,bottomMargin=40,title='v14独立评审逐项回应').build(story,onFirstPage=footer,onLaterPages=footer)
d=fitz.open(f)
for i,page in enumerate(d):
 assert '\ufffd' not in page.get_text();page.get_pixmap().save(T/f'response_{i+1:02d}.png')
print('RESPONSE_PAGES',len(d))
# Correct duplicate-content path matching: exact paths take precedence.
r=json.loads((W/'analysis/provenance/baseline_remote_sha256.json').read_text(encoding='utf-8'));local={f.relative_to(W/'baseline_v13').as_posix():hashlib.sha256(f.read_bytes()).hexdigest() for f in (W/'baseline_v13').rglob('*') if f.is_file()};rows=[]
for rp,sha in r.items():
 lp=rp if local.get(rp)==sha else next(p for p,h in local.items() if h==sha);rows.append(dict(remote_path=rp,local_path=lp,sha256=sha,name_changed=rp!=lp))
pd.DataFrame(rows).to_csv(O/'provenance/v14/baseline_filename_map.tsv',sep='\t',index=False,encoding='utf-8');v=O/'provenance/v14/local_numeric_validation.json';j=json.loads(v.read_text());j['baseline_filename_mapping']=sum(z['name_changed'] for z in rows);v.write_text(json.dumps(j,indent=2),encoding='utf-8')