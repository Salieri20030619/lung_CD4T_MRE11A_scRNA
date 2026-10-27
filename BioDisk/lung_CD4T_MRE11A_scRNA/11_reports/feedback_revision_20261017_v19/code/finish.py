from pathlib import Path
import json,shutil,hashlib,zipfile
import pandas as pd
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,PageBreak,Table,TableStyle
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from xml.sax.saxutils import escape
O=Path(__file__).parent;H=O.parent;V=H/'lung_append_analysis_20261017_v18';B=H/'lung_three_cohort_20261017_v16/data';S=H/'lung_dual_benchmark_20261017_v17'
D=H/'Downloads/线粒体焦亡_反馈定向修订_v19_20261017';D.mkdir(exist_ok=True);K=O/'review_bundle';K.mkdir(exist_ok=True)
def rd(p):return pd.read_csv(p,sep='\t',low_memory=False)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def cp(p,rel):
 t=K/rel;t.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,t)
a=rd(O/'tables/R05_current_all90.tsv');f=rd(O/'tables/R05_current_focus_conclusions.tsv');q=json.loads((O/'reconstruction_QA.json').read_text(encoding='utf-8'))
assert len(a)==90 and a.gene.nunique()==90
assert [a.A_status_current.str.startswith(x).sum() for x in ['A1','A2','A3']]==[52,37,1]
assert len(f)==10 and f.gene.nunique()==9 and q['max_error']<1e-12
assert len(rd(O/'tables/R02_all18_routes_evidence.tsv'))==18
assert not rd(O/'tables/R02_corrected_PROS1_alternatives.tsv').group_expression_pass.any()
sections=[('R01','交集定义与筛选追溯','原第1-3、6-9页','R01_selection_flow.tsv；R01_candidate_DE_all_tested_records.tsv；R01_focus_selection_record.tsv','完整受检9091基因/52448条主要DE记录 → 原始P<0.05的2765基因/6126记录 → 功能111/274 → 原v16网络交集55/155。55与历史42取并集，重叠7，得到90；此步为扩展，不是嵌套淘汰。90中89具名义DE（344记录），MRE11单列A3；功能筛后57/172，32入A2；同状态入口后56/168，IL10另入A2；来源复核后54/164，GSDMD/S100A8入A2；自身损伤与输出作用分开后52/160，GZMB/TNFSF10入A2。A2=32+1+2+2=37。记录单位为基因-状态-对比，90并集为基因行，重点组合为基因-状态。52个A1是计算支持入池，不是52个跨队列一致差异基因。重点9基因源于此前5个探索性优先/备选加4个历史补评；FOXP1有2个状态，共10组合，不全来自A1。v18冻结于2026-10-17 04:46:56 UTC，即北京时间12:46:56，属于发现后冻结，不是项目起始预设。完整DE保留原家族FDR、来源行号和模型。功能集合和NicheNet不是额外独立队列；交集不产生新的联合P或集合FDR。'),
('R02','撤回PROS1-AXL的本项目功能支持','原第17-18页','R02_all18_routes_evidence.tsv；TAM_receiver_sample_coverage.tsv；R02_upstream_statistics_annotated.tsv','实际PROS1资格使用了AXL。经典小鼠实验显示PROS1激活TYRO3/MERTK，在所测条件下不激活AXL；GAS6可激活AXL（Lew等2014，PMID25265470，https://pmc.ncbi.nlm.nih.gov/articles/PMC4206827/）。数据库PROS1_AXL边引用PMID29531161，原研究在间质型胶质瘤培养情境观察PROS1与AXL作用，因此不能说该边毫无文献依据；问题是不能外推为本项目小鼠CD4功能入口。已撤回该功能支持。ActCD4六样本细胞数36、121、35、284、100、225；TYRO3阳性数0、4、0、0、0、2；MERTK均0。要求合格样本至少10细胞、组内至少一半样本达到5%阳性，Raw/BG1均未通过。零转录检出不等于受体蛋白不存在。历史AMlike PROS1-GSN rho=0.943、P=0.0167、FDR=0.150，疾病调整P=0.278，保留为历史协变，不再充当功能入口。撤回PROS1后原55基因和52个A1仍有其他既有模型入口，计数不变，GSN继续备选。18路线分列受体必要组分、物种、来源和引用；部分条目只有KEGG来源，如实保留原始功能文献缺口。数据库边、受体表达、功能激活不等价，不把其余17条视为已验证功能路线。'),
('R03','评分顺序与权重可唯一重建','原第11、14、19页','R03_gene_counts_with_library.tsv.gz；R03_score_reconstruction_QA.tsv；R03_equal_state_weights_actual.tsv','令x_gc为细胞原始计数：C_gs=Σ_c x_gc；L_s=Σ_c L_c；E_gs=10000×C_gs/L_s。分母是完整原始文库UMI，包含未映射/歧义基因，不是仅模块成员计数和。候选表达G_gs=ln(1+E_gs)；模块M_ms=mean_g ln(1+E_gs)，均值维度为保留且可映射模块基因。先按样本汇总计数再归一化，不是逐细胞log后平均。零表达参与均值，缺失映射不填零；至少3个映射成员且覆盖≥50%才能评价。9候选从模块移除，实际重叠者DDIT4位于MitoISR/Glycolysis；所有成员、覆盖与共同成员家族随附。Core等权状态标准化：Naive/ActCD4/Treg各1/3；候选先平均三状态E再log，模块先计算每状态M再平均，AUCell平均每状态样本AUC。v14旧权重来自六样本平均比例，分别0.466639、0.309979、0.223382，与1/3不等价。v14逐细胞匹配对照评分后聚合，与v18样本汇总程序表达尺度不同，不能直接比较数值。2852评分值重建最大绝对误差8.88×10^-16；实际权重和细胞数逐样本随附。'),
('R04','联合模型及置换解释补全','原第11-15、19页','R04_actual_sample_rosters.tsv；R04_design_matrices.tsv；R04_permutation_blocks.tsv；code/v18/analyze.py','联合关联n=10，排列6!×2!×2!=2880；联合组间比较C(6,3)×2×2=80；疾病调整关联3!×3!=36；时间组内排列(3!)^4=1296，去一个样本为2!×(3!)^3=432。183682每时间点仅1+1，只报告方向和效应；疾病调整时四个单例单元不贡献组内交换。联合n=10不能称为独立复现。34个设计的实际名单、顺序、n及排列数与既有结果一致，截距、队列/时间层、疾病项、残差化矩阵与置换块均导出；组间比较标签分配也随附。MDFIC调整后rho=0.950、P=0.0833，表述为效应较强、精确检验证据有限；不能断言未显著完全由样本量或离散性造成。近似t检验保持单列，不替代预设置换P。主分析、疾病调整、组成、时间及共同成员各保留自己的FDR家族。分析端设计核对不等于独立审阅已验证可交换性或全部P值。'),
('R05','按具体证据更新候选定位','原第6、10、13-16页','R05_current_focus_conclusions.tsv；R05_current_all90.tsv','当前总表并列疾病DE、状态联系、损伤程序、上游入口、第二队列各时间方向、时间队列及未决问题；保留原完整家族FDR。旧统一优先名单标为历史信息，不作为当前排名。后文列出重点结论，完整数值和来源索引见10行组合表。'),
('R06','A/B命名修订','原第1-2、5-9、20页','R05_current_all90.tsv','A表示候选资格，B改称完整机制支持程度。本轮尚无候选达到完整机制支持标准，不等于没有选题价值。P<0.05可作为名义信号，不新设所有队列、所有敏感性均显著或统一FDR、网络排名淘汰门槛。受体适用性、来源身份等实质问题仍影响分层。当前定位不是已经完成创新性、可操作性和实验投入比较的最终排名。'),
('R07','提供可复核材料','复核材料此前未提供','README_复核入口.txt；MANIFEST_SHA256.tsv；QA_v19.json','复核包包含冻结协议、候选选择/流程、完整DE家族与设计、样本评分、模块成员与覆盖、sender表达、18路线、受体资格、clean AUC与成员、完整关联家族、共同成员敏感性、权重、名单、数值QA和关键脚本。无嵌套旧压缩包，不重复原始细胞矩阵/大型网络对象。可复核样本级统计与评分，不能仅凭本包从细胞层重跑所有SCENIC。原审阅未独立复算P/FDR；176个BH家族/4项P核查仍明确为v18分析端记录，不冒充独立审核通过。')]
rr=pd.DataFrame([dict(ID=i,severity='必须修复/补充',status='已修订，待独立复核',original_location=loc,revision_location=title,evidence_files=files,basis_impact_current_conclusion=body,acceptance='按对应附件追溯并复核；未改变的原统计保留') for i,title,loc,files,body in sections]);rr.to_csv(O/'tables/R01_R07_response.tsv',sep='\t',index=False,encoding='utf-8')
for p in (O/'tables').iterdir():cp(p,p.name)
for p in (V/'tables').iterdir():
 if p.is_file():cp(p,('historical_v18/' if p.name=='frozen_route_LR_coverage.tsv' else '')+p.name)
for p in (V/'review_tables').glob('*.tsv'):cp(p,'historical_v18/'+p.name)
for n in ['protocol_frozen.json','coverage_sensitivity_record.json','QA.json','pdf_QA.json']:cp(V/n,n)
for n in ['reconstruction_QA.json','scope.json']:cp(O/n,'v19_'+n)
cp(O/'inputs/reviewer_feedback.txt','reviewer_feedback.txt')
for p in (O/'inputs').glob('PMID_*.json'):cp(p,'literature/'+p.name)
for p in (V/'inputs').iterdir():
 if p.is_file() and p.suffix in ['.json','.tsv']:cp(p,'inputs/'+p.name)
for n in ['resource.R','export_scores.py','aucell_clean.py','analyze.py','source_checks.py','tf_effect_and_time.py','common_gene_sensitivity.py','verify_numbers.py','normalization_audit.py']:cp(V/n,'code/v18/'+n)
for p in O.glob('*.py'):cp(p,'code/v19/'+p.name)
for p in (B/'code').glob('*'):
 if p.is_file():cp(p,'code/v16/'+p.name)
for n in ['assess.py','network.R','expression.py']:cp(S/n,'code/v17/'+n)
cp(B/'protocol_frozen.json','protocols/v16_protocol_frozen.json');cp(S/'protocol_v17.json','protocols/v17_protocol.json');cp(B/'tables/DE_all_models.tsv','original_statistics/DE_all_models.tsv')
for p in (B/'tables').glob('*__design.tsv'):cp(p,'model_designs/'+p.name)
for n in ['candidate_complete_routes.tsv','function_membership_extended.tsv','function_literature_extensions.tsv','three_axis_candidates_all_contrasts.tsv','nominal_DE_axis_membership.tsv','sample_inventory.tsv','model_diagnostics.tsv']:cp(B/'tables'/n,'v16_trace/'+n)
for p in (S/'tables').iterdir():
 if p.is_file() and p.suffix in ['.tsv','.gz','.json']:cp(p,'v17_trace/'+p.name)
p=H/'lung_review_20261017_v14/revision_v14/supplement/S18_fixed_composition_weights.tsv';cp(p,'historical_v14/'+p.name)
for p in (V/'logs').glob('*environment*.json'):cp(p,'environment/'+p.name)
qa={'status':'analysis-side targeted QA passed; independent review pending','A_counts':[52,37,1],'focus':[9,10],'routes':18,**q,'statistics_unchanged':{n:sha(V/'tables'/n)==sha(K/n) for n in ['all_associations.tsv','programme_group_effects.tsv','common_gene_sensitivity.tsv','sample_gene_and_program_scores.tsv']}}
assert all(qa['statistics_unchanged'].values())
(O/'QA_v19.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8');cp(O/'QA_v19.json','QA_v19.json')
(K/'README_复核入口.txt').write_text('v19定向修订，待独立复核；非客户版、非最终投入排名。\n先读R01_R07_response.tsv及R05两张当前结论表。\n原统计在all_associations等表；完整DE在original_statistics/，对应设计在model_designs/。\nR03计数分母+模块成员可重建评分，R04保存实际设计和置换块。历史PROS1结果在标注表及historical_v18/，不得作为当前功能支持。\ncode/保留原项目绝对路径，异地复算需映射本包与原项目路径。原始细胞矩阵和大型对象不在本包。\nQA.json及numerical_QA.tsv是v18分析端记录，QA_v19.json为本次记录，均不是独立审阅通过证明。\nMANIFEST_SHA256.tsv核查所有附件完整性。\n',encoding='utf-8')
pdfmetrics.registerFont(TTFont('CN',r'C:\Windows\Fonts\ARIALUNI.TTF'))
sty=ParagraphStyle('p',fontName='CN',fontSize=10.5,leading=18,spaceAfter=11,wordWrap='CJK');title=ParagraphStyle('h',parent=sty,fontSize=18,leading=27,spaceAfter=18,textColor=colors.HexColor('#153c50'));small=ParagraphStyle('s',parent=sty,fontSize=8.5,leading=14)
story=[]
def p(t,s=sty):story.append(Paragraph(escape(t).replace('\n','<br/>'),s))
def page():story.append(PageBreak())
p('矽肺线粒体损伤/焦亡项目\nv19 定向修订与当前结论',title)
p('2026-10-17｜内部独立复核材料｜已修订，待独立复核')
p('本轮保留52个计算支持入池基因，A2=37、A3=1。重点9基因、10组合形成凋亡转录程序、Treg状态及代谢等具体研究问题，不是统一最终选题排名。')
p('实质修正：撤回PROS1-AXL在本项目小鼠CD4中的功能入口支持，保留历史相关统计。经典受体TYRO3/MERTK未达到既定转录覆盖门槛。GSN仍有其他既有模型入口，因此继续备选；不得将替代模型边视为功能已验证。')
p('此次仅处理R01-R07及必要复核，不重跑全套QC、注释或SCENIC。2852个评分值重建一致，最大误差8.88×10^-16；34关联设计与保存的样本顺序、n和排列数一致。原P/FDR保留。')
p('v16/v17/v18原件保留。本报告更新受影响解释；历史PDF中的旧优先名单、PROS1功能入口和旧B类命名不再代表当前结论。完整TSV和代码见配套复核包。')
p('当前没有候选达到完整机制支持标准，不等于没有选题价值。分析端QA不等于独立评审通过；本版等待按原问题编号复核。')
for i,t,loc,files,body in sections:
 page();p(i+'｜'+t,title);p('状态：已修订，待独立复核。原位置：'+loc,small)
 # Break long CJK body into meaningful sentence groups without losing text.
 sentences=body.split('。');chunk=''
 for s in sentences:
  if not s:continue
  chunk+=s+'。'
  if len(chunk)>180:p(chunk);chunk=''
 if chunk:p(chunk)
 p('修订证据与验收入口：'+files,small)
page();p('重点候选当前结论（1/2）',title)
items=[('TTC39C / ActCD4','凋亡相关转录状态候选。ApoExec联合rho=0.943，P=0.00347，FDR=0.0556；疾病调整rho=0.922，P=0.0278，FDR=0.556，属于名义信号。D7方向随归一化口径变化，D56反向；上游未补齐。TLike的MitoISR时间联系不算同一ActCD4凋亡关联复现。'),('MDFIC / ActCD4','损伤转录关联备选。ApoExec联合rho=0.947，P=0.00278，FDR=0.0556；调整rho=0.950，P=0.0833，效应较强、精确检验证据有限。第二队列D7/D56疾病方向相反；时间队列ApoExec日龄调整P=0.444，支持不足。'),('SKIL / Treg','Treg活化/状态调控。TAct联合rho=0.965，P=0.00208，FDR=0.0744；调整rho=0.897，P=0.0278，FDR=0.868。p53/内源性凋亡段未补齐。AMlike CCL2-SKIL P=0.0167，FDR=0.150；调整P=0.278。不能与Broad Mac来源优势跨范围拼接。'),('IL1R2 / Treg','代谢与反馈调节。ENO1 RNA联合rho=0.912，P=0.000694，FDR=0.0285；疾病调整P=0.111。PyroExec P=0.988，当前数据不支持该焦亡关联。RNA正相关不能证明酶活抑制或抗焦亡。已发表IL1R2-ENO1机制来自巨噬细胞，迁移到Treg待验证（PMID40704655）。')]
for t,b in items:p(t);p(b)
page();p('重点候选当前结论（2/2）',title)
for t,b in [('FOXP1 / Core及Treg','状态重塑。Core PyroExec联合P=0.0219，等权状态标准化后P=0.494；损伤关联中的组成解释须单列。RNA与自身AUCell不概括为稳定同向。两个状态分开，不拼接为连续链。'),('DDIT4 / ActCD4','文献驱动的应激/损伤假说。PyroExec P=0.751，AllMac TNF-DDIT4 P=0.658，当前数据不支持这些目标关联。外部髓核细胞ROS/TXNIP/NLRP3实验（PMID38771134）不等于本项目完整链支持。'),('GSN / ActCD4','继续备选。ApoExec P=0.264及来源竞争问题保留。PROS1-AXL功能入口撤回，其他既有模型入口保留资格，不视为巨噬来源特异功能证明。'),('AHNAK / Naive；ST6GALNAC3 / ActCD4','AHNAK保留T细胞表型/OXPHOS线索，OXPHOS联合P=0.0417、调整P=0.0833。ST6GALNAC3保留MitoROS响应转录线索，联合P=0.00972、FDR=0.0972，调整P=0.222。均不升级为已经观察到线粒体损伤。')]:p(t);p(b)
p('完整疾病DE、状态、损伤、入口、外部方向和未决问题见R05_current_focus_conclusions.tsv；90基因资格见R05_current_all90.tsv。所有名义信号保留原家族FDR，不在重点子集重算。',small)
p('复核边界：原始细胞矩阵和大型网络对象未重复打包；未完成蛋白功能、目标细胞损伤实测及连续因果验证。后续按R01-R07返回已解决/部分解决/未解决，不将此完成记录当作独立审核通过。',small)
pdf=D/'v19_定向修订与当前结论.pdf'
def footer(c,d):c.setFont('CN',8);c.drawString(45,26,'v19 | 内部复核材料 | 2026-10-17');c.drawRightString(550,26,str(d.page))
SimpleDocTemplate(str(pdf),pagesize=(595.28,841.89),leftMargin=45,rightMargin=45,topMargin=42,bottomMargin=45).build(story,onFirstPage=footer,onLaterPages=footer)
cp(pdf,pdf.name)
manifest=[dict(path=str(p.relative_to(K)).replace('\\','/'),bytes=p.stat().st_size,sha256=sha(p)) for p in K.rglob('*') if p.is_file() and p.name!='MANIFEST_SHA256.tsv'];pd.DataFrame(manifest).to_csv(K/'MANIFEST_SHA256.tsv',sep='\t',index=False,encoding='utf-8')
z=D/'v19_独立复核材料.zip'
with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as zz:
 for p in K.rglob('*'):
  if p.is_file():zz.write(p,str(p.relative_to(K)))
with zipfile.ZipFile(z) as zz:assert zz.testzip() is None and not any(n.endswith('.zip') for n in zz.namelist())
info={'pdf':str(pdf),'zip':str(z),'pdf_sha256':sha(pdf),'zip_sha256':sha(z),'files':len(manifest)+1,'zip_MB':z.stat().st_size/1e6};(O/'delivery.json').write_text(json.dumps(info,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(info,ensure_ascii=False))