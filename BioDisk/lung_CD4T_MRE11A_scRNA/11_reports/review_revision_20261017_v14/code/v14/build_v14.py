from pathlib import Path
import json,re,html,shutil,hashlib
import pandas as pd
from bs4 import BeautifulSoup,Tag
import fitz
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,PageBreak,Image,Table,TableStyle,KeepTogether
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors
W=Path('C:/Users/31937/lung_review_20261017_v14');O=W/'revision_v14';T=W/'render';T.mkdir(exist_ok=True);Q=O/'supplement_v14';BASE=W/'baseline_v13'
soup=BeautifulSoup((BASE/'report_v13.html').read_text(encoding='utf-8'),'html.parser');secs=soup.find_all('section');cover=secs[0];sections=secs[2:];changes=[]
def sec(prefix):return next(s for s in sections if s.find(['h1','h2']).get_text().startswith(prefix))
def para(s,text):
 p=soup.new_tag('p');p.string=text;s.append(p);return p
def table(s,heads,rows):
 t=soup.new_tag('table');r=soup.new_tag('tr')
 for h in heads:c=soup.new_tag('th');c.string=str(h);r.append(c)
 t.append(r)
 for row in rows:
  r=soup.new_tag('tr')
  for v in row:c=soup.new_tag('td');c.string=str(v);r.append(c)
  t.append(r)
 s.append(t);return t
def add_after(prefix,title,body):
 ss=soup.new_tag('section');ss['id']='v14-'+str(len(sections));ss['data-chapter']='no';h=soup.new_tag('h2');h.string=title;ss.append(h)
 for p in body:para(ss,p)
 sections.insert(sections.index(sec(prefix))+1,ss);return ss
def replace(old,new):
 hits=[p for s in [cover]+sections for p in s.find_all('p') if old in p.get_text()]
 assert len(hits)==1,(old,len(hits));p=hits[0];before=p.get_text();p.string=before.replace(old,new);changes.append({'before':before,'after':p.get_text()})
cover.find('p').string='独立评审后补充修订版 v14｜2026年10月17日｜待再次独立评审'
replace('STAT1与GSDMD的转录关联在多个CD4范围中方向一致。','STAT1与GSDMD的转录关联在多个CD4范围中方向一致，但仅固定组成模型通过相应BH校正，且低丰度GSDMD的CD4内源性仍需确认。')
replace('接收端调控明确后','确定接收端候选后')
replace('固定组成后的Core STAT1–GSDMD关系为rho=1，P=0.002778、BH=0.02778，在该补充检验家族中通过BH校正；控制疾病组别后P=0.05556、BH=0.5556。三个算法种子方向一致。完整基因分析及组成结果见表S09–S10。','固定组成后的Core STAT1–GSDMD关系为rho=1，P=0.002778、BH=0.02778，仅该补充模型通过其完整家族BH。固定组成后再控制疾病组别，P=0.05556、BH=0.5556；原始组成直接控制组别则为rho=0.75、P=0.27778、BH=0.5556。不能排除共同疾病状态驱动。三个算法种子一致属于计算稳定性，Core、ActCD4、Treg来自同一批六只小鼠，均不构成独立验证队列。完整家族与输入见S09–S10、S21–S22。')
# Resolution is a property of the specified exact test, not a new test or threshold.
para(sec('1.8 '),'上述3对3双侧精确检验的20种分配包含标签完全对调的成对结果；绝对均值差相同，故理论最小P为2/20=0.10。使用原始检验规则及P<0.05口径，不更换检验或采用单侧P。')
para(sec('2 '),'这些方向和效应量保留为探索性结果。当前3对3双侧精确检验最小可达P=0.10，未达P<0.05不能单独解释为不存在变化；本限制只适用于该枚举组间检验，不适用于edgeR模型检验，也不适用于720种配对的样本相关检验。')
para(sec('4.4 '),'GSDMD来源核查直接影响本节解释：Saline_2的424个Core CD4细胞全部检出SFTPC和SCGB1A1；31个GSDMD阳性细胞均仅1 UMI。严格异源标记排除后，该样本剩余0个细胞，原六样本关联敏感性无法评价。患者中3个STAT1/GSDMD双阳性效应CD4也均伴髓系混合标记。结果保留为低丰度转录关联，尚不能确认纯CD4内源性GSDMD调控。去除背景基因或改变归一化分母不等于确定GSDMD来源。原始计数与保留样本表见S22。')
# Explicit four-model table follows the focal result, values are read from recomputation.
key=pd.read_csv(Q/'S21_key_exact_P_recalculation.tsv',sep='\t');labels={('Original','Pooled'):'原始组成 / 合并',('Original','ConditionAdjusted'):'原始组成 / 调整组别',('Fixed','Pooled'):'固定组成 / 合并',('Fixed','ConditionAdjusted'):'固定组成 / 调整组别'}
table(sec('4.4 '),['Core模型（均n=6）','rho','P','原完整家族BH'],[[labels[(r.composition,r.model)],f'{r.rho:.4f}',f'{r.P:.6f}',f'{r.BH:.6f}'] for r in key.itertuples()])
m=add_after('1.8 ','1.9 评审后限定补充分析与统计核对',[
'新增分析在读取结果前冻结六条配体轴、两个巨噬细胞范围、四个CD4范围和seed17。IL15对应STAT5A/B；TNF对应NFKB1/2、RELA/B；IL18对应NFKB1、RELA、JUN/JUNB、FOS/FOSB；IL27、IFNB1和IFNG对应STAT1。每个来源、接收范围和模型构成一个包含15项计划配对的家族，仅对可评估P进行BH校正；常量或缺少主种子regulon保留NA及原因。合并与组别调整分别沿用720及36种精确配对。',
'复合物资格重新读取冻结原始计数：IL15发送端要求IL15与IL15RA同细胞检出，IL27发送端要求IL27与EBI3；接收端分别核对IL2RB/IL2RG、IL18R1/IL18RAP、IL27RA/IL6ST、IFNAR1/IFNAR2及IFNGR1/IFNGR2。TNF采用TNFRSF1A或TNFRSF1B任一阳性规则。每个独立样本保留细胞分母、阳性数和比例。',
'全肺竞争来源复用GSE314943冻结Broad标签，以原始计数和原始总UMI计算各样本、各标签的CXCL16与IFNG平均CP10k及检出率。跨样本组均值对供者等权；未对细胞进行推断检验。原Broad标签与精细AMlike/AllMac来自不同冻结范围，分别记录。',
'基因集重叠核查按历史canon别名映射重建DropOverlap。定向敏感性从既有匹配对照中剔除对应STAT1或IKZF2调控子成员及TF自身，保留原模块成员、归一化、样本和细胞，不重新匹配对照。七个程序在每个TF、范围和模型内构成新敏感性BH家族；固定组成单独计算。原始评分仅作重现参考，不以缩小家族重标其BH。'])
num=json.loads((O/'provenance/v14/local_numeric_validation.json').read_text(encoding='utf-8'))
para(m,f'对{num["BH_families_checked"]}个已保存完整检验家族独立复算BH，并从六样本输入复算Core STAT1–GSDMD四种精确检验。新增S21–S27保存完整家族、样本名单、原始统计、重叠成员、全部计划配对及未评估项。所有新分析属于评审后定向探索性补充。')
ov=add_after('4.4 ','4.5 去除重叠匹配对照后，IKZF2关联仍依赖组成',[
'原DropOverlap已排除全部预设功能模块成员，GSDMD也不在STAT1调控子中；但匹配对照并未同步去重叠。seed17中，STAT1与MitoROS对照重叠13个基因，IKZF2与ApoExec对照重叠15个。该构造需要单独进行评分敏感性检查。',
'重建的336个原始样本评分与冻结值一致后，按1.9节去除重叠对照。STAT1–MitoROS在Core中仍为rho=0.9429、P=0.01667，新敏感性家族BH=0.05833；固定组成后仍为rho=0.8857、P=0.03333、BH=0.23333。固定组成并调整组别后P=0.16667。其方向保留，证据继续定位为探索性关联。',
'IKZF2–ApoExec在Core中变为rho=0.9429、P=0.01667、BH=0.11667；固定组成后变为rho=0.0286、P=1、BH=1。因此，当前结果不支持把IKZF2作为直接调控凋亡的主证据。该敏感性未改写原S08及图9，也不涉及以单基因CP10k为终点的STAT1–GSDMD检验。',
'S25列出三个种子全部重叠成员，S27列出全部新评分、原值重现与相关检验。广义GO自噬集合与regulon的重叠另行列出；它不在原七程序TF关联检验家族内，不将该重叠误记为现有TF检验循环证据。'])
para(sec('5.2 '),'配体上调与TF关联并非同一发送端内的连续证据：AMlike支持CXCL16名义上调但与ActCD4 STAT1关联P=1；AllMac的对应关联P=0.002778但疾病差异P=0.149。两条线索分别支持来源表达和接收端候选，不拼接为已建立的因果链。')
para(sec('5.1.1 '),'CXCL16对ActCD4上调目标的来源内排名为AMlike第199、AllMac第280；下调目标为第154及216。其优先级依据为来源表达、CXCR6覆盖和可验证的募集/滞留问题，不能称为多种通讯方法一致支持的首选驱动配体。完整来源排名见S19。')
prox=add_after('5.3 ','5.4 受体与近端响应限定各条候选通路',[
'本次240项计划配对中200项可评估，40项保留NA；7项原始P<0.05，均未通过相应家族BH。该汇总包括同一批小鼠的重复范围和模型，不能当作240次独立实验。以下逐轴说明同一发送端与接收范围中的证据及缺口。'])
table(prox,['候选轴','原始计数与近端结果','当前定位'],[
['IL15','病例AMlike中IL15/IL15RA同细胞为4/72、3/94、1/67；ActCD4 IL2RB/IL2RG为116/284、48/100、90/225。主种子17缺少STAT5A/B调控子，32项计划关联记NA。','具有转录组分支持；跨呈递蛋白及STAT5活化未确认，不把NA当通过。'],
['TNF','ActCD4任一TNF受体为92/284、34/100、54/225。AllMac TNF–ActCD4 NFKB2 rho=0.8857，P=0.03333、BH=0.43333；Naive的对应组别调整P=0.02778、BH=0.24074。','状态特异的名义近端线索；不是STAT1/IKZF2连接或蛋白活化证据。'],
['IL18','ActCD4完整受体为20/284、11/100、33/225。AMlike IL18–ActCD4 NFKB1为负相关，rho=-0.8857、P=0.03333、BH=0.4000；AP1计划关联未达P<0.05。','保留负向结果；当前数据不支持正向IL18驱动该近端响应链。'],
['IL27','ActCD4完整受体为14/284、9/100、12/225；AMlike病例IL27/EBI3同细胞均为0，AllMac为0/388、2/358、3/259。对应ActCD4 STAT1合并P分别0.41944和0.24167。','接收端转录组分存在；完整发送端和连续关联支持不足。'],
['IFNB1','ActCD4 IFNAR1/IFNAR2为14/284、3/100、19/225。AMlike IFNB1为常量，8项记NA；AllMac–ActCD4 STAT1合并P=0.61667。','低丰度/常量限制来源评价；未建立定向STAT1支持。'],
['IFNG','ActCD4 IFNGR1/IFNGR2为43/284、14/100、28/225。巨噬来源与ActCD4 STAT1未达P<0.05；全肺T/NK表达分布见5.5节。','T/NK来源为竞争解释，未确认其驱动STAT1。']])
para(prox,'S23保留全部正向、负向、名义及NA结果，S24提供逐样本完整复合物计数。SCENIC在此提供近端转录响应代理，不能代替磷酸化、受体占据、IL15反式呈递或分泌蛋白测量。尚缺IL18后的IFNG诱导时序及这些近端信号到IKZF2/STAT1的同细胞功能连接，不以已运行NicheNet或SCENIC替代。')
srcsec=add_after('5.4 ','5.5 CXCL16具有多来源表达，T/NK提供IFNG竞争解释',[
'冻结全肺Broad标签的逐小鼠表达比较显示，CXCL16并非只见于巨噬细胞。下表为病例三个独立小鼠的等权平均CP10k；数值反映测得转录表达，不等同组织分泌贡献或来源特异的功能效应。完整逐样本细胞数、检出率及全部标签见S26。'])
src=pd.read_csv(Q/'S26_whole_lung_competing_sources.tsv',sep='\t');rr=[]
for g,labs in [('CXCL16',['Mac','DC','AT1','AT2Like','Fib']),('IFNG',['NKLike','T','Mac'])]:
 for label in labs:
  a=src[(src.gene==g)&(src.source==label)&(src.group=='Silica')];rr.append([g,label,str(len(a)),f'{a.mean_CP10k.mean():.4f}',f'{a.positive_pct.mean():.2f}%'])
table(srcsec,['基因','冻结Broad标签','独立样本数','等权平均CP10k','等权平均检出率'],rr)
para(srcsec,'CXCL16在Mac、DC及AT1等标签中均有表达；因此应把后续问题限定为特定巨噬亚群在特定阶段是否具有额外功能贡献，并以空间定位或来源特异干预与上皮/DC等比较。IFNG在NKLike和T标签中高于Mac，T/NK是接收端STAT1的合理竞争来源；此处尚未建立配体分泌、来源到接收端的空间接触或因果联系。Broad标签沿用旧注释，未在本轮重新做谱系注释审核；少细胞标签的数值需结合S26分母使用。')
para(sec('7 '),'创新性收窄：You等2024年已报告硅尘相关CD4组织驻留记忆T细胞及CXCR6等滞留标志；Wein等2019年在肺CD8背景中研究了气道上皮CXCL16与CXCR6定位。[2,11] 因而，CXCR6与肺部T细胞募集/滞留有关不是本项目的新发现。当前待检验的新问题是：特定矽肺阶段、特定巨噬亚群的CXCL16，是否对特定CD4状态具有来源可区分的募集/滞留效应，以及该效应是否连接后续细胞内表型。')
para(sec('8 '),'小样本检验与模型边界：本项目主要推断来自六只小鼠，按探索性证据解释。模块/AUC组间枚举检验的P下限为0.10，不以其未显著抹去方向和效应。样本关联的原始组成、固定组成及组别调整具有不同估计目标，仅在特定补充家族通过BH的结果不能概括为跨模型稳定显著。去重叠对照敏感性与竞争来源分析已补足可计算部分，但不能修复环境RNA来源不明或替代蛋白与功能验证。')
para(sec('9 '),'评审后补充进一步限定了候选的用途：STAT1–GSDMD的转录关联、CXCL16–CXCR6的表达定位、TNF近端名义信号与T/NK来源IFNG竞争解释分别记录。当前结论仍是形成可继续验证的候选方向，尚未建立完整巨噬细胞—CD4损伤/焦亡机制。')
t=sec('9 ').find('table')
for code,desc in [('S21','完整BH家族、106家族复算、四模型精确P及输入来源'),('S22','六样本名单、背景身份、GSDMD UMI及原模型输入'),('S23','240项计划近端配对、原始AUC疾病效应、可评估性'),('S24','发送端与接收端完整复合物逐样本共检出'),('S25','三个种子regulon、模块与匹配对照重叠明细'),('S26','CXCL16与IFNG全肺竞争来源逐样本描述'),('S27','匹配对照去重叠敏感性：全部样本评分、成员与统计')]:
 tr=soup.new_tag('tr')
 for v in [code,desc]:td=soup.new_tag('td');td.string=v;tr.append(td)
 t.append(tr)
# Two literature entries already existed; preserve their precise primary-source links.
for p in sec('参考').find_all('p'):
 if p.get_text().startswith('[2]'):p.clear();p.append('[2] You等，2024。Crystalline silica-induced recruitment and immuno-imbalance of CD4+ tissue resident memory T cells promote silicosis progression. ');a=soup.new_tag('a',href='https://www.nature.com/articles/s42003-024-06662-z');a.string='原始研究';p.append(a)
 if p.get_text().startswith('[11]'):p.clear();p.append('[11] Wein等，2019。CXCR6 regulates localization of tissue-resident memory CD8 T cells to the airways. J Exp Med. DOI:10.1084/jem.20181308. ');a=soup.new_tag('a',href='https://rupress.org/jem/article/216/12/2748/132521/CXCR6-regulates-localization-of-tissue-resident');a.string='原始研究';p.append(a)
# Put focal risk and the four-model table on a named continuation page.
focal=sec('4.4 ');risk=next(p for p in focal.find_all('p') if p.get_text().startswith('GSDMD来源核查'));model_table=focal.find_all('table')[-1]
continuation=add_after('4.4 ','4.4.1 GSDMD来源风险与四种模型分别评价',[]);continuation.append(risk.extract());continuation.append(model_table.extract())
# Give the long supplement index its own page, with a short first column.
conclusion=sec('9 ');index=conclusion.find('table');idxsec=add_after('9 ','9.1 完整补充统计与复核入口',[]);idxsec.append(index.extract())
para(idxsec,'S01–S20保留于supplement；S21–S27位于supplement_v14。附表完整保留原始P、BH、样本范围与NA，不只提供显著子集。')
nav=BeautifulSoup('<section><h1>阅读路径</h1><p>数据与方法 → 功能表型 → 候选 → TF与来源核查 → 通讯与竞争来源 → 人类证据 → 验证定位</p></section>','html.parser').section
for ss in sections:
 if ss.get('data-chapter')=='yes':
  p=soup.new_tag('p');a=soup.new_tag('a',href='#'+ss['id']);a.string=ss.find(['h1','h2']).get_text();p.append(a);nav.append(p)
para(nav,'v14保留原18幅图与S01–S20数值，新增1.9、4.5、5.4、5.5节及S21–S27。逐项评审回应另附；本版等待再次独立评审。')
css=soup.find('style').get_text();(O/'report_v14.html').write_text('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>矽肺评审后补充报告v14</title><style>'+css+'</style>'+str(cover)+str(nav)+''.join(map(str,sections))+'</html>',encoding='utf-8')
fig_no={Path(img['src']).name:int(re.match(r'图(\d+)',img.find_next_sibling('p',class_='caption').get_text()).group(1)) for ss in sections for img in ss.find_all('img')}
engine=Path('C:/codex_test/tmp/silicosis_v10/build_integrated.py').read_text(encoding='utf-8').split('# Generate a continuous typeset report')[1].split('# Content conservation checks')[0].split('\n',1)[1].replace('v10','v14')
engine=engine.replace('{1:300,2:300,3:240,4:300,5:140,8:200,10:160,14:250,15:250}','{1:230,2:230,3:300,4:240,5:300,6:140,9:295,11:160,14:260,15:290,17:250,18:250}')
engine=engine.replace('if nc==3: widths=', 'if nc==2: widths=[90,width-90]'+chr(10)+'   elif nc==3: widths=')
engine=engine.replace("('TOPPADDING',(0,0),(-1,-1),5)","('TOPPADDING',(0,0),(-1,-1),3 if nc==2 else 5)").replace("('BOTTOMPADDING',(0,0),(-1,-1),5)","('BOTTOMPADDING',(0,0),(-1,-1),3 if nc==2 else 5)")
exec(engine,globals())
(O/'provenance/v14/page_map.json').write_text(json.dumps(bookmarks,ensure_ascii=False,indent=2),encoding='utf-8');(O/'provenance/v14/text_replacements.json').write_text(json.dumps(changes,ensure_ascii=False,indent=2),encoding='utf-8')
d=fitz.open(O/'silicosis_integrated_report_v14.pdf');assert len(fig_no)==18
for i,p in enumerate(d):
 assert '\ufffd' not in p.get_text();assert all(b[0]>=0 and b[1]>=0 and b[2]<=842 and b[3]<=595 for b in p.get_text('blocks')),(i,p.get_text())
 p.get_pixmap(matrix=fitz.Matrix(1,1)).save(T/f'page_{i+1:02d}.png')
print('MAIN_PDF_PAGES',len(d));print('PAGE_MAP',json.dumps(bookmarks,ensure_ascii=False))