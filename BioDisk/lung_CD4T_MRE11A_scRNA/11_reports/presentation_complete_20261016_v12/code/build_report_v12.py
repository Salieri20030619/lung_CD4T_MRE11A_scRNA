from pathlib import Path
import json,re,html,shutil
from collections import Counter
import pandas as pd
from bs4 import BeautifulSoup,Tag
import fitz
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,PageBreak,Image,Table,TableStyle,KeepTogether
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors
R=Path('C:/codex_test');O=R/'output/silicosis_v12';T=R/'tmp/silicosis_v12';OLD=R/'output/silicosis_v11'
hs=BeautifulSoup((OLD/'report_v11.html').read_text(encoding='utf-8'),'html.parser');cover=hs.find_all('section')[0];sections=hs.find_all('section')[2:];cover.find('p').string='评分与NicheNet展示完善版 v12｜2026年10月16日'
def find(prefix):return next(s for s in sections if s.find(['h1','h2']).get_text().startswith(prefix))
def sec(text):return BeautifulSoup(text,'html.parser').section
methods=sec('''<section id="section-score-nichenet-methods" data-chapter="no"><h2>1.2 死亡程序评分与NicheNet预测的具体定义</h2>
<h3>表达匹配的模块评分</h3><p>每个细胞先按原始总UMI计算CP10k并取log1p。模块评分 = 模块基因平均log1p(CP10k) − 匹配对照基因平均log1p(CP10k)。随后在每个独立样本、每个CD4范围内取细胞评分均值，组均值再对样本等权计算。该方法是项目预设的自定义表达匹配评分。</p>
<table><tr><th>程序</th><th>完整基因集</th><th>主队列覆盖/对照基因数</th></tr><tr><td>焦亡相关 PyroExec</td><td>PYCARD、CASP1、GSDMD、IL1B、IL18</td><td>5/5成员可测且有表达；100个匹配对照</td></tr><tr><td>凋亡执行 ApoExec</td><td>BAX、BAK1、APAF1、CASP9、CASP3、CASP7、CYCS</td><td>7/7成员可测且有表达；140个匹配对照</td></tr><tr><td>炎症小体 Inflam</td><td>NLRP3、AIM2、NLRC4、PYCARD、CASP1</td><td>5/5成员可测且有表达；100个匹配对照</td></tr></table>
<p>参考表达取各样本平均表达的等权均值，按非零参考表达划分25个分位箱。对每个模块基因从同箱选择最多20个表达最接近的对照；同箱不足20个时从全部合格池选取，模块内去重。合格池排除历史预设候选IL1R2、PENK、TTC39C、FOXP1、TNFRSF9、LAG3、MRE11，以及全部预设程序成员、上皮/髓系背景标记、线粒体编码、核糖体和免疫球蛋白基因。主分析与BG1敏感性沿用同一对照集。全部成员、对照和覆盖见S18。</p>
<p>3对3组间评分检验枚举20种标签分配，以绝对组均值差为统计量；BH分别用于Core三项死亡程序和三个亚群的九项死亡程序，All与BG1各自校正。保留全部非空样本，未以细胞数代替生物学样本量。</p>
<h3>NicheNet输入、指标与来源排名</h3><p>NicheNet复用小鼠v2配体-靶基因先验。各CD4范围的背景为受检基因、先验靶基因和表达输入的交集；上/下调目标分别按原始P&lt;0.05及log2FC正/负定义。活动值由接收端目标集预测，随后按AMlike或AllMac来源表达资格过滤并重新排名。corrected AUPR为原AUPR减去目标数/背景数；完整输入、活动值和排名见S19。</p></section>''')
sections.insert(sections.index(find('1.1'))+1,methods)
pyro=sec('''<section id="section-pyroptosis" data-chapter="yes"><h1>2 焦亡相关评分在各CD4范围均向上变化</h1>
<p>Core、Naive、ActCD4和Treg的焦亡相关转录评分均为硅尘组较高，但样本级检验未达到名义显著。各范围差值依次为+0.05643、+0.04488、+0.05163和+0.03957；原始P依次为0.20、0.40、0.50和0.60，BH为0.30、1、1和1。</p>
<img src="figures/V12_01_pyroptosis_scores.png"/><p class="caption">图1 焦亡相关转录评分。蓝色C1-C3为对照，红色S1-S3为硅尘，每点为独立小鼠，横线为组均值；四面板使用相同纵轴。统计单位为3对3小鼠，Core与三个亚群不是独立队列。</p>
<p>这五个基因构成的是焦亡相关转录程序。当前数据不支持将任一亚群称为疾病中显著焦亡最强的亚群。炎症小体评分、全部样本值及BG1敏感性保留在S18；后续GSDMD逐基因分析与此模块结果分别解释。</p></section>''')
old=find('2 ');sections[sections.index(old)]=pyro
apo=sec('''<section id="section-apoptosis" data-chapter="no"><h2>2.1 凋亡执行评分的变化方向随CD4状态不同</h2>
<p>凋亡执行相关评分在Core中升高，差值+0.00957、P=0.10、BH=0.30；ActCD4差值约+0.000056、P=1、BH=1。Naive与Treg则分别下降−0.00758和−0.01056，P为0.30和0.70，BH均为1。当前数据不支持一致的疾病相关凋亡程序升高。</p>
<img src="figures/V12_02_apoptosis_scores.png"/><p class="caption">图2 凋亡执行相关转录评分。七基因集合、评分单位和样本汇总方式与1.2节一致；各点为样本均值，未将细胞作为独立重复。所有范围均使用相同纵轴。</p>
<p>Core均值与亚群内变化并不完全一致，因此4.2节进一步检查IKZF2与凋亡程序的关联是否受到亚群组成影响。凋亡评分的组间变化和TF-程序样本关联是两个不同问题。</p></section>''');sections.insert(sections.index(pyro)+1,apo)
find('2.1 线粒体').find('h2').string='2.2 线粒体与自噬程序呈现亚群特异方向'
comp=sec('''<section id="section-composition-scatter" data-chapter="no"><h2>4.2 STAT1与IKZF2的关联对组成调整反应不同</h2>
<p>每一行比较同一TF-程序关系，左列使用实际亚群比例，右列统一为六只小鼠的平均亚群比例。蓝色为对照、红色为硅尘；两列保留同一批样本编号及相同坐标范围，直接展示组成调整改变了哪些样本的位置。</p>
<img src="figures/V12_09_composition_scatter.png"/><p class="caption">图9 实际比例与统一比例下的样本级关联。每个范围均n=6，图中ρ、P、BH沿用原结果；统一权重为Naive 46.66%、ActCD4 31.00%、Treg 22.34%。坐标和精确统计见S18。</p>
<p>Core IKZF2-ApoExec的ρ由1.0000（P=0.002778、BH=0.2743）变为−0.7714（P=0.1028、BH=0.5688）。Core STAT1-MitoROS由0.9429（P=0.01667、BH=0.3678）变为0.8857（P=0.03333、BH=0.4467）。前者方向改变，后者保留正向联系。</p>
<p>Treg STAT1-MitoROS原关联为ρ=0.8857、P=0.03333、BH=0.3790。另行控制疾病组别后，Core固定组成的STAT1-MitoROS P=0.1667；这一模型与图中的组成统一不同。完整关系和种子敏感性见S08。</p></section>''');old=find('4.2');sections[sections.index(old)]=comp
nn_activity=sec('''<section id="section-nichenet-activity" data-chapter="no"><h2>5.1.1 NicheNet预测活性为配体排名提供定量依据</h2>
<p>在AMlike来源表达资格内，ActCD4的上调目标集包括182个基因、背景5236个。TNF与IL15的来源内排名分别为4和6，corrected AUPR分别为0.02043和0.01916。Naive、Treg使用各自的目标和背景集合，图中同时展示前列配体及实际指标。</p>
<img src="figures/V12_14_nichenet_activity.png"/><p class="caption">图14 NicheNet预测活动值。AMlike中每个亚群取来源内排名前8；目标集合、背景数及corrected AUPR均来自冻结结果。红色为现有重点配体。完整Core/亚群、上/下调及AMlike/AllMac来源排名见S19。</p>
<p>预测排名回答“哪些先验配体能较好解释接收端差异基因集合”。来源表达支持、受体完整性、配体疾病差异和同样本关联仍由不同证据分别评价；预测活动值不作为生物学效应的P值。</p></section>''')
nn_target=sec('''<section id="section-nichenet-death-targets" data-chapter="no"><h2>5.1.2 NicheNet先验与死亡程序基因的连接可逐项核对</h2>
<p>将既有配体-靶基因先验映射到五个焦亡相关和七个凋亡执行基因，并叠加ActCD4的差异表达。12个模块基因中9个进入既有ActCD4差异背景，只有GSDMD进入其上调目标集合。IL18、CASP9和CASP7未进入该差异背景。</p>
<img src="figures/V12_15_nichenet_death_targets.png"/><p class="caption">图15 配体-死亡程序基因的原始先验权重与ActCD4差异。配体取既有AMlike/ActCD4上调排名前5并加入重点配体，未按死亡基因权重筛选。条带标记目标/背景资格，灰色差异格为未受检。原始P、BH和全部权重见S19。</p>
<p>本图补足既有预测的模块对应关系，没有将焦亡或凋亡基因集重新作为目标拟合NicheNet。先验无激活/抑制符号，因此权重较高不直接表示促焦亡、促凋亡或本队列的因果作用。</p></section>''')
idx=sections.index(find('5.1 '))+1;sections[idx:idx]=[nn_activity,nn_target]
last=find('9 ');tb=last.find('table')
for label,desc in [('S18','焦亡/凋亡评分定义、匹配对照、样本值与固定组成散点坐标'),('S19','NicheNet完整预测指标、目标/背景成员及配体-死亡程序基因先验权重')]:
 tr=hs.new_tag('tr')
 for text in [label,desc]:td=hs.new_tag('td');td.string=text;tr.append(td)
 tb.append(tr)

# Number figures in reading order and preserve the source filenames for traceability.
fig_no={};idx=0
for ss in sections:
 for img in ss.find_all('img'):
  idx+=1;fig_no[Path(img['src']).name]=idx;cap=img.find_next_sibling('p',class_='caption');assert cap is not None
  cap.string=re.sub(r'^图\d+\s*',f'图{idx}  ',cap.get_text())
assert idx==18 and fig_no['V12_09_composition_scatter.png']==9 and fig_no['V12_15_nichenet_death_targets.png']==15
nav=sec('<section><h1>阅读路径</h1><p>数据与评分定义 → 焦亡/凋亡表型 → 候选 → TF调控 → NicheNet与通讯 → 人类证据 → 实验选择</p></section>')
for ss in sections:
 if ss.get('data-chapter')=='yes':
  p=hs.new_tag('p');link=hs.new_tag('a',href='#'+ss['id']);link.string=ss.find(['h1','h2']).get_text();p.append(link);nav.append(p)
p=hs.new_tag('p');p.string='重点更新：1.2节给出评分与NicheNet定义；第2节分别展示焦亡与凋亡；4.2节改为样本散点；5.1.1-5.1.2节展开NicheNet活动值和靶基因。S17候选选择依据继续保留。';nav.append(p)
css=hs.find('style').get_text();(O/'report_v12.html').write_text('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>矽肺研究报告 v12</title><style>'+css+'</style>'+str(cover)+str(nav)+''.join(map(str,sections))+'</html>',encoding='utf-8')

engine=(R/'tmp/silicosis_v10/build_integrated.py').read_text(encoding='utf-8').split('# Generate a continuous typeset report')[1].split('# Content conservation checks')[0].split('\n',1)[1].replace('v10','v12')
oldmap='{1:300,2:300,3:240,4:300,5:140,8:200,10:160,14:250,15:250}'
newmap='{1:230,2:230,3:300,4:240,5:300,6:140,9:295,11:160,14:260,15:290,17:250,18:250}'
assert oldmap in engine;engine=engine.replace(oldmap,newmap)
exec(engine,globals())
final=BeautifulSoup((O/'report_v12.html').read_text(encoding='utf-8'),'html.parser');nums=[int(re.match(r'图(\d+)',p.get_text()).group(1)) for p in final.find_all('p',class_='caption')];assert nums==list(range(1,19))
d=fitz.open(O/'silicosis_integrated_report_v12.pdf');text='\n'.join(p.get_text() for p in d)
for key in ['−0.7714','0.5688','0.4467','0.1667','182','5236','0.02043','0.01916','140个','100个']:
 assert key in text,key
for i,p in enumerate(d):
 assert '\ufffd' not in p.get_text();assert all(b[0]>=-1 and b[2]<=843 and b[1]>=0 and b[3]<=596 for b in p.get_text('blocks')),(i,p.get_text())
 p.get_pixmap(matrix=fitz.Matrix(1,1)).save(T/f'page_{i+1:02d}.png')
(O/'provenance/page_map.json').write_text(json.dumps(bookmarks,ensure_ascii=False,indent=2),encoding='utf-8')
(O/'provenance/presentation_review_resolution.json').write_text(json.dumps({'R01':{'status':'resolved_in_presentation_scope','changes':'NicheNet metrics, source ranks, target/background definitions and prior-to-death-gene mapping added; no death-target model refit claimed.'},'R02':{'status':'resolved','changes':'Pyroptosis and apoptosis separately plotted and narrated; full panels, matched controls, scoring formula and sample aggregation supplied.'},'R03':{'status':'resolved','changes':'Section 4.2 coefficient-line plot replaced by four sample-coordinate scatter panels, same scales and sample IDs.'},'new_inferential_tests':0,'figures':18,'pages':len(d)},ensure_ascii=False,indent=2),encoding='utf-8')
shutil.copy2(__file__,O/'code/build_report_v12.py');print('REPORT_COMPLETE pages',len(d));print(json.dumps(bookmarks,ensure_ascii=False))
