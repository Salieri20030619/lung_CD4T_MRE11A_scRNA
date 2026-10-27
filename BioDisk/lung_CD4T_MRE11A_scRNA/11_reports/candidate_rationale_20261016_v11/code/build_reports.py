from pathlib import Path
import json,re,html,shutil,hashlib
from collections import Counter
import pandas as pd
from bs4 import BeautifulSoup,Tag
import fitz
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,PageBreak,Image,Table,TableStyle,KeepTogether
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors
R=Path('C:/codex_test');T=R/'tmp/silicosis_v11';O=R/'output/silicosis_v11';OLD=R/'output/silicosis_v10'
for folder in ['figures','assets','plot_data','source_tables']:
 shutil.copytree(OLD/folder,O/folder,dirs_exist_ok=True)
for p in (OLD/'supplement').iterdir():shutil.copy2(p,O/'supplement'/p.name)
hs=BeautifulSoup((OLD/'report_v10.html').read_text(encoding='utf-8'),'html.parser');original_text=hs.get_text();original_cells=Counter(x.get_text() for x in hs.find_all(['td','th']));changes=[]
for p in hs.find_all('p'):
 text=p.get_text()
 if '每批审阅10个' in text:
  new='历史初筛以样本级差异表达结果建立候选池：对每张差异表按原始P<0.05和|log2FC|≥0.5，分别导出上调和下调前20个。随后结合亚群效应、样本一致性、文献创新性及背景负担敏感性确定重点候选。该效应门槛用于复原历史流程，不作为后续分析统一准入条件。'
 elif '分批文献评估与统计筛选' in text:new=text.replace('分批文献评估与统计筛选','文献评估、样本级核查与统计筛选')
 elif '逐批文献评估' in text:new=text.replace('逐批文献评估','候选文献评估')
 else:continue
 changes.append({'old':text,'new':new});p.string=new
cover=hs.find_all('section')[0];cover.find('p').string='候选选择依据补充版 v11｜2026年10月16日'
sections=hs.find_all('section')[2:]
insert=BeautifulSoup('''<section id="section-s17" data-chapter="no"><h2>3.3 高效应基因的选择还取决于细胞来源与机制适配</h2>
<p>效应量排序用于定位变化幅度，不能单独决定实验主靶点。补充表S17逐项说明Core和ActCD4名义差异绝对效应前20中、未列入11个既有重点基因的27个基因，并给出同状态log2FC、原始P、原BH、低于其效应的现有候选、文献发现与后续处置。</p>
<table><tr><th>依据类型</th><th>可核实内容</th><th>本报告的处理</th></tr>
<tr><td>历史明确记录</td><td>TNFRSF9、LAG3已有相近疾病/功能研究，创新优先级下降；FOXP1作为辅助方向。IL1R2、PENK优先，TTC39C作为探索候选。</td><td>保留原候选的角色差异，不将“未作为创新主靶点”写成统计淘汰。</td></tr>
<tr><td>历史提出来源核查</td><td>IGKC、IGHG2B、IGHG2C、JCHAIN在前列差异结果中出现；原记录要求核查细胞来源。</td><td>不把高检出率直接等同高污染负担，也不声称已有最终排除结论。</td></tr>
<tr><td>本次补充文献评估</td><td>对27个高效应基因逐项列出已有研究的疾病、细胞和作用层次；没有历史决定的条目明确标注。</td><td>分为来源待核实、功能证据待补和保留比较候选等情况，不以“已有研究”作为统一排除标准。</td></tr></table>
<p>例如，ActCD4中NME2、CCL5和GSN的变化为下调，其绝对效应应与同为下调的FOXP1比较；不能因不符合上调叙事而略去。对于仍与当前问题匹配、尚无充分排除依据的基因，S17保留为比较候选。</p>
<p>完整效应比较覆盖四个CD4范围的1821个基因-状态条目；其中未逐项文献评估的条目明确标注。历史依据来自2026年9月14日原始筛选对话，S17另列记录编号H01-H06及来源，区分当时决定与本次解释。</p>
<p><a href="supplement/S17_candidate_selection_rationale.xlsx">补充表S17（Excel）</a>；<a href="supplement/S17_candidate_selection_rationale.pdf">补充表S17（阅读版）</a>。</p></section>''','html.parser').section
idx=next(i for i,s in enumerate(sections) if s.find(['h1','h2']).get_text().startswith('3.2'))+1
sections.insert(idx,insert)
conclusion=next(s for s in sections if s.find(['h1','h2']).get_text().startswith('9 '));tb=conclusion.find('table');tr=hs.new_tag('tr')
for txt in ['S17','高效应基因的历史筛选依据、逐基因文献评估及完整效应比较']:
 td=hs.new_tag('td');td.string=txt;tr.append(td)
tb.append(tr)
nav=BeautifulSoup('<section><h1>阅读路径</h1><p>数据设计 → CD4功能表型 → 候选与选择依据 → TF调控 → 巨噬细胞通讯 → 人类证据 → 实验选择与讨论</p></section>','html.parser').section
for s in sections:
 if s.get('data-chapter')=='yes':
  p=hs.new_tag('p');a=hs.new_tag('a',href='#'+s['id']);a.string=s.find(['h1','h2']).get_text();p.append(a);nav.append(p)
p=hs.new_tag('p');p.string='新增3.3节及补充表S17。图1-15和原统计结果保持不变，候选筛选历史表述根据原始对话修正。';nav.append(p)
css=hs.find('style').get_text();(O/'report_v11.html').write_text('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>矽肺研究报告 v11</title><style>'+css+'</style>'+str(cover)+str(nav)+''.join(map(str,sections))+'</html>',encoding='utf-8')
fig_no={Path(img['src']).name:int(re.match(r'图(\d+)',img.find_next_sibling('p',class_='caption').get_text()).group(1)) for s in sections for img in s.find_all('img')}
engine=(R/'tmp/silicosis_v10/build_integrated.py').read_text(encoding='utf-8').split('# Generate a continuous typeset report')[1].split('# Content conservation checks')[0]
engine=engine.split('\n',1)[1].replace('v10','v11');exec(engine,globals())

# Supplemental reading table: same 27 rows and sources as the workbook.
data=json.loads((T/'workbook_data.json').read_text(encoding='utf-8'));m=data['sheets'][0];records=[dict(zip(m['columns'],r)) for r in m['rows']]
hist=data['sheets'][1];refs=data['references'];ref_by_url={x['url']:x['id'] for x in refs}
small=ParagraphStyle('s17cell',fontName='CJK',fontSize=9,leading=13,wordWrap='CJK')
th=ParagraphStyle('s17head',parent=small,textColor=colors.white)
def P(text,style=small):return Paragraph(html.escape(str(text)).replace('\n','<br/>'),style)
def tbl(rows,widths):
 t=Table([[P(c,th if i==0 else small) for c in row] for i,row in enumerate(rows)],colWidths=widths,repeatRows=1,hAlign='LEFT')
 t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#183448')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f2f5f7')]),('VALIGN',(0,0),(-1,-1),'TOP'),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#183448')),('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6)]));return t
st=[Paragraph('补充表S17｜高效应基因的选择依据',styles['h1']),P('统计与文献比较截至2026年10月16日。详细评估覆盖27个独立基因；Excel另提供四个CD4范围的完整效应比较、候选基准和原始文献链接。'),Spacer(1,9),P('阅读规则：按同一状态、原始P<0.05及|log2FC|比较；另列同方向的低效应候选。仅复用冻结P和BH，不重算统计。表中“本次补充评估”不冒充历史逐基因排除决定。'),Spacer(1,12)]
st.append(tbl([['对象','历史明确记录','解释边界']]+[[r[0],r[1],r[3]] for r in hist['rows']],[125,365,276]))
st.append(Spacer(1,10));st.append(P('原始对话：数据集质量比较分析。筛选相关记录H01-H06保存在资料包中。历史上“检出率高不等于UMI负担高”的修正被保留。'))
for start in range(0,len(records),3):
 st.append(PageBreak());st.append(Paragraph(f'S17 逐基因说明（{start+1}-{min(start+3,len(records))} / 27）',styles['h1']))
 rows=[['基因及代表效应','文献具体发现与场景','未优先/保留依据','历史依据及来源']]
 for r in records[start:start+3]:
  vals=f"{r['基因']}\n{r['代表状态']}\nlog2FC={r['log2FC']:+.4f}\nP={r['原始P']:.4g}\nBH={r['原BH']:.4g}"
  finding=r['文献具体发现']+'\n场景：'+r['研究场景']
  reason=r['当前未优先/保留依据']+'\n当前建议：'+r['当前处置建议']
  histtxt=r['历史记录']+'\n文献编号：'+ref_by_url[r['来源链接']]
  rows.append([vals,finding,reason,histtxt])
 st.append(tbl(rows,[123,220,253,170]));st.append(Spacer(1,12));st.append(P('以上效应均为对应代表状态的原始值。完整状态比较、低效应候选名单和再次考虑条件见Excel。仅有其他细胞或疾病的文献，不足以判定矽肺CD4机制已被解决。'))
st.append(PageBreak());st.append(Paragraph('S17 原始文献与数据库来源',styles['h1']))
for r in refs:
 st.append(KeepTogether([P(f"{r['id']}  {r['year']}  {r['title']}"),Paragraph(f'<link href="{html.escape(r["url"],quote=True)}" color="#167e78">{html.escape(r["url"])}</link>',small),Spacer(1,9)]))
sup=O/'supplement/S17_candidate_selection_rationale.pdf'
def foot(c,d):
 c.saveState();c.setFont('CJK',8);c.setFillColor(colors.HexColor('#183448'));c.drawString(38,20,'补充表S17｜高效应基因选择依据');c.drawRightString(804,20,str(d.page));c.restoreState()
SimpleDocTemplate(str(sup),pagesize=(842,595),leftMargin=38,rightMargin=38,topMargin=30,bottomMargin=38,title='补充表S17 高效应基因选择依据').build(st,onFirstPage=foot,onLaterPages=foot)
(O/'provenance/narrative_changes.json').write_text(json.dumps(changes,ensure_ascii=False,indent=2),encoding='utf-8')
final=BeautifulSoup((O/'report_v11.html').read_text(encoding='utf-8'),'html.parser');assert '每批审阅10个' not in final.get_text();now=Counter(x.get_text() for x in final.find_all(['td','th']));assert all(now[x]>=n for x,n in original_cells.items())
for name,path in [('report',O/'silicosis_integrated_report_v11.pdf'),('S17',sup)]:
 d=fitz.open(path)
 for i,p in enumerate(d):
  assert '\ufffd' not in p.get_text();assert all(b[0]>=-1 and b[2]<=843 and b[1]>=0 and b[3]<=596 for b in p.get_text('blocks')),(name,i)
  p.get_pixmap(matrix=fitz.Matrix(1,1)).save(T/f'{name}_{i+1:02d}.png')
 print(name,'pages',len(d))
print('narrative changes',len(changes));shutil.copy2(__file__,O/'code/build_reports.py')
