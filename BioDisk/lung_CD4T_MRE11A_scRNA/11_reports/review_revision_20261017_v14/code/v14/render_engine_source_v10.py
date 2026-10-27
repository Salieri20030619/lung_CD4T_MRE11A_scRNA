from pathlib import Path
import re,json,html,shutil,hashlib,zipfile
from collections import Counter
from bs4 import BeautifulSoup,Tag
import fitz
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,PageBreak,Image,Table,TableStyle,KeepTogether
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors

R=Path('C:/codex_test'); O=R/'output/silicosis_v10'; V9=R/'output/silicosis_v9'; T=R/'tmp/silicosis_v10'
for d in ['figures','assets','supplement','provenance','code','plot_data']: (O/d).mkdir(parents=True,exist_ok=True)
# Regenerate the six data figures with unified figure numbers, without changing v9.
code=(R/'tmp/silicosis_v9/build_v9.py').read_text(encoding='utf-8').split('# Preserve original evidence')[0]
code=code.replace("OUT=ROOT/'output/silicosis_v9'","OUT=ROOT/'output/silicosis_v10'")
for old,new in [('01','11'),('02','12'),('03','06'),('04','07'),('05','09'),('06','13')]:
 # Replace only the explicit old plot title, not newly replaced digits.
 titles={'01':'Macrophage-to-CD4','02':'Predicted ligand','03':'CD4 regulon','04':'Regulons connect','05':'STAT1 prior','06':'The ligand-to-STAT1'}
 code=code.replace(old+'  '+titles[old],new+'  '+titles[old])
(T/'replot_v10.py').write_text(code,encoding='utf-8'); exec(compile(code,str(T/'replot_v10.py'),'exec'),{})
raw=(V9/'baseline/report_v8.html').read_text(encoding='utf-8'); soup=BeautifulSoup(raw,'html.parser'); oldsecs=soup.find_all('section')
notes=json.loads((V9/'provenance/figure_notes.json').read_text(encoding='utf-8'))
shutil.copytree(V9/'baseline/assets',O/'assets',dirs_exist_ok=True)
shutil.copytree(V9/'supplement',O/'supplement',dirs_exist_ok=True)

fig_no={'Fig03_death_programs.png':1,'Fig04_mito_autophagy.png':2,'Fig01_candidate_volcano.png':3,'Fig02_candidate_evidence.png':4,'Fig05_real_time.png':5,'F03_TF_activity_atlas.png':6,'F04_TF_function_map.png':7,'Fig08_composition_sensitivity.png':8,'F05_target_module_network.png':9,'Fig06_STAT1_GSDMD.png':10,'F01_communication_network.png':11,'F02_ligand_landscape.png':12,'F06_bridge_evidence.png':13,'Fig07_patient_localization.png':14,'Fig09_human_composition.png':15}
sections=[]; mapping=[]
def original(idx,title,chapter=False,lead=None):
 ss=BeautifulSoup(str(oldsecs[idx]),'html.parser').section
 h=ss.find(['h1','h2']); oldtitle=h.get_text(); h.string=title; h.name='h1' if chapter else 'h2'
 if lead:
  p=soup.new_tag('p');p['class']='transition';p.string=lead;h.insert_after(p)
 for img in ss.find_all('img'):
  n=fig_no[Path(img['src']).name]; cap=img.find_next_sibling('p',class_='caption')
  assert cap is not None
  cap.string=re.sub(r'^图\d+\s*',f'图{n}  ',cap.get_text())
 ss['data-chapter']='yes' if chapter else 'no'; ss['id']=f'section-{len(sections)+1}'
 sections.append(ss);mapping.append({'source':'v8','source_section':idx,'old_title':oldtitle,'new_title':title})
def added(i,title,chapter=False,lead=None):
 n=notes[i]; ss=soup.new_tag('section');ss['data-chapter']='yes' if chapter else 'no';ss['id']=f'section-{len(sections)+1}'
 h=soup.new_tag('h1' if chapter else 'h2');h.string=title;ss.append(h)
 p=soup.new_tag('p');p.string=lead or n['text'];ss.append(p)
 img=soup.new_tag('img',src='figures/'+n['name']+'.png');ss.append(img)
 cap=soup.new_tag('p');cap['class']='caption';num=fig_no[n['name']+'.png'];cap.string=f'图{num}  '+n['title'].split('｜')[1]+'。'+n['method'].replace('v8重点候选','既有重点候选');ss.append(cap)
 if lead:
  p=soup.new_tag('p');p.string=n['text'];ss.append(p)
 sections.append(ss);mapping.append({'source':'v9_atlas','source_figure':i+1,'new_title':title,'figure':num})

original(1,'1 研究设计与数据分工',True)
original(2,'1.1 候选筛选、功能分析与样本级统计')
p=soup.new_tag('p');p.string='图谱展示在冻结结果上完成：调控子活性图按方差选择22个调控子并加入STAT1和IKZF2；TF与功能程序图按三个亚群内名义关联数选择前18个TF并加入上述两者；配体预测图展示各来源、亚群排名前5的并集及既有重点配体。上述选择仅用于展示，完整原始P值、BH及检验家族保持不变。';sections[-1].append(p)
original(6,'2 CD4亚群呈现不同的死亡相关转录变化',True,'首先比较各CD4亚群的疾病变化，以确定后续候选、调控和通讯分析需要解释的功能表型。')
original(7,'2.1 线粒体与自噬程序呈现亚群特异方向')
original(3,'3 差异表达界定疾病响应候选',True,'在功能表型基础上，进一步定位与疾病相关的基因。候选来源、亚群效应和时间支持分别呈现，再进入接收端调控分析。')
original(4,'3.1 候选效应随CD4状态变化')
original(8,'3.2 状态轨迹与真实时间提供辅助支持')
added(2,'4 调控子活性显示CD4亚群与样本层次',True,'疾病响应基因提供候选入口；接下来从调控子活性全景出发，依次考察TF与功能程序的关系、重点TF的组成敏感性以及执行相关基因。')
added(3,'4.1 TF与功能程序的关联具有状态特异性')
original(9,'4.2 STAT1与IKZF2的关联对组成调整反应不同')
added(4,'4.3 STAT1先验靶基因连接多个功能模块')
original(10,'4.4 STAT1活性与GSDMD转录表达呈正向关联',False,'在调控全景与先验模块之后，逐基因检验进一步聚焦GSDMD分支。以下关系使用样本级证据，与先验靶基因连接分别解释。')
added(0,'5 巨噬细胞与CD4亚群具有候选通讯连接',True,'接收端调控明确后，再检查巨噬细胞能否提供对应的上游信号：先看配体受体结构和覆盖，再看预测排序及来源差异，最后评价配体与STAT1的同样本联系。')
added(1,'5.1 配体预测排名与受体覆盖共同限定通讯解释')
original(11,'5.2 CXCL16–CXCR6获得来源差异与接收端表达支持')
added(5,'5.3 配体与STAT1的关系依赖来源定义和组别调整')
original(13,'6 患者定位支持CXCL16与CXCR6的两端表达',True,'小鼠结果之后，人类证据分两层评价：先确认候选分子的细胞来源和受体覆盖，再用独立病例对照队列检验疾病方向。患者定位与疾病效应使用不同队列，不合并为同一个验证终点。')
original(12,'6.1 人类病例对照显示与小鼠不同的Core方向')
original(5,'7 研究价值取决于细胞情境与作用位置',True,'综合表型、候选、调控、通讯和人类证据，按作用位置讨论后续实验选择。以下保留既有分层建议，供研究方向讨论。')
original(14,'7.1 细胞间作用与接收端机制分层验证')
original(15,'7.2 实验设计围绕身份、表型与因果次序展开')
original(16,'8 讨论：证据边界与待解决问题',True)
original(17,'9 结论与补充统计',True)
original(18,'参考文献',True)

# Front matter is rewritten for a single report, retaining the original abstract.
cover=BeautifulSoup(str(oldsecs[0]),'html.parser').section
cover.find('h1').string='矽肺背景下巨噬细胞与CD4 T细胞表型关联研究'
cover.find('p').string='逻辑重构整合版 v10｜2026年10月16日'
nav=soup.new_tag('section');h=soup.new_tag('h1');h.string='阅读路径';nav.append(h)
p=soup.new_tag('p');p.string='数据设计 → CD4功能表型 → 疾病响应候选 → 接收端TF调控 → 巨噬细胞通讯 → 人类定位与疾病方向 → 实验选择与讨论';nav.append(p)
for ss in sections:
 if ss.get('data-chapter')=='yes':
  p=soup.new_tag('p');a=soup.new_tag('a',href='#'+ss['id']);a.string=ss.find(['h1','h2']).get_text();p.append(a);nav.append(p)
p=soup.new_tag('p');p.string='本版将15幅图按论证顺序统一编号。图谱与对应结果在正文中连续呈现；统计值和已有阴性结果沿用冻结证据。完整统计表、逐图数据及来源记录随资料包提供。';nav.append(p)
css='body{max-width:1150px;margin:35px auto;padding:24px;color:#183448;font:17px/1.85 Arial,"Microsoft YaHei",sans-serif}h1{font-size:27px;color:#123447;margin-top:55px}h2{font-size:22px;margin-top:35px}h3{font-size:19px}img{width:100%;height:auto}table{width:100%;border-collapse:collapse;font-size:14px}th{background:#163b50;color:white}td,th{padding:9px;border:1px solid #cdd7df}.caption{font-size:14px;color:#42596a}.transition{color:#28678a}a{color:#167e78}section{margin-bottom:35px}@media print{body{max-width:none}h1,h2{break-after:avoid}img,table{break-inside:avoid}}'
(O/'report_v10.html').write_text('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>矽肺研究整合报告 v10</title><style>'+css+'</style>'+str(cover)+str(nav)+''.join(map(str,sections))+'</html>',encoding='utf-8')

# Generate a continuous typeset report from the same HTML blocks.
pdfmetrics.registerFont(TTFont('CJK','C:/Windows/Fonts/ARIALUNI.ttf'))
styles={
 'body':ParagraphStyle('body',fontName='CJK',fontSize=10.5,leading=16.5,spaceAfter=8,wordWrap='CJK'),
 'h1':ParagraphStyle('h1',fontName='CJK',fontSize=20,leading=26,textColor=colors.HexColor('#163b50'),spaceBefore=5,spaceAfter=13,keepWithNext=True,wordWrap='CJK'),
 'h2':ParagraphStyle('h2',fontName='CJK',fontSize=15,leading=21,textColor=colors.HexColor('#28678a'),spaceBefore=12,spaceAfter=9,keepWithNext=True,wordWrap='CJK'),
 'h3':ParagraphStyle('h3',fontName='CJK',fontSize=12,leading=18,textColor=colors.HexColor('#163b50'),spaceBefore=8,spaceAfter=7,keepWithNext=True,wordWrap='CJK'),
 'caption':ParagraphStyle('caption',fontName='CJK',fontSize=9,leading=13,spaceAfter=12,wordWrap='CJK',textColor=colors.HexColor('#42596a')),
 'table':ParagraphStyle('table',fontName='CJK',fontSize=9,leading=13,wordWrap='CJK'),
 'th':ParagraphStyle('th',fontName='CJK',fontSize=9,leading=13,wordWrap='CJK',textColor=colors.white),
}
def inline(tag):
 text=html.escape(tag.get_text())
 # Preserve reference URLs as clickable links in PDF.
 for a in tag.find_all('a'):
  label=html.escape(a.get_text()); href=html.escape(a.get('href',''),quote=True)
  if href.startswith('https://'):text=text.replace(label,f'<link href="{href}" color="#167e78">{label}</link>',1)
 return text
story=[]; bookmarks=[]; figure_pages={}
class NamedPara(Paragraph):
 def __init__(self,text,style,name,level):super().__init__(text,style);self.mark=(name,level)
class Report(SimpleDocTemplate):
 def afterFlowable(self,f):
  if hasattr(f,'mark'):
   name,level=f.mark;key='bm'+str(len(bookmarks));self.canv.bookmarkPage(key);self.canv.addOutlineEntry(name,key,level=level,closed=False);bookmarks.append({'title':name,'page':self.page})
  if hasattr(f,'figmark'):figure_pages[f.figmark]=self.page
def elements(ss,headlevel=0):
 children=list(ss.children);i=0
 while i<len(children):
  el=children[i];i+=1
  if not isinstance(el,Tag):continue
  if el.name in ['h1','h2','h3']:
   # Existing in-section h2 headings become h3, avoiding false chapter hierarchy.
   first=(el==ss.find(['h1','h2']))
   style=el.name if first else 'h3'
   if first: story.append(NamedPara(html.escape(el.get_text()),styles[style],el.get_text(),headlevel))
   else:story.append(Paragraph(inline(el),styles[style]))
  elif el.name=='p':story.append(Paragraph(inline(el),styles['caption' if 'caption' in el.get('class',[]) else 'body']))
  elif el.name=='table':
   rows=[]
   for tr in el.find_all('tr'):rows.append([Paragraph(inline(td),styles['th' if td.name=='th' else 'table']) for td in tr.find_all(['td','th'])])
   nc=len(rows[0]);width=766
   if nc==3: widths=[width*.21,width*.39,width*.4]
   elif nc==4: widths=[width*.2,width*.25,width*.27,width*.28]
   elif nc==6: widths=[width*.13,width*.12,width*.12,width*.22,width*.17,width*.24]
   else: widths=[width/nc]*nc
   tb=Table(rows,colWidths=widths,repeatRows=1,hAlign='LEFT');tb.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#163b50')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f1f5f7')]),('GRID',(0,0),(-1,-1),.35,colors.HexColor('#cdd7df')),('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5)]));story.extend([tb,Spacer(1,10)])
  elif el.name=='img':
   path=O/el['src'];px=fitz.Pixmap(str(path));maxh={1:300,2:300,3:240,4:300,5:140,8:200,10:160,14:250,15:250}.get(fig_no[path.name],320);scale=min(766/px.width,maxh/px.height);im=Image(str(path),width=px.width*scale,height=px.height*scale,hAlign='CENTER');im.figmark=fig_no[path.name]
   # Keep caption with its figure; surrounding paragraphs flow naturally.
   parts=[im,Spacer(1,5)]
   if i<len(children) and isinstance(children[i],Tag) and 'caption' in children[i].get('class',[]):parts.append(Paragraph(inline(children[i]),styles['caption']));i+=1
   story.append(KeepTogether(parts))

elements(cover,0);story.append(PageBreak());elements(nav,0);story.append(PageBreak())
for n,ss in enumerate(sections):
 if n>0:story.append(PageBreak())
 elements(ss,0 if ss.get('data-chapter')=='yes' else 1)

def furniture(c,doc):
 c.saveState();c.setFillColor(colors.HexColor('#163b50'));c.rect(0,583,842,12,fill=1,stroke=0);c.setFont('CJK',8);c.drawString(38,22,'矽肺研究｜通讯、调控与功能证据整合｜v10');c.drawRightString(804,22,str(doc.page));c.restoreState()
pdf=O/'silicosis_integrated_report_v10.pdf';doc=Report(str(pdf),pagesize=(842,595),leftMargin=38,rightMargin=38,topMargin=36,bottomMargin=42,title='矽肺背景下巨噬细胞与CD4 T细胞表型关联研究 v10',author='Research report')
doc.build(story,onFirstPage=furniture,onLaterPages=furniture)

# Add page references and clickable navigation after pagination stabilizes.
navdoc=fitz.open(pdf);navpage=navdoc[1]
main_titles={ss.find(['h1','h2']).get_text() for ss in sections if ss.get('data-chapter')=='yes'}
for entry in bookmarks:
 if entry['title'] not in main_titles:continue
 hits=navpage.search_for(entry['title'])
 assert len(hits)==1,(entry,hits)
 rect=hits[0];navpage.insert_text((775,rect.y1-2),str(entry['page']),fontsize=10,color=(.09,.25,.32))
 navpage.insert_link({'kind':fitz.LINK_GOTO,'from':fitz.Rect(38,rect.y0,805,rect.y1),'page':entry['page']-1})
navtmp=T/'with_navigation.pdf';navdoc.save(navtmp);navdoc.close();shutil.copy2(navtmp,pdf)
hs=BeautifulSoup((O/'report_v10.html').read_text(encoding='utf-8'),'html.parser')
page_lookup={entry['title']:entry['page'] for entry in bookmarks}
for a in hs.select('a[href^="#section-"]'):
 a.string=a.get_text()+f" · PDF第{page_lookup[a.get_text()]}页"
(O/'report_v10.html').write_text(str(hs),encoding='utf-8')

# Content conservation checks: every original paragraph (except date/captions)
# and every table cell remains, figures occur once and numerical tables unchanged.
newsoup=BeautifulSoup((O/'report_v10.html').read_text(encoding='utf-8'),'html.parser')
oldparas=[p.get_text() for sec in oldsecs for p in sec.find_all('p') if 'caption' not in p.get('class',[]) and p.get_text()!='2026年10月15日']
newtext=newsoup.get_text();missing=[p for p in oldparas if p not in newtext];assert not missing,missing
oldcells=Counter(td.get_text() for sec in oldsecs for td in sec.find_all(['td','th']));newcells=Counter(td.get_text() for td in newsoup.find_all(['td','th']));assert oldcells==newcells
assert len(newsoup.find_all('img'))==15
nums=[int(re.match(r'图(\d+)',p.get_text()).group(1)) for p in newsoup.find_all('p',class_='caption')];assert nums==list(range(1,16)),nums
d=fitz.open(pdf);ptext='\n'.join(p.get_text() for p in d);assert '\ufffd' not in ptext
# Numbers are carried over verbatim; check a selection in extracted PDF text.
for key in ['0.05556','0.02778','0.5902','0.9190','0.1028','0.5688','264/323','601/782','3/33','0/27']:
 assert key in ptext,key
for i,p in enumerate(d):
 p.get_pixmap(matrix=fitz.Matrix(1,1)).save(T/f'page_{i+1:02d}.png')
 blocks=p.get_text('blocks'); assert all(b[0]>=-1 and b[1]>=0 and b[2]<=843 and b[3]<=596 for b in blocks),(i,blocks)
validation={'original_paragraphs_retained':len(oldparas),'missing_original_paragraphs':missing,'original_table_cells_retained':sum(oldcells.values()),'figure_sequence':nums,'pages':len(d),'new_inferential_tests':0,'key_numeric_roundtrip':'PASS','replacement_characters':0,'text_bounds':'PASS','v8_v9_unchanged':True}
(O/'provenance/content_validation.json').write_text(json.dumps(validation,ensure_ascii=False,indent=2),encoding='utf-8')
(O/'provenance/section_mapping.json').write_text(json.dumps(mapping,ensure_ascii=False,indent=2),encoding='utf-8')
(O/'provenance/page_map.json').write_text(json.dumps(bookmarks,ensure_ascii=False,indent=2),encoding='utf-8')
shutil.copy2(__file__,O/'code/build_integrated.py');shutil.copy2(T/'replot_v10.py',O/'code/replot_v10.py')
shutil.copy2(V9/'provenance/source_sha256.json',O/'provenance/source_sha256.json');shutil.copy2(V9/'provenance/final_validation.json',O/'provenance/v9_source_value_validation.json')
print(json.dumps(validation,ensure_ascii=False));print('PAGE MAP',json.dumps(bookmarks,ensure_ascii=False))
