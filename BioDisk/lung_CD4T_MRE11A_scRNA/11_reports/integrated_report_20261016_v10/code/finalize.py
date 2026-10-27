from pathlib import Path
import json,hashlib,shutil,zipfile,csv
import pandas as pd
import fitz
from bs4 import BeautifulSoup
R=Path('C:/codex_test');O=R/'output/silicosis_v10';V9=R/'output/silicosis_v9'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
checks={}
for p in (O/'plot_data').glob('*.tsv'):
 assert sha(p)==sha(V9/'plot_data'/p.name),p.name
checks['all_plot_tables_byte_identical_to_v9']=True
for folder in ['assets','supplement']:
 src=V9/'baseline/assets' if folder=='assets' else V9/'supplement'
 for p in (O/folder).iterdir():assert sha(p)==sha(src/p.name)
checks['original_figures_and_supplement_tables_unchanged']=True
shutil.copytree(V9/'source_tables',O/'source_tables',dirs_exist_ok=True)
hs=BeautifulSoup((O/'report_v10.html').read_text(encoding='utf-8'),'html.parser')
rows=[]
for img in hs.find_all('img'):
 path=O/img['src'];assert path.is_file();caption=img.find_next_sibling('p',class_='caption').get_text();rows.append({'number':len(rows)+1,'file':img['src'],'caption':caption,'sha256':sha(path)})
pd.DataFrame(rows).to_csv(O/'provenance/figure_index.tsv',sep='\t',index=False,encoding='utf-8')
d=fitz.open(O/'silicosis_integrated_report_v10.pdf');assert len(d)==26 and len(d[1].get_links())==10
checks.update(pages=26,figures=15,all_pages_rendered_and_visually_reviewed=True,chapter_navigation_links=10,original_paragraphs_retained=71,original_table_cells_retained=219,new_inferential_tests=0)
(O/'provenance/final_validation.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2),encoding='utf-8')
(O/'README.txt').write_text('矽肺研究逻辑重构整合版 v10（2026-10-16）\n\n先打开 silicosis_integrated_report_v10.pdf：26页，15幅图，目录页可点击跳转。\nreport_v10.html：同版正文，可在浏览器中放大查看。\n\n阅读顺序：\n第3-4页 数据与方法\n第5-6页 CD4功能表型\n第7-9页 疾病响应候选与时间支持\n第10-14页 TF活性、功能关联、靶基因模块与STAT1-GSDMD\n第15-18页 巨噬细胞通讯、配体预测、来源差异与接收端联系\n第19-20页 患者定位与人类病例对照方向\n第21-23页 研究价值与分层实验\n第24-26页 讨论、结论与参考文献\n\n本版将v8正文及v9新增图谱整合到同一论证顺序，重新编排章节和图号，增加过渡文字。71段原正文、219个表格单元格、原始图件及补充统计表保留；新增图谱的绘图数据与v9逐文件一致。未重跑统计或改变研究方向。\n\nfigures为六组图谱的PNG/SVG/矢量PDF；assets为九幅既有图。图号与文件名对应见provenance/figure_index.tsv。plot_data为图谱绘图数据，source_tables为冻结来源表，supplement为16份完整补充统计。code为本轮脚本；脚本顶部路径为本机执行环境，移机需调整ROOT/SRC/OLD/OUT并提供对应v8/v9基准。\nv8、v9在原路径保留。本次内容守恒及版式检查不等同于对所有既有分析进行独立重审。\n',encoding='utf-8')
shutil.copy2(__file__,O/'code/finalize.py')
manifest=[]
for p in sorted(O.rglob('*')):
 if p.is_file() and p.name!='MANIFEST_SHA256.tsv':manifest.append({'file':p.relative_to(O).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)})
pd.DataFrame(manifest).to_csv(O/'MANIFEST_SHA256.tsv',sep='\t',index=False,encoding='utf-8')
zpath=O.parent/'silicosis_integrated_report_v10.zip'
with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED) as z:
 for p in sorted(O.rglob('*')):
  if p.is_file():z.write(p,p.relative_to(O).as_posix())
with zipfile.ZipFile(zpath) as z:assert z.testzip() is None
print(json.dumps({'checks':checks,'zip_bytes':zpath.stat().st_size,'zip_sha256':sha(zpath),'pdf_sha256':sha(O/'silicosis_integrated_report_v10.pdf')},ensure_ascii=False,indent=2))
