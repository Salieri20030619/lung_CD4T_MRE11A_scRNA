from pathlib import Path
import json,hashlib,shutil,zipfile,re
from collections import Counter
import pandas as pd
from bs4 import BeautifulSoup
import fitz
R=Path('C:/codex_test');O=R/'output/silicosis_v13';T=R/'tmp/silicosis_v13';OLD=R/'output/silicosis_v12'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
old=BeautifulSoup((OLD/'report_v12.html').read_text(encoding='utf-8'),'html.parser');new=BeautifulSoup((O/'report_v13.html').read_text(encoding='utf-8'),'html.parser')
oldr=old.find_all('section')[5:];newr=new.find_all('section')[11:];assert len(oldr)==len(newr)
allowed='评分单位和样本汇总方式与1.2节一致'
for a,b in zip(oldr,newr):
 assert a.find(['h1','h2']).get_text()==b.find(['h1','h2']).get_text()
 if a.find(['h1','h2']).get_text().startswith('9 '):
  rows=b.find('table').find_all('tr');assert rows[-1].find('td').get_text()=='S20'
  copy=BeautifulSoup(str(b),'html.parser').section;copy.find('table').find_all('tr')[-1].decompose()
  assert a.get_text()==copy.get_text()
 else:assert a.get_text().replace(allowed,'评分单位和样本汇总方式与1.3节一致')==b.get_text()
oldimgs=[img['src'] for a in oldr for img in a.find_all('img')];newimgs=[img['src'] for a in newr for img in a.find_all('img')];assert oldimgs==newimgs and len(oldimgs)==18
for src in oldimgs:assert sha(OLD/src)==sha(O/src)
for p in (OLD/'supplement').iterdir():assert sha(p)==sha(O/'supplement'/p.name)
sources=pd.read_csv(O/'provenance/S20_method_source_manifest.tsv',sep='\t')
for r in sources.itertuples():assert sha(O/r.package_path)==r.sha256
pdf=fitz.open(O/'silicosis_integrated_report_v13.pdf');assert len(pdf)==37
txt='\n'.join(page.get_text() for page in pdf)
for term in ['Scrublet','Harmony','filterByExpr','glmQLFit','GRNBoost2','DropOverlap','DoRothEA','K15_R1','64种','720种','36种','GSE174725','GSE296682','7204','S20']:assert term in txt,term
assert '\ufffd' not in txt
maps=json.loads((O/'provenance/page_map.json').read_text(encoding='utf-8'));methods=[m for m in maps if re.match(r'^1\.[1-8] ',m['title'])];assert len(methods)==8
shutil.copytree(OLD/'provenance',O/'provenance/baseline_v12',dirs_exist_ok=True)
pd.DataFrame([{'figure':i+1,'file':src,'sha256':sha(O/src)} for i,src in enumerate(newimgs)]).to_csv(O/'provenance/figure_index.tsv',sep='\t',index=False,encoding='utf-8')
checks={'pages':37,'methods_pages':[4,5,6,7,8,9,10,11],'methods_subsections':8,'result_sections_preserved':len(oldr),'result_changes':'Only score-method cross-reference 1.2 -> 1.3 and S20 index row added','figures_preserved':18,'prior_supplement_files_preserved':len(list((OLD/'supplement').iterdir())),'verified_method_source_files':len(sources),'new_inferential_tests':0,'replacement_characters':0,'visual_review':'All eight Methods pages inspected; unchanged result layouts preserved.'}
(O/'provenance/final_validation.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2),encoding='utf-8')
(O/'README.txt').write_text('矽肺报告 v13：完整方法学修订版\n\n主报告 silicosis_integrated_report_v13.pdf 共37页。第4-11页为完整方法学：\n1.1 预处理、QC、细胞注释\n1.2 edgeR伪批量差异与候选筛选\n1.3 功能评分、完整基因集、匹配对照\n1.4 SCENIC、GRNBoost2、cisTarget、AUCell及TF先验\n1.5 DPT、样本级斜率、真实时间edgeR与交集\n1.6 配体受体资格、NicheNet与发送端关联\n1.7 人类注释复核、供者比较、患者定位及组成分解\n1.8 统计家族、置换、敏感性及绘图规则\n\n图1-18、全部既有统计表及结果正文保持原值。仅更新评分方法交叉引用和补充表索引。S20提供方法-参数-结果-来源对照，另附完整预设模块及广义自噬成员、匹配对照。方法来源脚本及SHA256在provenance/method_sources和S20_method_source_manifest.tsv。\n本轮未重跑分析。历史阶段方法按最终结果实际采用的规则说明，旧门槛与当前汇总分开记录。\n面向用户的完整报告复制到Windows下载目录；源数据、中间结果、代码与审计证据不自动复制到下载目录。\n',encoding='utf-8')
shutil.copy2(__file__,O/'code/finalize_v13.py')
rows=[]
for p in sorted(O.rglob('*')):
 if p.is_file() and p.name!='MANIFEST_SHA256.tsv' and not p.name.endswith('.inspect.ndjson'):rows.append({'file':p.relative_to(O).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)})
pd.DataFrame(rows).to_csv(O/'MANIFEST_SHA256.tsv',sep='\t',index=False,encoding='utf-8')
zp=O.parent/'silicosis_methods_complete_v13.zip'
with zipfile.ZipFile(zp,'w',zipfile.ZIP_DEFLATED) as z:
 for row in rows:z.write(O/row['file'],row['file'])
 z.write(O/'MANIFEST_SHA256.tsv','MANIFEST_SHA256.tsv')
with zipfile.ZipFile(zp) as z:assert z.testzip() is None
print(json.dumps(checks,ensure_ascii=False));print('ZIP_SHA256',sha(zp));print('PDF_SHA256',sha(O/'silicosis_integrated_report_v13.pdf'))
