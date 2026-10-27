from pathlib import Path
import json,hashlib,shutil,zipfile,re,collections
import pandas as pd
import fitz
from bs4 import BeautifulSoup
W=Path('C:/Users/31937/lung_review_20261017_v14');O=W/'revision_v14';B=W/'baseline_v13';Q=O/'supplement_v14';D=Path('C:/Users/31937/Downloads/线粒体焦亡_评审修订_v14_20261017');D.mkdir(parents=True,exist_ok=True)
# Keep previous validation as historical evidence, replace current validation below.
old=O/'provenance/final_validation.json'
if old.exists():shutil.copy2(old,O/'provenance/v14/baseline_v13_final_validation.json')
for name in ['report_v13.html','silicosis_integrated_report_v13.pdf','MANIFEST_SHA256.tsv']:
 f=O/name;assert f.resolve().parent==O.resolve()
 if f.exists():f.unlink()
code=O/'code/v14';code.mkdir(exist_ok=True)
for f in W.glob('*.py'):shutil.copy2(f,code/f.name)
shutil.copy2('C:/codex_test/tmp/silicosis_v10/build_integrated.py',code/'render_engine_source_v10.py')
shutil.copy2(W/'review_original.txt',O/'provenance/v14/independent_review_original.txt')
oldhtml=BeautifulSoup((B/'report_v13.html').read_text(encoding='utf-8'),'html.parser');newhtml=BeautifulSoup((O/'report_v14.html').read_text(encoding='utf-8'),'html.parser')
assert [i['src'] for i in oldhtml.find_all('img')]==[i['src'] for i in newhtml.find_all('img')]
for i in oldhtml.find_all('img'):assert hashlib.sha256((B/i['src']).read_bytes()).digest()==hashlib.sha256((O/i['src']).read_bytes()).digest()
oldrows=collections.Counter(tuple(c.get_text() for c in tr.find_all(['td','th'])) for tr in oldhtml.find_all('tr'));newrows=collections.Counter(tuple(c.get_text() for c in tr.find_all(['td','th'])) for tr in newhtml.find_all('tr'));assert not (oldrows-newrows)
for f in (B/'supplement').iterdir():
 if f.is_file() and not f.name.endswith('.inspect.ndjson'):assert hashlib.sha256(f.read_bytes()).digest()==hashlib.sha256((O/'supplement'/f.name).read_bytes()).digest()
newtext=newhtml.get_text();assert '接收端调控明确后' not in newtext
pdf=fitz.open(O/'silicosis_integrated_report_v14.pdf');text='\n'.join(p.get_text() for p in pdf)
for v in ['0.27778','0.05556','0.02778','240','200','424','31','106','0.10']:assert v in text,v
assert len(pdf)==43
for docname,count in [('silicosis_integrated_report_v14.pdf',43),('review_response_v14.pdf',9)]:
 d=fitz.open(O/docname);assert len(d)==count
 for p in d:assert '\ufffd' not in p.get_text();assert all(b[0]>=0 and b[1]>=0 and b[2]<=842 and b[3]<=595 for b in p.get_text('blocks'))
for f in Q.iterdir():
 if f.suffix=='.tsv':
  s=f.read_text(encoding='utf-8');assert '\ufffd' not in s;pd.read_csv(f,sep='\t')
checks=pd.read_csv(Q/'S21_BH_recalculation_audit.tsv',sep='\t');assert len(checks)==106 and (checks.status=='PASS').all()
prox=pd.read_csv(Q/'S23_targeted_proximal_associations.tsv',sep='\t');assert len(prox)==240 and prox.P.notna().sum()==200 and (prox.P<.05).sum()==7 and (prox.BH<.05).sum()==0
validation={'version':'v14','main_pages':43,'response_pages':9,'original_figures_unchanged':18,'original_supplement_files_unchanged':len([f for f in (B/'supplement').iterdir() if not f.name.endswith('.inspect.ndjson')]),'original_table_rows_retained':sum(oldrows.values()),'BH_families_recomputed':106,'key_exact_tests_recomputed':4,'new_proximal_pairs':240,'evaluable_proximal_pairs':200,'new_nominal_proximal_pairs':7,'new_BH_significant_proximal_pairs':0,'visual_QA':'all pages rendered and inspected; focal risk and supplement index pagination repaired','status':'READY_FOR_INDEPENDENT_REVIEW','not_verified':['full upstream QC and annotation rerun','edgeR refit','SCENIC network rerun','NicheNet model rerun','pure CD4 endogenous GSDMD','protein and causal mechanism']}
(O/'provenance/final_validation.json').write_text(json.dumps(validation,ensure_ascii=False,indent=2),encoding='utf-8')
readme='''矽肺/线粒体焦亡项目 v14｜2026-10-17｜待再次独立评审

主报告：silicosis_integrated_report_v14.pdf（43页）
逐项回应：review_response_v14.pdf（9页，R01–R08）
原统计附件：supplement/S01–S20，原始值不变。
新增统计附件：supplement_v14/S21–S27，含全部原始P、BH、样本输入、阴性及NA。

本次完成：106完整BH家族复算、4种关键精确P复算、240项限定近端配对、完整受体/发送端组分、GSDMD原始UMI核查、regulon/匹配对照重叠与定向敏感性、全肺竞争来源描述。
未建立完整巨噬细胞—CD4损伤/焦亡因果链。GSDMD细胞来源、蛋白活化及功能连接仍需验证。旧v13保留在原目录，不覆盖。

TSV使用UTF-8无BOM，制表符分隔；NA为不可评估，不应替换成0。gz表可解压后读取。原S17另有xlsx阅读文件。
原始评分在S27仅为重现参照，不用缩小的敏感性家族重标原BH。
完整代码、来源、检查日志保留服务器：
/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/review_revision_20261017_v14

统计计算可用服务器batch_bench环境执行code/v14/analyze_v14.py及overlap_v14.py；后续交付层核验与排版使用本地validate_numeric.py、build_v14.py、response_v14.py和finalize_v14.py。脚本含冻结路径，复现需使用对应服务器及本地工作目录。完整报告生成不代表已获独立评审通过或获准外部发送。
'''
(O/'README.txt').write_text(readme,encoding='utf-8')
# Document original primary papers used for R06.
(O/'provenance/v14/literature_verification.json').write_text(json.dumps({'access_date':'2026-10-17','sources':[{'doi':'10.1038/s42003-024-06662-z','url':'https://www.nature.com/articles/s42003-024-06662-z','pubmed':'https://pubmed.ncbi.nlm.nih.gov/39122899/','use':'silica CD4 resident-memory and CXCR6 context; not novel here'},{'doi':'10.1084/jem.20181308','url':'https://rupress.org/jem/article/216/12/2748/132521/CXCR6-regulates-localization-of-tissue-resident','use':'airway epithelial CXCL16 in CD8 context; cannot directly validate CD4 mechanism'}]},indent=2),encoding='utf-8')
# Build full server record and reviewer-facing bundle separately.
manifest=[]
for f in sorted(O.rglob('*')):
 if f.is_file() and f.name!='MANIFEST_SHA256.tsv':manifest.append(dict(file=f.relative_to(O).as_posix(),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
pd.DataFrame(manifest).to_csv(O/'MANIFEST_SHA256.tsv',sep='\t',index=False,encoding='utf-8')
server_zip=W/'v14_server_record.zip'
with zipfile.ZipFile(server_zip,'w',zipfile.ZIP_DEFLATED) as z:
 for f in O.rglob('*'):
  if f.is_file():z.write(f,f.relative_to(O).as_posix())
for name in ['silicosis_integrated_report_v14.pdf','review_response_v14.pdf']:shutil.copy2(O/name,D/name)
reviewfiles=[O/'silicosis_integrated_report_v14.pdf',O/'review_response_v14.pdf',O/'review_response.tsv',O/'README.txt',O/'provenance/final_validation.json']+list((O/'supplement').glob('*'))+list(Q.glob('*'))
with zipfile.ZipFile(D/'v14_review_materials.zip','w',zipfile.ZIP_DEFLATED) as z:
 for f in reviewfiles:
  if f.is_file():z.write(f,f.relative_to(O).as_posix())
for zpath in [server_zip,D/'v14_review_materials.zip']:
 with zipfile.ZipFile(zpath) as z:assert z.testzip() is None
for name in ['silicosis_integrated_report_v14.pdf','review_response_v14.pdf']:assert hashlib.sha256((O/name).read_bytes()).digest()==hashlib.sha256((D/name).read_bytes()).digest()
result={'server_archive_sha256':hashlib.sha256(server_zip.read_bytes()).hexdigest(),'server_archive_bytes':server_zip.stat().st_size,'delivery':[{ 'file':str(f),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in D.iterdir() if f.is_file()]};(W/'delivery_validation.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(result,ensure_ascii=False));print(json.dumps(validation,ensure_ascii=False))