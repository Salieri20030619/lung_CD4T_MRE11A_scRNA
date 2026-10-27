from pathlib import Path
import json,hashlib,zipfile,shutil,pandas as pd,fitz
W=Path('C:/Users/31937/lung_four_axes_20261017_v15');T=W/'tables';O=W/'deliverables';D=Path('C:/Users/31937/Downloads/线粒体焦亡_四线筛选_v15_20261017');D.mkdir(exist_ok=True,parents=True)
assert hashlib.sha256((W/'PLAN_FROZEN.json').read_bytes()).hexdigest()=='1852715bf9c6a5a38aae4e0dae903833a40b7562df26902cef06989894073ac5'
report=fitz.open(O/'four_axis_screen_report_v15.pdf');plan=fitz.open(O/'four_axis_frozen_plan_v15.pdf');assert len(report)==15 and len(plan)==2;text='\n'.join(p.get_text() for p in report)
for n in ['15,913','1,297','1,161','686,924','667,420','0.175','0.36111','0.044609','0.317495','34','DUSP5']:assert n in text,n
for doc in [report,plan]:
 for p in doc:assert '\ufffd' not in p.get_text();assert all(b[0]>=0 and b[1]>=0 and b[2]<=842 and b[3]<=595 for b in p.get_text('blocks'))
final=pd.read_csv(T/'05_intersection/final_candidate_decisions.tsv',sep='\t');assert len(final)==1 and final.gene.iloc[0]=='DUSP5' and final.frozen_four_axis_pass.all() and not final.recommend_as_settled_mechanistic_main.any()
nov=pd.read_csv(T/'06_novelty/final_novelty_assessment.tsv',sep='\t');assert nov.grade.iloc[0]=='B_context_specific'
val=json.loads((T/'provenance/numeric_validation.json').read_text());assert val['status']=='PASS' and val['BH_families']==48
for p in W.glob('*.py'):compile(p.read_text(encoding='utf-8'),str(p),'exec')
# Validate text table decoding, excluding gzip which pandas already read in computational QA.
for p in T.rglob('*.tsv'):
 s=p.read_text(encoding='utf-8');assert '\ufffd' not in s,p
readme='''四证据线联合筛选 v15｜2026-10-17

先读 four_axis_screen_report_v15.pdf（15页），再看 four_axis_frozen_plan_v15.pdf（2页）。
冻结计划 SHA256: 1852715bf9c6a5a38aae4e0dae903833a40b7562df26902cef06989894073ac5

结论分层：
1. 按冻结规则，四线交集为 DUSP5–Treg（1个）。宽松交集34个，包含上述1个。
2. 文献核查后，CDH1接触关系有适用基础；BMP1/OCLN/H2-Q1不能等价地解释为可靠上游。
3. 追加decontX诊断使DUSP5与焦亡程序关联减弱（合并P=0.175；组别调整P=0.36111）。该诊断在筛选之后追加，不回写冻结交集。
4. 当前没有可建议直接投入的成熟机制主线。DUSP5保留为筛选阳性但诊断不稳定的待核查对象。

tables/01_quality：全基因质量与样本值；02_DE：原完整差异；03_damage：完整损伤关联；04_macrophage：来源内完整联系；05_intersection：全交集与排除原因；06_novelty：创新性、受体资格、来源身份和模型诊断。
TSV为UTF-8无BOM，制表符分隔；NA不是零。gz文件可解压读取。校正count是模型估计，可为小数，并非新增测量的UMI。
原始BH在完整家族内保留。诊断来源与损伤分别处理，不将原候选重新在小家族中校正以升级证据。

完整代码、原始/中间矩阵、运行日志与环境在服务器：
/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/four_axis_screen_20261017_v15
本地工作目录：C:/Users/31937/lung_four_axes_20261017_v15
既有v13和v14未覆盖。本次交付为中文内部分析与独立复核材料，不是客户机制主线推荐，也未向外部接收方发送。
'''
(O/'README.txt').write_text(readme,encoding='utf-8');shutil.copy2(W/'PLAN_FROZEN.json',O/'PLAN_FROZEN.json')
qa={'status':'COMPLETE_READY_FOR_REVIEW','report_pages':15,'plan_pages':2,'plan_hash_unchanged':True,'frozen_four_axis_candidates':1,'settled_mechanistic_recommendations':0,'inclusive_exploratory_candidates':34,'BH_families_checked':48,'raw_count_exact_tests_checked':32,'historical_rosters_reproduced':48,'decontx_samples':6,'visual_review':'all 15 report pages and both plan pages inspected; tables and figures readable','numeric_roundtrip':'PASS','UTF8_tables':'PASS','not_verified':['independent-cohort validation','biochemical or causal mechanism','absolute endogenous source proof']};(O/'DELIVERY_VALIDATION.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
for f in O.glob('*.pdf'):shutil.copy2(f,D/f.name)
# Reviewer statistics only: no raw cell matrices, scripts or working logs in Downloads.
zipout=D/'v15_complete_report_and_statistics.zip'
with zipfile.ZipFile(zipout,'w',zipfile.ZIP_DEFLATED) as z:
 for f in O.iterdir():
  if f.is_file():z.write(f,f.name)
 for folder in ['01_quality','02_DE','03_damage','04_macrophage','05_intersection','06_novelty']:
  for f in (T/folder).glob('*'):
   if f.is_file() and ('cell_' not in f.name) and f.suffix in ['.tsv','.gz','.json']:z.write(f,'tables/'+folder+'/'+f.name)
with zipfile.ZipFile(zipout) as z:assert z.testzip() is None
# Small server authoring record; existing computation outputs remain where generated.
server=W/'v15_authoring_record.zip'
with zipfile.ZipFile(server,'w',zipfile.ZIP_DEFLATED) as z:
 for f in O.iterdir():
  if f.is_file():z.write(f,'report/final/'+f.name)
 for f in W.glob('*.py'):z.write(f,'code/final_authoring/'+f.name)
 for f in (W/'figures').iterdir():
  if f.is_file():z.write(f,'report/figures/'+f.name)
 z.write(W/'report_page_map.json','report/final/page_map.json')
result={'authoring_archive_sha256':hashlib.sha256(server.read_bytes()).hexdigest(),'authoring_archive_bytes':server.stat().st_size,'deliverables':[{'path':str(f),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in D.iterdir() if f.is_file()]};(W/'delivery_checksums.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(result,ensure_ascii=False));print(json.dumps(qa,ensure_ascii=False))