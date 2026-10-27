from pathlib import Path
import sys,json,hashlib,zipfile,csv,re,datetime
O=Path(sys.argv[1]);q=json.loads((O/'provenance/report_validation.json').read_text(encoding='utf-8'));assert q['status']=='PASS'
for name in ['literature_search_v3.py','search_supplement_v3.py','direction_v3.py','depth_v3.py','report_v3.py']:
 assert json.loads((O/'provenance'/(name+'.exit.json')).read_text(encoding='utf-8'))['exit_code']==0
s=(O/'research_memo_zh.html').read_text(encoding='utf-8')
for p in re.findall(r'<img src="([^"]+)"',s):assert (O/p).is_file()
(O/'provenance/visual_QA.json').write_text(json.dumps({'status':'PASS','figures_inspected':3,'image_links_verified':3,'UTF8_report_readback':'PASS'},indent=2),encoding='utf-8')
status={'status':'COMPLETE_WITH_EXPLICIT_UNRESOLVED_BIOLOGICAL_QUESTIONS','literature_review':'33 reproducible queries;18 curated studies;gene-family identity and model limits documented','direction_diagnostics':'completed within available data','partly_explained':['TTC39C D7 estimand sensitivity','cohort sampling 20weeks versus D7/D56'],'unresolved':['MDFIC external direction','TTC39C D56 direction','AHNAK human NaiveCD4 direction','biological cause of cross-stage differences'],'candidate_priority':['ST6GALNAC3 first exploratory validation','AHNAK parallel mechanistic-context candidate','TTC39C matched-stage validation','MDFIC deprioritize currently'],'v2_unchanged_files':125,'new_cell_count_gate':False,'external_cell_level_inferential_tests':0,'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(O/'RESEARCH_COMPLETE.json').write_text(json.dumps(status,indent=2),encoding='utf-8')
(O/'README.txt').write_text('''候选创新性与方向矛盾诊断 v3（内部研究）

入口：research_memo_zh.html
literature/curated_evidence.tsv：18项关键原始研究与预印本，区分功能、表达关联及模型。
literature/search_log.tsv + supplemental_search_log.tsv：33条可复核检索式。
literature/novelty_priority_assessment.tsv：创新空间与综合优先级。
direction/findings_zh_en.tsv：D01-D04稳定编号、中文解释及可转发英文反馈。
direction/cohort_time_evidence.tsv：发现队列20周与外部D7/D56的原文依据。
direction/：归一化估计量、检出率/强度分解、身份/深度标准化及人类删一供者记录。
figures/：3张诊断图。

完成范围：文献检索与现有数据内诊断；不是所有生物学矛盾已被消除。
TTC39C D7部分由估计量解释；D56、MDFIC及人类AHNAK的方向问题仍需独立验证。
原v2结果、P/BH和候选规则保持不变。未按细胞数删除样本，未做外部单样本细胞级P检验。

服务器上的 literature/cache 保存核读缓存；为避免分发原文，缓存与含长原文摘录的日志不纳入此包。
代码和完成状态保留；原始数据、文献原文和部分查询结果仍依赖外部源。
This is an internal research memo, not a final customer report or external communication.
''',encoding='utf-8')
files=[]
for f in sorted(O.rglob('*')):
 if not f.is_file():continue
 rel=f.relative_to(O)
 if 'cache' in rel.parts:continue
 if f.name.endswith('_public_metadata.txt'):continue
 if f.name in ['OUTPUT_MANIFEST_SHA256.tsv','package_v3.py.log','package_v3.py.exit.json']:continue
 if f.suffix=='.log' and f.name not in ['direction_v3.py.log','depth_v3.py.log','report_v3.py.log']:continue
 if f.suffix in ['.json','.py','.txt','.tsv','.html']:
  f.read_text(encoding='utf-8');assert not f.read_bytes().startswith(b'\xef\xbb\xbf')
 files.append({'path':str(rel),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
with (O/'OUTPUT_MANIFEST_SHA256.tsv').open('w',encoding='utf-8',newline='') as h:
 w=csv.DictWriter(h,fieldnames=['path','bytes','sha256'],delimiter='\t');w.writeheader();w.writerows(files)
zpath=O.parent/(O.name+'.zip')
with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED) as z:
 for r in files:z.write(O/r['path'],r['path'])
 z.write(O/'OUTPUT_MANIFEST_SHA256.tsv','OUTPUT_MANIFEST_SHA256.tsv')
with zipfile.ZipFile(zpath) as z:assert z.testzip() is None
print(json.dumps({'zip':str(zpath),'bytes':zpath.stat().st_size,'sha256':hashlib.sha256(zpath.read_bytes()).hexdigest(),'files':len(files)}))