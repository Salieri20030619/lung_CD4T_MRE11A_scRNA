"""Build a non-destructive, code-only project delivery. Never execute analysis code."""
from pathlib import Path
import os, re, csv, json, hashlib, shutil, ast, subprocess, collections, datetime
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
O=P/'11_reports/review_revision_20261006_v6'
C=O/'package'
for d in ['code','documentation','provenance','report']:(C/d).mkdir(parents=True,exist_ok=True)
CODE={'.py','.r','.sh','.ps1'}
groups={
 '01_qc_annotation':('原始数据质控与注释','方法1.2–1.3；结果2.1；图1–2'),
 '02_cd4_pseudobulk':('CD4与供者伪批量差异','方法1.4、1.8；结果2.1–2.2'),
 '03_function_programs':('功能GSEA/GSVA/UCell及12程序','方法1.8–1.9；结果2.10–2.13；图7–10、S2–S4；表S2–S8、S13'),
 '04_scenic':('SCENIC及调控子评分','方法1.6；结果2.5–2.8、2.12；图4–5、10、S1；表S5–S6'),
 '05_macrophage_communication':('巨噬注释、通讯与NicheNet','方法1.5；结果2.3–2.4；讨论3.2；图3；表2、S9–S11'),
 '06_source_composition':('来源、组成与背景敏感性','方法1.4；讨论3.1–3.3；表3–4'),
 '07_bootstrap':('预定关联及供者重采样','方法1.7；结果2.9；图6；表S12'),
 '08_candidate_diagnostics':('独立候选、匹配背景与删除敏感性','方法1.6；结果2.6–2.7；讨论3.3；表3'),
 '09_mouse_crosscohort':('小鼠与既有补充队列方向对照','方法1.2、1.6、1.8–1.9；结果2.4、2.8、2.13；图S1–S2'),
 '10_report_figures':('整合报告、冻结UMAP与GSEA绘图','全文、图1–10、S1–S4及所有表格'),
 '11_delivery_revision':('本次审阅修订与交付核验','方法1.5–1.7；结果2.6；讨论3.2–3.4；表2–4')}
selected={}; exclusions=[]
def add(f,group,role='analysis',version_note=''):
 if not f.is_file() or f.is_symlink():return
 assert f.resolve().is_relative_to(P.resolve())
 selected[str(f.relative_to(P))]=(group,role,version_note)
def tree(rel,group,role='analysis',version_note=''):
 root=P/rel
 if not root.exists():raise FileNotFoundError(root)
 for f in root.rglob('*'):
  r=f.relative_to(root)
  if any(x in {'__pycache__','.library_helpers','reproduce','source','context'} or x.startswith('.') for x in r.parts):continue
  if f.suffix.lower() in CODE:add(f,group,role,version_note)
  elif f.suffix.lower() in {'.txt','.json','.yaml','.yml'} and any(x in f.name.lower() for x in ['spec','sessioninfo','methods','environment']):
   if f.stat().st_size<200000 and 'complete' not in f.name.lower():add(f,group,'specification_environment',version_note)

for f in (P/'scripts').iterdir():
 if f.suffix.lower() not in CODE:continue
 n=f.name
 group='01_qc_annotation'
 if n.startswith(('09','10','11')):group='02_cd4_pseudobulk'
 if n.startswith(('12','14','15','16','17','18','lung_mechanism')):group='03_function_programs'
 add(f,group,'analysis_or_original_launcher','QC v3 is the corrected audit; older QC and repair launchers document the historical input-recovery path, not an alternative current endpoint.' if n.startswith(('01','02')) else '')
analysis={
 'mechanism_step1_20260927_152218_s4V5eF':'03_function_programs',
 'mechanism_step2_20260927_154919_EsN4gk':'03_function_programs',
 'mechanism_step3_20260928_103057_lIRL38':'05_macrophage_communication',
 'mechanism_step4_20260928_150337_J58beC':'08_candidate_diagnostics',
 'mechanism_step5_20260928_153825_Q4sB1x':'05_macrophage_communication',
 'mechanism_step6_20260928_163015_VoRQqA':'08_candidate_diagnostics',
 'mechanism_step7_20260929_100321_6hm4dZ':'08_candidate_diagnostics',
 'mechanism_step8_20260929_101407_jWFgLl':'05_macrophage_communication',
 'scenic1_20260929_103513_xRoEnr':'04_scenic',
 'scenic2_20261001_092655_ZKNWD6':'04_scenic',
 'scenic3_20261001_123452_Kvp9h5':'05_macrophage_communication',
 'scenic4_20261002_091545_iazXf1':'04_scenic',
 'GSE314943':'09_mouse_crosscohort','GSE183682':'09_mouse_crosscohort',
 'GSE274372':'09_mouse_crosscohort','GSE174725':'09_mouse_crosscohort',
 'screen2_TNK_20260911_103239_nfW0yu':'01_qc_annotation',
 'screen2_Trefine_nuzVQu':'01_qc_annotation','screen2_CD4_target_wz1lfx':'01_qc_annotation',
 'screen6_G174_recover_8E2LLp':'09_mouse_crosscohort',
 'background_sensitivity_20260914_110753_tQJZlp':'06_source_composition'}
for rel,g in analysis.items():tree('02_analysis/'+rel,g,version_note='Frozen historical analysis; retain distinct cohorts, subtypes and scoring definitions. Late mechanism6/7 branches selected; earlier superseded branches excluded.')
reports={
 'human_first_20261004_v1':'05_macrophage_communication',
 'human_alignment_20261005_v1':'08_candidate_diagnostics',
 'human_gap_diagnostics_20261005_v1':'08_candidate_diagnostics',
 'IKZF2_IL7R_target_deletion_20261005_v1':'08_candidate_diagnostics',
 'phenotype_candidate_discovery_20261005_v1':'08_candidate_diagnostics',
 'independent_candidates_revision_20261005_v1':'08_candidate_diagnostics',
 'bootstrap_TNF_IL18_20261004_v1':'07_bootstrap',
 'breakpoint_TNF_IL18_20261004_v1':'07_bootstrap',
 'dual_line_TNF_IL18_20261004_v1':'07_bootstrap',
 'mouse_sender_20261002_v1':'05_macrophage_communication',
 'missing_links_20261003_v1':'05_macrophage_communication',
 'aligned_evidence_20261002_v1':'06_source_composition',
 'stepwise_diagnostic_20261002_v1':'06_source_composition',
 'front_half_20261004_v1':'06_source_composition'}
for rel,g in reports.items():tree('11_reports/'+rel,g)
tech='11_reports/technical_delivery_20261002_ikzf2_stat1_v1/revision_v2/final_release'
tree(tech+'/code','10_report_figures','historical_report_assembly')
tree(tech+'/prior_review_repaired/scenic3_review','05_macrophage_communication','reviewed_sender_annotation')
tree(tech+'/provenance','10_report_figures')
for rel in ['phase_delivery_IKZF2_STAT1_20261006_v2','manuscript_integrated_20261006_v3','manuscript_integrated_20261006_v4','manuscript_graphics_20261006_v5']:
 tree('11_reports/'+rel+'/code','10_report_figures','report_assembly','Different report-generation layers are retained because the final builders import earlier document, figure and table functions.')
# The current report snapshot carries the source-code dependencies of the earlier builders.
for f in (P/'11_reports/manuscript_graphics_20261006_v5/reproduce').glob('*.py'):
 if any(f.name.startswith(x) for x in ['manuscript_','revision_','narrative_','graphics_']) and not any(x in f.name for x in ['library','upload','transfer','inventory','locate','inspect','server','package','close','diagnose']):add(f,'10_report_figures','report_dependency')
for root in [P/'11_reports/manuscript_graphics_20261006_v5/provenance',P/'11_reports/human_first_20261004_v1/provenance']:
 for f in root.glob('*'):
  if f.is_file() and ('sessionInfo' in f.name or f.name=='GSEA_display_spec.json'):add(f,'10_report_figures','environment_specification')
for f in (O/'code').glob('*'):
 if f.suffix.lower() in CODE:add(f,'11_delivery_revision','current_delivery')
# Non-sensitive definitions only: no sample, donor or cell result tables.
for rel in ['02_analysis/mechanism_step1_20260927_152218_s4V5eF/tables/module_definition.tsv','00_project/PROJECT_SCOPE.md']:
 if (P/rel).exists():add(P/rel,'03_function_programs','definition')

reject_names=['library','upload','transfer_status','inventory','inspect','locate','word_status','word_windows']
for rel in list(selected):
 n=Path(rel).name.lower()
 if any(x in n for x in reject_names):exclusions.append({'original_path':rel,'reason':'delivery transport/discovery utility, not scientific workflow'});del selected[rel]

patterns={
 'private_key':re.compile(r'-----BEGIN (?:RSA |OPENSSH |EC |DSA )?PRIVATE KEY-----'),
 'credential_url':re.compile(r'(?:https?|socks5)://[^\s/\"\']+:[^\s/@\"\']+@'),
 'known_token':re.compile(r'\b(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{35,}|sk-[A-Za-z0-9_-]{30,}|AKIA[0-9A-Z]{16})\b'),
 'literal_secret':re.compile(r'(?i)\b(?:password|passwd|api_key|access_token|client_secret|proxy_password)\s*[:=]\s*[\"\']([^\"\'\n]{8,})[\"\']'),
 'signed_url':re.compile(r'(?i)(?:[?&](?:sig|signature|X-Amz-Signature)=)[A-Za-z0-9%+/=_-]{15,}')}
records=[];byhash={};scan=[];statics=[];refs=[];imports=set()
for rel,(group,role,note) in sorted(selected.items(),key=lambda x:(x[1][0],x[0])):
 f=P/rel;raw=f.read_bytes();sha=hashlib.sha256(raw).hexdigest();s=raw.decode('utf-8-sig')
 for kind,rx in patterns.items():
  for m in rx.finditer(s):
   # Variable placeholders and Python dictionary keys are not literal credentials.
   if kind=='literal_secret' and any(x in m.group(1).lower() for x in ['environ','placeholder','your_','getenv','none','example','replace']):continue
   scan.append({'path':rel,'line':s[:m.start()].count('\n')+1,'kind':kind})
 if any(x['path']==rel for x in scan):continue
 duplicate=sha in byhash
 dest=byhash[sha] if duplicate else 'code/'+group+'/'+rel
 if not duplicate:
  target=C/dest;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw);byhash[sha]=dest
  if f.suffix.lower()=='.py':
   try:
    tree_=ast.parse(s,filename=rel)
    for n in ast.walk(tree_):
     if isinstance(n,ast.Import):imports.update(a.name.split('.')[0] for a in n.names)
     if isinstance(n,ast.ImportFrom) and n.module:imports.add(n.module.split('.')[0])
    statics.append({'file':dest,'check':'python_ast','status':'PASS'})
   except SyntaxError as e:statics.append({'file':dest,'check':'python_ast','status':'FAIL','line':e.lineno,'message':e.msg})
  elif f.suffix.lower()=='.sh':
   r=subprocess.run(['bash','-n',str(target)],capture_output=True,text=True)
   statics.append({'file':dest,'check':'bash_n','status':'PASS' if r.returncode==0 else 'FAIL','message':r.stderr[:500]})
  for i,line in enumerate(s.splitlines(),1):
   # Save dependencies as paths only; do not export whole source-data rows.
   for literal in re.findall(r'[\"\']([^\"\'\n]{3,400})[\"\']',line):
    if re.search(r'(?:02_analysis/|11_reports/|/mnt/|\.rds$|\.h5ad$|\.feather$|\.loom$|\.gmt$|\.py$|\.R$|\.sh$)',literal) and '\t' not in literal:
     refs.append({'code_file':dest,'line':i,'path_literal':literal,'classification':'external_resource' if literal.startswith('/mnt/') and not literal.startswith(str(P)) else 'project_or_relative_dependency'})
 records.append({'category':group,'role':role,'original_path':str(f),'package_path':dest,'bytes':len(raw),'original_sha256':sha,'package_sha256':sha,'duplicate_alias':duplicate,'version_note':note})
if scan:
 (O/'sensitive_review_needed.json').write_text(json.dumps(scan,indent=2))
 print(json.dumps({'status':'REVIEW_REQUIRED','findings':scan}));raise SystemExit(2)
def tsv(rows,name):
 if not rows:return
 keys=list(dict.fromkeys(k for d in rows for k in d))
 with (C/'provenance'/name).open('w',encoding='utf-8',newline='') as h:
  w=csv.DictWriter(h,fieldnames=keys,delimiter='\t');w.writeheader();w.writerows(rows)
tsv(records,'code_manifest.tsv');tsv(exclusions,'excluded_utilities.tsv');tsv(refs,'literal_dependency_index.tsv');tsv(statics,'static_checks.tsv')
counts=collections.Counter(r['category'] for r in records if not r['duplicate_alias'] and Path(r['package_path']).suffix.lower() in CODE)
summary={'source_project':str(P),'output':str(O),'source_files':len(records),'unique_files':len(byhash),'duplicate_aliases':sum(r['duplicate_alias'] for r in records),'unique_code_files':sum(counts.values()),'category_code_counts':dict(counts),'credential_scan_findings':scan,'static_failures':[x for x in statics if x['status']!='PASS'],'new_analysis_runs':0,'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(C/'provenance/package_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
(C/'documentation/python_import_inventory.txt').write_text('\n'.join(sorted(imports))+'\n',encoding='utf-8')
tsv([{'category':g,'purpose_zh':v[0],'report_sections_figures':v[1],'unique_code_files':counts.get(g,0)} for g,v in groups.items()],'code_report_mapping.tsv')
tsv([{'original_path':r['original_path'],'package_path':r['package_path'],'category':r['category'],'report_sections_figures':groups[r['category']][1]} for r in records],'per_file_report_mapping.tsv')
for f in (O/'documentation').glob('*'):
 if f.is_file():shutil.copy2(f,C/('README_zh.md' if f.name=='review_README_zh.md' else 'documentation/'+('运行顺序与依赖.md' if f.name=='review_RUN_ORDER_zh.md' else f.name)))
print(json.dumps(summary,ensure_ascii=True))
