from pathlib import Path
import json,shutil,hashlib,zipfile,datetime
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');V=P/'11_reports/manuscript_integrated_20261006_v3';O=P/'11_reports/manuscript_integrated_20261006_v4'
assert O.resolve().is_relative_to(P.resolve())
shutil.copytree(V/'source_tables',O/'source_tables',dirs_exist_ok=True)
shutil.copytree(V/'reproduce',O/'reproduce',dirs_exist_ok=True)
shutil.copy2(V/'deliverables/CD4_IKZF2_STAT1_integrated_manuscript_zh.docx',O/'reproduce/manuscript_delivery/CD4_IKZF2_STAT1_integrated_manuscript_zh.docx')
shutil.copytree(V/'figures',O/'reproduce/manuscript_delivery/figures',dirs_exist_ok=True)
for f in (V/'provenance').iterdir():
 if f.is_file() and ('source_manifest' in f.name or 'methods' in f.name or f.name in ['TF_program_BH_sizes.json','enrichment_inventory.json','added_evidence_aggregate.json']):shutil.copy2(f,O/'provenance'/f.name)
for f in (O/'code').iterdir():
 if f.is_file():shutil.copy2(f,O/'reproduce'/f.name)
N=O/'reproduce/narrative_delivery';N.mkdir(exist_ok=True)
for f in (O/'provenance').iterdir():
 if f.is_file() and f.name in ['old_blocks.json','baseline.json','added_evidence_aggregate.json']:shutil.copy2(f,N/f.name)
shutil.copytree(O/'figures',N/'figures',dirs_exist_ok=True)
qa=json.loads((O/'provenance/QA_result.json').read_text(encoding='utf8'));audit=json.loads((O/'provenance/editorial_audit.json').read_text(encoding='utf8'))
assert not qa['out_of_bounds'] and not qa['internal_term_hits'] and not qa['empty_BH_columns']
assert audit['tables_identical']==15 and audit['source_data_unchanged'] and not audit['old_prose_numerals_for_review']
for name,sha in qa['sha256'].items():assert hashlib.sha256((O/'deliverables'/name).read_bytes()).hexdigest()==sha
(O/'README_中文.txt').write_text('尘肺CD4整合文稿：叙述顺序修订\n\n结果按免疫全景与CD4表型、巨噬通讯、IKZF2/STAT1独立候选、功能程序组织；阴性结果和限制集中讨论，完整数值保留补充表。\n此轮未重新分析。15张表与Library版本2相同，既有统计源数据哈希保持不变。\n\n正式文稿位于deliverables，9幅使用的图位于figures，QA及版本回执位于provenance，原有完整科学源表位于source_tables。\n在Windows及原有Python/Word环境，进入reproduce目录运行python narrative_build.py，再运行narrative_render.ps1；python narrative_qa.py、python narrative_audit.py、python narrative_numeric_audit.py复核排版、数据保留和数值。\n图重建使用既有聚合结果；新组合图由python narrative_figures.py生成，其余图沿用前一版本。\n原始表达矩阵与正式先验留在原项目位置，此ZIP不上传Library。\n',encoding='utf8')
status={'date_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'server':'osv01','output':str(O),'pages':qa['pages'],'figures':9,'tables':15,'references':15,'new_scientific_models':0,'numeric_checks':'34/34','table_data_preserved':15,'prior_prose_numerals_preserved':True,'library_status':'pending'}
(O/'provenance/FINAL_STATUS.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf8')
z=O/'deliverables/CD4_narrative_reproduction_server_only.zip'
with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as zz:
 for f in sorted(O.rglob('*')):
  if f.is_file() and f!=z and f.suffix!='.zip':zz.write(f,f.relative_to(O))
status['zip_bytes']=z.stat().st_size;status['zip_sha256']=hashlib.sha256(z.read_bytes()).hexdigest();(O/'provenance/FINAL_STATUS.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf8');print(json.dumps(status))
