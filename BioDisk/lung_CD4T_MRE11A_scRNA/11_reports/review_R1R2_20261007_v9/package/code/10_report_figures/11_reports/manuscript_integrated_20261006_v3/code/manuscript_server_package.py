from pathlib import Path
import json,hashlib,shutil,zipfile,datetime
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/manuscript_integrated_20261006_v3';V=P/'11_reports/phase_delivery_IKZF2_STAT1_20261006_v2';H=P/'11_reports/human_first_20261004_v1/results'
assert O.resolve().is_relative_to(P.resolve())
for sub in ['source_tables/human_communication','source_tables/human_pathways','source_tables/mouse_pathways','source_tables/bootstrap','source_tables/regulons','reproduce/manuscript_delivery','reproduce/revision_delivery','provenance','deliverables','code','figures']:(O/sub).mkdir(parents=True,exist_ok=True)
manifest=[]
def copy(f,to):
 f=Path(f);to=Path(to)
 if not f.is_file():return
 to.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(f,to)
 manifest.append({'source':str(f),'archive_path':str(to.relative_to(O)),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
for f in H.iterdir():
 if f.is_file() and f.suffix in ['.gz','.tsv','.json','.txt'] and not f.name.startswith('stage3_all_mapped_gene_directions'):copy(f,O/'source_tables/human_communication'/f.name)
for root in [P/'10_tables/CD4_final/pathway_v1',P/'10_tables/CD4_final/pathway_sensitivity_v1',P/'10_tables/CD4_final/module_audit_v1',P/'10_tables/CD4_final/subtype_module_audit_v1']:
 for f in root.glob('*.tsv'):copy(f,O/'source_tables/human_pathways'/root.name/f.name)
for name in ['stage7_20260905_234200_aiQYL0','stage8_20260906_003153_6mFKmg']:
 root=P/'02_analysis/GSE314943'/name
 for f in (root/'tables').glob('*'):
  if f.is_file() and f.suffix in ['.tsv','.csv','.gz']:copy(f,O/'source_tables/mouse_pathways'/('enrichment' if name.startswith('stage7') else 'rank_scores')/f.name)
 copy(root/'analysis.R',O/'code'/('mouse_enrichment_original.R' if name.startswith('stage7') else 'mouse_rank_scores_original.R'))
for f in (P/'11_reports/bootstrap_TNF_IL18_20261004_v1/results').glob('*.tsv'):copy(f,O/'source_tables/bootstrap'/f.name)
R=next((P/'02_analysis').glob('scenic2_*'))
for name in ['TF_associations.tsv','TF_group_effects.tsv','regulon_exclusion_coverage.tsv','run_summary.tsv','gene_alias_audit.tsv']:copy(R/'tables'/name,O/'source_tables/regulons'/name)
for name in ['analyze.py','stats.py','methods.txt']:copy(R/name,O/'code'/('original_regulon_'+name))
copy(P/'02_analysis/mechanism_step5_20260928_153825_Q4sB1x/tables/geneset_members.tsv',O/'source_tables/mouse_full_target_background_members.tsv')
copy(V/'provenance/report_source_manifest.json',O/'provenance/prior_source_manifest.json')
copy(V/'provenance/methods_readonly_extract.txt',O/'provenance/earlier_methods_readonly_extract.txt')
copy(V/'provenance/report_evidence_aggregate.json',O/'reproduce/revision_delivery/report_evidence_aggregate.json')
copy(V/'provenance/document_manifest.json',O/'reproduce/revision_delivery/document_manifest.json')
copy(O/'provenance/added_evidence_aggregate.json',O/'reproduce/manuscript_delivery/added_evidence_aggregate.json')
for f in (O/'code').glob('manuscript_*'):
 if f.suffix in ['.py','.ps1']:copy(f,O/'reproduce'/f.name)
(O/'provenance/reproduction_source_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf8')
(O/'README_中文.txt').write_text('尘肺CD4整合研究文稿\n\n正式报告位于deliverables，主图及补充图位于figures。\n此次工作只复用并重写既有结果；未新增科学模型、聚类、SCENIC或Bootstrap。\nsource_tables保留完整疾病模型、通讯候选、靶集/背景和配体靶成员、富集与评分结果、原调控子及重采样统计。\nprovenance保留源文件路径、SHA256及质量核验。原始表达矩阵、正式先验和完整大对象仍位于原项目路径，未复制进报告文件。\n\n重建图表与文档：在Windows、Python环境（pandas/numpy/matplotlib/python-docx/PyMuPDF/Pillow）及宋体/黑体/微软雅黑字体可用时，进入reproduce目录，运行python manuscript_figures.py，再运行python manuscript_build.py；Word已安装时运行manuscript_render.ps1导出PDF，随后运行python manuscript_qa.py和python manuscript_numeric_audit.py。\n此复现包重建既有结果的图表和文稿，不重新拟合科学模型；原始科学计算代码和对象位置见源清单。\nLibrary只更新同一PDF及DOCX，不上传个体数据或此复现ZIP。\n',encoding='utf8')
z=O/'deliverables/CD4_integrated_reproduction_server_only.zip'
with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as zz:
 for f in sorted(O.rglob('*')):
  if f.is_file() and f!=z and f.suffix!='.zip':zz.write(f,f.relative_to(O))
status={'date_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'server':'osv01','output':str(O),'source_files_archived':len(manifest),'zip_bytes':z.stat().st_size,'zip_sha256':hashlib.sha256(z.read_bytes()).hexdigest(),'new_scientific_models':0,'library_status':'pending'}
(O/'provenance/FINAL_STATUS.json').write_text(json.dumps(status,indent=2),encoding='utf8');print(json.dumps(status))
