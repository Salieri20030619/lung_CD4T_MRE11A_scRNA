from pathlib import Path
import json,hashlib,shutil,datetime
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');V=P/'11_reports/manuscript_integrated_20261006_v4';O=P/'11_reports/manuscript_graphics_20261006_v5'
assert O.resolve().is_relative_to(P.resolve())
for f in (V/'figures').iterdir():
 if f.is_file():shutil.copy2(f,O/'figures'/f.name)
shutil.copytree(V/'source_tables',O/'source_tables',dirs_exist_ok=True)
shutil.copytree(V/'reproduce',O/'reproduce',dirs_exist_ok=True)
N=O/'reproduce/narrative_delivery';N.mkdir(exist_ok=True)
for name in ['CD4_IKZF2_STAT1_integrated_manuscript_zh.docx','CD4_IKZF2_STAT1_integrated_manuscript_zh.pdf']:
 shutil.copy2(V/'deliverables'/name,N/name)
shutil.copytree(V/'figures',N/'figures',dirs_exist_ok=True)
for name in ['document_manifest.json','old_blocks.json','added_evidence_aggregate.json']:
 if (V/'provenance'/name).exists():shutil.copy2(V/'provenance'/name,N/name)
for f in (O/'code').iterdir():
 if f.is_file():shutil.copy2(f,O/'reproduce'/f.name)
G=O/'reproduce/graphics_delivery'
for d in ['provenance','figures','QA']:(G/d).mkdir(parents=True,exist_ok=True)
for name in ['old_blocks.json','added_evidence_aggregate.json']:
 shutil.copy2(N/name,G/name)
for f in (O/'provenance').glob('*.json'):
 if f.name in ['GSEA_figure_manifest.json','GSEA_display_spec.json','UMAP_QA.json','GSEA_input_audit.json','baseline_v3.json']:shutil.copy2(f,G/'provenance'/f.name)
shutil.copytree(O/'figures',G/'figures',dirs_exist_ok=True)
inputs=[]
for f in sorted((O/'tables').iterdir()):
 if f.is_file():inputs.append({'path':str(f.relative_to(O)),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
orig=[]
for path in ['04_preprocessing/GSE288887_broad_clustered_scored.rds','06_CD4T_analysis/CD4_final/GSE288887_CD4_final.rds','07_DEG/CD4_pseudobulk_v1/pathway_v1/GSE288887_CD4_pathway_results.rds','10_tables/CD4_final/pathway_v1/GO_KEGG_GSEA.tsv','02_analysis/mechanism_step2_20260927_154919_EsN4gk/input/human_cells.tsv']:
 f=P/path;orig.append({'path':str(f),'bytes':f.stat().st_size,'mtime_utc':datetime.datetime.fromtimestamp(f.stat().st_mtime,datetime.timezone.utc).isoformat(),'sha256':hashlib.sha256(f.read_bytes()).hexdigest() if f.stat().st_size<50_000_000 else None})
(O/'provenance/plot_input_manifest.json').write_text(json.dumps({'frozen_plot_inputs':inputs,'original_objects':orig,'large_original_object_hash_note':'Large Seurat objects referenced by path, size and mtime; extracted frozen metadata/coordinates used for plotting are fully hashed and archived.'},indent=2),encoding='utf8')
(O/'README_复现.txt').write_text('本次仅UMAP与GSEA图形修订，无新降维、聚类、疾病模型或富集检验。\n原始冻结对象位置与绘图输入哈希见provenance/plot_input_manifest.json。tables含实际坐标、完整绘图对象和ES核验，仅保存在服务器。\n服务器绘图：既有R_4.5.2环境执行code/graphics_umap.R和code/graphics_gsea.R。前者读取冻结坐标，后者读取原gseaResult对象；Top5展示规则见GSEA_display_spec.json。\nWindows报告：reproduce目录依次运行graphics_build.py、graphics_pdf_render.py、graphics_pdf_finalize.py、graphics_audit.py和graphics_qa.py。PDF由本机Edge无头渲染精确提取的DOCX段落、表格和内嵌图片，再核验1490文本块一致性。graphics_render.ps1为可选Word路线，本次因Word COM异常耗时未采用；未修改Word或系统全局配置。\nGSEA示例脚本仅供阅读和风格参考，位于provenance，禁止作为本项目分析脚本执行。\n原分析表、负结果与正结果均保留；本次报告及关键高清图上传Library，细胞明细、矩阵和ZIP不上传。\n',encoding='utf8')
print(json.dumps({'plot_inputs':len(inputs),'original_references':len(orig),'archive_prepared':True}))
