"""Assemble delivery provenance without modifying any of the original 285 scripts."""
from pathlib import Path
import json,hashlib,shutil,zipfile,csv,subprocess,sys,ast
import pymupdf as fitz
B=Path(__file__).resolve().parent;D=B/'r1r2_revision_20261007';K=D/'package';C=K/'report_build/current';H=K/'report_build/history/library_v7'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,value):p.write_text(json.dumps(value,ensure_ascii=False,indent=2),encoding='utf-8')
pdf=K/'report/CD4T_review_R1R2_20261007.pdf';doc=fitz.open(pdf)
obsolete=K/'report/CWP_Silica_CD4_IKZF2_STAT1_20261006_review.pdf'
if obsolete.exists():
    assert obsolete.resolve().parent==(K/'report').resolve() and sha(obsolete)==sha(H/'outputs/library_v7.pdf')
    obsolete.unlink()
for name in ['r1r2_prepare_inputs.py','r1r2_build_pdf.py','r1r2_package.py','r1r2_validate_package.py']:shutil.copy2(B/name,C/'code'/name)
parity=json.loads((C/'QA/PDF_content_parity.json').read_text(encoding='utf-8'))
assert parity['pages']==len(doc) and parity['pdf_sha256']==sha(pdf)
visual={'pdf_sha256':sha(pdf),'pages':len(doc),'all_pages_reviewed':list(range(1,len(doc)+1)),'contact_sheets_reviewed':4,'enlarged_target_page':15,'clipped_text_or_table':False,'non_target_pixel_comparison':'29 complete pages identical; page 15 identical outside Table 2','new_statistical_calculations':0}
write(C/'QA/visual_review.json',visual)
for name in ['PDF_content_parity.json','table2_source_validation.json','non_target_invariance.json']:shutil.copy2(C/'QA'/name,K/'provenance'/name)
shutil.copy2(C/'revision_manifest.json',K/'provenance/revision_manifest.json')
shutil.copy2(C/'QA/table2_source_validation.json',K/'provenance/numeric_revision_checks.json')
write(K/'provenance/QA_result.json',{'pdf_sha256':sha(pdf),'pages':len(doc),'tables':18,'figures':14,'text_blocks_and_cells_verified':parity['html_text_blocks_and_cells_verified'],'missing_text_blocks':parity['missing_text_blocks'],'visual_review':visual,'non_target_content_unchanged':True,'scope':'Current R1/R2 PDF; old 29-page QA is explicitly historical'})
receipt=D/'library_finalization_receipt.json'
assert receipt.exists(),'Wait for guarded Library finalization before final packaging'
library=json.loads(receipt.read_text(encoding='utf-8'))['structuredContent']['results'][0]
assert library['status']=='succeeded' and library['library_file_id']=='libfile_51bec083de2c81919ea20e968043ae82' and library['current_version_number']==8
shutil.copy2(receipt,C/'QA/library_finalization_receipt.json')
release={'release':'20261007_R1R2','pdf':parity,'library_file_id':library['library_file_id'],'library_current_version':library['current_version_number'],'based_on_library_version':7,'repair_scope':['R1 ZIP provenance','R2 Table 2 rho'],'unchanged':['R3 diagnostic P/BH policy','G288 interpretation sentence','title','Word','all non-target manuscript content','original 285 scripts'],'current_build_root':'report_build/current','historical_v7_root':'report_build/history/library_v7','historical_v5_root':'report_build/history/library_v5','new_scientific_models':0}
write(K/'provenance/CURRENT_RELEASE.json',release)
readme='''# CD4 T细胞报告：复审R1/R2交付

本包对应2026-10-07复审R1/R2修正版，原Library PDF第8版。正式报告为 `report/CD4T_review_R1R2_20261007.pdf`；准确页数、字节数及SHA256见 `provenance/CURRENT_RELEASE.json`。PDF仅在表2增加偏秩相关ρ列，其余正文、标题、图件和原n/P/BH均保持。没有追加24条诊断的BH，没有改动G288零检出解释，没有更新Word。

## 当前验收入口

- `provenance/CURRENT_RELEASE.json`：当前发布身份和版本关系。
- `provenance/PDF_content_parity.json`、`QA_result.json`、`package_validation.json`：均对应当前PDF；页数由实际PDF读取。
- `report_build/current/`：本次真实编辑/构建脚本、源表六行聚合摘录、最终HTML文字源、表2局部PDF、逐页校样及核验结果。
- `report_build/history/library_v7/`：30页基础PDF、实际最终编辑/渲染脚本、94条文字变更、输入/输出HTML及1643文本块核验等既有记录。
- `report_build/history/library_v5/`：29页版本的QA和说明，只用于历史溯源，不作为当前验收。
- `documentation/最终文稿构建与版本关系.md`：版本到输入、代码、输出与QA的完整对应及重建入口。
- `FILE_MANIFEST_SHA256.tsv`：本包所有其他文件的大小和SHA256。

## 原科学代码保持

`code/`内原285份脚本及全部配置逐字节保留。报告构建和核验补充放在独立的 `report_build/`，不计入原285份脚本。原代码来源索引450行及去重关系仍在 `provenance/code_manifest.tsv`；原静态解析记录明确属于历史科学代码核验，不冒充本次完整流程运行。`documentation/运行顺序与依赖.md`区分科学流程、29页文稿构建与当前PDF构建。

原项目：`/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA`。人类病例为CWP；小鼠主队列为GSE314943硅尘模型。IKZF2与STAT1仍为独立候选。

## 复现边界

本包是可追溯交付，不是携带全部数据和数据库的一键运行镜像。原始计数、FASTQ、Seurat/RDS/H5AD/loom、SCENIC基序库、NicheNet先验、冻结调控子和科学图表输入仍依赖服务器；不含个体级数据和凭据。当前PDF的局部修订输入已随包保存，可按文稿构建说明重建表2及核验。

GSE314943早期stage1–5部分独立执行脚本未在既有源码盘点中找到；stage5 CD4对象及元数据仍为外部输入，未补写假冒历史脚本。动物/混样—文库映射、人类MSigDB确切发布号及早期UCell默认maxRank的追溯限制仍保留。没有重跑聚类、SCENIC、NicheNet、Bootstrap或科学模型。

原代码可能含创建目录、安装依赖或覆盖输出的历史启动器，应先阅读并配置独立运行目录。静态解析、文件完整性与当前PDF核验不能代替完整科学流程复现。

只读验收：

```text
python report_build/current/code/r1r2_validate_package.py --package <解压后的包根目录>
```

该命令核对当前文件清单、PDF身份、表2行值、285份原脚本及路径/敏感模式，不执行分析脚本。
'''
(K/'README_zh.md').write_text(readme,encoding='utf-8')
oldorder=(K/'report_build/history/library_v5/运行顺序与依赖.md').read_text(encoding='utf-8')
oldorder=oldorder.replace('| 11 本次修订 |','| 11 历史29页文稿（Library v5） |').replace('29页PDF、18表；终点及证据分层修订；代码包，无新科学计算','历史29页PDF、18表；仅解释该版本来源，不作为当前验收')
(K/'documentation/运行顺序与依赖.md').write_text('# 当前入口\n\n当前R1/R2 PDF构建与验收见 `最终文稿构建与版本关系.md`。以下保留原科学流程和历史文稿入口；第11项属于29页版本。原285份脚本不改写。\n\n'+oldorder,encoding='utf-8')
guide='''# 最终文稿构建与版本关系

| 版本 | 输入与实际代码 | 输出与验收 |
|---|---|---|
| 历史Library v5（29页） | `code/11_delivery_revision/.../review_revision.py`等原脚本 | `report_build/history/library_v5/QA/`，不用于当前验收 |
| Library v7（30页） | `report_build/history/library_v7/inputs/`；同目录`code/logic_edit_report.py`、`logic_render_pdf.py`、`logic_finalize_package.py` | `outputs/logic_render_source.html`、`outputs/library_v7.pdf`、`outputs/editorial_changes.json`（94条）；`QA/PDF_content_parity.json`（1643块）及逐页校样 |
| 当前Library v8（R1/R2） | `report_build/current/inputs/IL15_table2_source_rows.json`与v7基础PDF；`current/code/r1r2_build_pdf.py` | `report/CD4T_review_R1R2_20261007.pdf`；`current/sources/final_report.html`、`table2_overlay.pdf`及`current/QA/` |

原v7脚本按当时工作区布局原样归档，不能将其目录位置误当作已适配安装路径。重建v7文字时，在独立工作目录放置这三份脚本；将`inputs/review_render_source.html`复制至`review_delivery/provenance/`，将`inputs/original_blocks.json`复制至`logic_revision_20261007/provenance/`，建立`logic_revision_20261007/QA/`，运行`logic_edit_report.py`。该编辑器的完整输出与94条变更均已提供。`logic_render_pdf.py`使用Windows Edge及系统字体；PDF元数据和不同渲染器版本可能影响字节哈希，正式发布以保存的PDF哈希为准。`logic_finalize_package.py`是当时打包脚本，依赖上一版ZIP；它不是当前验收入口。

当前重建只更新表2，不重渲染正文：

```text
python report_build/current/code/r1r2_build_pdf.py --package <独立解压副本>
```

需要Python、lxml、PyMuPDF、numpy、Pillow、ReportLab，以及Windows的SimSun/Times New Roman字体；使用字体及输入哈希列于`current/revision_manifest.json`。此命令会更新副本内当前PDF、HTML及QA，不应直接对归档原件运行。它验证v7基础PDF哈希，按完整键匹配IL15六行、核对原n/P/BH，仅生成表2局部页面并替换原区域；随后核对1650个文本块/表格单元格及逐页像素。当前最终PDF无需重新运行科学模型。

`r1r2_prepare_inputs.py`记录本次从既有工作区归档代码/QA并经SSH只读取得源表摘录的实际操作；依赖原工作区和服务器，不属于脱离原环境的一键下载器。`r1r2_package.py`记录实际打包过程，同样依赖工作区、基线ZIP和已完成Library回执。独立包内只读检查使用`r1r2_validate_package.py`。

旧29页的断言仍是对应历史脚本的一部分，为保持原285份脚本字节不变而保留；当前验证使用实际PDF读出的页数和哈希。当前与历史验收记录均通过`version_file_inventory.tsv`及`FILE_MANIFEST_SHA256.tsv`绑定文件身份。

仅处理R1、R2；24条诊断的P/FDR口径和G288句子均未改。四个新增ρ从既有结果摘录，不是重新拟合；MacAM的ρ写“未计算”。
'''
(K/'documentation/最终文稿构建与版本关系.md').write_text(guide,encoding='utf-8')
write(K/'provenance/package_summary.json',{'source_project':'/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA','current_release':'20261007_R1R2','original_unique_code_files':285,'original_source_index_rows':450,'original_code_unchanged':True,'historical_v7_build_scripts':3,'current_build_scripts':4,'new_scientific_analysis_runs':0,'historical_package_summary':'report_build/history/library_v5/QA/package_summary.json'})
baselinezip=zipfile.ZipFile(next((B/'logic_revision_20261007').glob('*.zip')))
codefiles=[p for p in (K/'code').rglob('*') if p.is_file()]
assert all(p.read_bytes()==baselinezip.read(p.relative_to(K).as_posix()) for p in codefiles)
originalscripts=[p for p in codefiles if p.suffix.lower() in ['.py','.r','.sh','.ps1']];assert len(originalscripts)==285
word=B/'graphics_delivery/CD4_IKZF2_STAT1_integrated_manuscript_zh.docx';assert sha(word)=='1bd05f6e28a5c60da9a012b5397f237e0999919d43e0890b893074314882476c'
validation={'pdf_relative_path':pdf.relative_to(K).as_posix(),'pdf_pages':len(doc),'pdf_sha256':sha(pdf),'pdf_bytes':pdf.stat().st_size,'code_files':285,'all_original_code_files_byte_identical':True,'current_report_build_scripts':4,'historical_v7_build_scripts':3,'source_index_rows':450,'all_pdf_pages_reviewed':len(doc),'new_scientific_models':0,'new_P_or_BH_computations':0,'table2_values_matched_by_full_key':6,'old_word_unchanged':True,'historical_static_checks':'Existing Python/Shell/R/PowerShell records retained for original scientific code; not a full pipeline run','remaining_reproduction_gaps':['GSE314943 early stage1-5 standalone scripts not all located','External raw data/frozen objects/model databases required for scientific recomputation','Animal-library mapping not closed'],'delivery_blockers':[]}
write(K/'provenance/package_validation.json',validation)
inventory=K/'provenance/version_file_inventory.tsv'
with inventory.open('w',encoding='utf-8',newline='') as f:
    w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['version_scope','relative_path','bytes','sha256'])
    for p in sorted((K/'report_build').rglob('*')):
        if p.is_file():
            scope='current_R1R2' if p.is_relative_to(C) else ('historical_Library_v7' if p.is_relative_to(H) else 'historical_Library_v5')
            w.writerow([scope,p.relative_to(K).as_posix(),p.stat().st_size,sha(p)])
manifest=K/'FILE_MANIFEST_SHA256.tsv'
with manifest.open('w',encoding='utf-8',newline='') as f:
    w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['relative_path','bytes','sha256'])
    for p in sorted(K.rglob('*')):
        if p.is_file() and p!=manifest:w.writerow([p.relative_to(K).as_posix(),p.stat().st_size,sha(p)])
r=subprocess.run([sys.executable,str(C/'code/r1r2_validate_package.py'),'--package',str(K)],capture_output=True,text=True,encoding='utf-8');assert r.returncode==0,r.stderr
print(r.stdout,flush=True)
out=D/'矽肺研究中巨噬细胞与CD4T细胞表型关联分析_PDF及代码_20261007_复审R1R2修正版.zip'
with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in sorted(K.rglob('*')):
        if p.is_file():z.write(p,p.relative_to(K).as_posix())
with zipfile.ZipFile(out) as z:
    assert z.testzip() is None and z.read(pdf.relative_to(K).as_posix())==pdf.read_bytes()
    for p in originalscripts:assert z.read(p.relative_to(K).as_posix())==baselinezip.read(p.relative_to(K).as_posix())
    result={'zip_name':out.name,'zip_bytes':out.stat().st_size,'zip_sha256':sha(out),'zip_members':len(z.namelist()),'max_member_path_length':max(map(len,z.namelist())),'pdf':parity,'original_scripts_unchanged':285,'package_validator':json.loads(r.stdout),'Library_version':8}
write(D/'delivery_package_result.json',result)
print(json.dumps(result,ensure_ascii=True),flush=True)
