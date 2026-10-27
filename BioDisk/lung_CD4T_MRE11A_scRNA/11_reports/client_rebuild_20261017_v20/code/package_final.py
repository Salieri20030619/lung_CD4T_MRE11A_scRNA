from pathlib import Path
import json,shutil,hashlib,zipfile,subprocess,time
import pandas as pd
O=Path(__file__).parent;H=O.parent;V=H/'lung_append_analysis_20261017_v18';W=H/'lung_feedback_revision_20261017_v19';B=H/'lung_three_cohort_20261017_v16/data';K=O/'review_bundle';K.mkdir(exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def cp(p,rel):
 t=K/rel;t.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,t)
qa=json.loads((O/'report_QA.json').read_text(encoding='utf-8'));assert qa['pages']==31 and qa['figure_images']==22 and not any(c['out_of_bounds'] for c in qa['checks']);qa.update(visual_review='all 31 rendered pages inspected; 22 panels and captions present; no clipping/overlap/missing glyphs observed',status='analysis-side delivery checks passed; new v20 work not independently reviewed');(O/'report_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
c=pd.read_csv(O/'tables/predelivery_checklist.tsv',sep='\t');c.loc[c.ID=='C05','status']='已完成本轮图文及数值验收';c.to_csv(O/'tables/predelivery_checklist.tsv',sep='\t',index=False,encoding='utf-8')
fm=json.loads((O/'figure_manifest.json').read_text(encoding='utf-8'));assert pd.Series([x['module'] for x in fm]).value_counts().eq(2).all()
for i,m in enumerate(fm):m['report_page']=i+5;m['file']='figures/'+m['id']+'.png';m['plot_values']='plot_data/'+m['id']+'.tsv'
pd.DataFrame(fm).to_csv(O/'tables/figure_index.tsv',sep='\t',index=False,encoding='utf-8')
for folder in ['tables','server_tables','plot_data','inputs']:
 for p in (O/folder).rglob('*'):
  if p.is_file():cp(p,folder+'/'+p.relative_to(O/folder).as_posix())
for p in (O/'figures').glob('*'):
 if p.suffix in ['.png','.pdf']:cp(p,'figures/'+p.name)
for p in O.glob('*.py'):cp(p,'code/v20/'+p.name)
cp(O/'enrichment_v20.R','code/v20/enrichment_v20.R')
for name in ['plan.json','review_scope.json','report_QA.json','report_page_index.json','figure_manifest.json','enrichment_complete.json','enrichment_v20.log','remote_export.log','report_text.txt']:
 if (O/name).exists():cp(O/name,name)
for p in (O/'literature').iterdir():
 if p.is_file():cp(p,'literature/'+p.name)
for p in (W/'review_bundle').rglob('*'):
 if p.is_file() and p.suffix!='.pdf' and p.name!='MANIFEST_SHA256.tsv':cp(p,'baseline_v19/'+p.relative_to(W/'review_bundle').as_posix())
TR=O/'sources/11_reports/customer_sequence_20261015_v2'
for p in (TR/'03A_state_pseudotime').iterdir():
 if p.is_file():cp(p,'trajectory/'+p.name)
for name in ['compute_v2.py','report_v2.py']:cp(TR/'code'/name,'code/trajectory/'+name)
for name in ['NicheNet_activity_all.tsv','NicheNet_run_inventory.tsv','expressed_LR_with_decoy_literature.tsv','candidate_complete_routes.tsv']:cp(B/'tables'/name,'v16_network/'+name)
readme='''v20 完整客户报告及复核材料（2026-10-17）
报告：31页，11个图示模块，每模块2张，共22张代表图。
客户完整PDF另存Windows下载目录；本包在工作目录保存，含可编辑代码、高清PNG/PDF、图值与统计来源，不含原始大型细胞矩阵。
优先级：lung_candidate_experiment_literature_20261017.txt优先于v19独立复核的实验排序意见。IL1R2主选、TTC39C损伤探索备选、SKIL状态备选。
本轮交付前核查：tables/predelivery_checklist.tsv、candidate_sample_detection.tsv、candidate_detection_summary.tsv、experimental_feasibility.tsv、literature_evidence_register.tsv。
图号/页码/图意/来源索引：tables/figure_index.tsv。原始图值在plot_data/。
全景/亚群/身份新提取：server_tables/。冻结注释和嵌入未改动，37,155全肺细胞、2,311 Core CD4；轨迹只含1,734 Naive/ActCD4细胞。
新增通路结果：server_tables/ORA_all.tsv和GSEA_all.tsv；完整成员及排名在同目录；tables/GSEA_with_precision_status.tsv标注估计下限。新增通路基于v16冻结joint_stratified DE，不改变既有候选资格。每个原始家族完整保留。
通路背景/共享成员解释：tables/displayed_pathway_member_context.tsv。疾病名称集合不能解释为感染，免疫球蛋白等背景成员不能直接当成T细胞内源作用。
当前候选和原始统计：baseline_v19/。旧PROS1相关输出只作历史结果；当前功能支持已撤回。背景修正不等于已消除全部环境RNA。
轨迹冻结结果与代码在trajectory/及code/trajectory/；本轮复算4个代表程序的64种符号翻转P及4个完整程序BH家族。细胞/bin不是重复。
新增数值核查：tables/v20_numeric_QA.tsv（18项）、trajectory_BH_QA.tsv（4家族），报告渲染检查report_QA.json。原v19独立评审范围见review_scope.json，不能扩展为v20已获独立评审。
关键复算脚本保留服务器原路径；异地复算需映射输入路径。细胞层全套SCENIC不在本包重建范围。
文献以原始研究页面/原文段落及厂家适用性核实；本地全文下载失败记录保留，未伪称全部下载成功。PM2.5研究仅摘要核查，不作全文图级验证。试剂存在不等于原代目标细胞已经验证有效。
推荐复核顺序：图1-4的冻结来源与比例分母；图7-10的背景/成员/完整家族；图17-18的状态范围与样本级统计；图20 TTC39C空缺原因；最终3候选的作用方向与实验可行性。
所有文本UTF-8无BOM，NA不以0代替。MANIFEST_SHA256.tsv核对本包完整性。未向客户或其他第三方发送。
'''
(K/'README_复核入口.txt').write_text(readme,encoding='utf-8')
rows=[dict(path=p.relative_to(K).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in K.rglob('*') if p.is_file() and p.name!='MANIFEST_SHA256.tsv'];pd.DataFrame(rows).to_csv(K/'MANIFEST_SHA256.tsv',sep='\t',index=False,encoding='utf-8')
z=O/'v20_图件统计与复核材料.zip'
with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as zz:
 for p in K.rglob('*'):
  if p.is_file():zz.write(p,p.relative_to(K).as_posix())
with zipfile.ZipFile(z) as zz:
 assert zz.testzip() is None and not any(n.endswith('.zip') for n in zz.namelist())
 for r in rows:assert hashlib.sha256(zz.read(r['path'])).hexdigest()==r['sha256']
d=json.loads((O/'report_delivery.json').read_text(encoding='utf-8'));d.update(support_zip=str(z),support_sha256=sha(z),support_files=len(rows)+1,support_MB=z.stat().st_size/1e6,report_pages=31)
(O/'delivery.json').write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(d,ensure_ascii=False))