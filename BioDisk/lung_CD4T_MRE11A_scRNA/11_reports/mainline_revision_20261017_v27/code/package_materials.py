from pathlib import Path
import json,shutil,pandas as pd,hashlib
O=Path(__file__).parent;H=O.parent;s=(O/'report.py').read_text(encoding='utf-8')
a=s.index("page();p('附录A｜筛选阈值、集合及候选去向'");b=s.index('screening_methods()',a)+len('screening_methods()')
s=s[:a]+"screening_methods()"+s[b:]
s=s.replace("p('候选筛选：证据维度、阈值与集合关系','h')","p('附录A｜筛选阈值、集合及候选去向','h')")
s=s.replace('PROS1相关路线不计入功能性受体支持。','').replace('；PROS1相关功能性受体证据不足','').replace('，PROS1相关功能性受体证据不足','')
s=s.replace('；不作为统一二元交集','').replace('，不作为统一二元交集','')
(O/'report.py').write_text(s,encoding='utf-8')
m=json.loads((O/'figure_manifest.json').read_text(encoding='utf-8'))
for x in m:
 x['caption']=x['caption'].replace('PROS1路线未纳入。','').replace('，PROS1不作为当前功能支持','').replace('，不作为统一二元交集','')
(O/'figure_manifest.json').write_text(json.dumps(m,ensure_ascii=False,indent=2),encoding='utf-8')
P=O/'supplement';(P/'plot_data').mkdir(exist_ok=True);(P/'tables').mkdir(exist_ok=True)
idx=[]
for item in m:
 id=item['id'];src=H/'lung_client_layout_20261017_v22/plot_data'/(id+'.tsv')
 if id in ['F05','F06']:src=H/'lung_cohort_bridge_20261017_v25/plot_data'/(id+'.tsv')
 if id in ['F11','F12','F17']:src=H/'lung_chord_clarity_20261017_v23/plot_data'/(id+'.tsv')
 if id=='F20':src=O/'plot_data/F20.tsv'
 if id in ['F23','F24']:src=H/'lung_client_layout_20261017_v22/plot_data'/('F21.tsv' if id=='F23' else 'F22.tsv')
 if id in ['F21','F22']:continue
 assert src.exists(),src
 shutil.copy2(src,P/'plot_data'/(id+'.tsv'));idx.append({'figure':id,'source':str(src),'sha256':hashlib.sha256(src.read_bytes()).hexdigest()})
# Explicit figure values, validated against existing screening QA.
q=json.loads((H/'lung_screening_display_20261017_v24/screening_QA.json').read_text(encoding='utf-8'));assert q['sequential_genes']==[89,57,56,54,52]
pd.DataFrame({'step':['tested','nominal_DE','DE_function','three_axes','union'],'genes':[9091,2765,111,55,90],'records':[52448,6126,274,155,None]}).to_csv(P/'plot_data/F21.tsv',sep='\t',index=False)
pd.DataFrame({'step':['nominal_DE','function','same_state_signal','source','role'],'genes':q['sequential_genes'],'records':q['sequential_records']}).to_csv(P/'plot_data/F22.tsv',sep='\t',index=False)
for folder,name in [('lung_append_analysis_20261017_v18/tables','associations_and_programs'),('lung_client_rebuild_20261017_v20/tables','report_tables'),('lung_client_rebuild_20261017_v20/server_tables','enrichment_and_figure_sources'),('lung_screening_display_20261017_v24/tables','screening')]:
 dest=P/'tables'/name;dest.mkdir(exist_ok=True)
 for f in (H/folder).glob('*'):
  if f.is_file():shutil.copy2(f,dest/f.name)
for name in ['DE_all_models.tsv','NicheNet_activity_all.tsv','NicheNet_top250_targets.tsv','NicheNet_run_inventory.tsv','expressed_curated_LR.tsv','expressed_LR_with_decoy_literature.tsv','communication_sample_coverage.tsv','candidate_complete_routes.tsv','sample_inventory.tsv']:
 src=H/'lung_three_cohort_20261017_v16/data/tables'/name;shutil.copy2(src,P/'tables'/name)
shutil.copytree(O/'evidence',P/'evidence',dirs_exist_ok=True)
for name in ['protocol_frozen.json','analyze.py']:shutil.copy2(H/'lung_append_analysis_20261017_v18'/name,P/'evidence'/name)
for f in (H/'lung_cohort_bridge_20261017_v25/plot_data').glob('S_*'):shutil.copy2(f,P/'plot_data'/f.name)
for f in (H/'lung_client_layout_20261017_v22/plot_data').glob('GSEA*'):shutil.copy2(f,P/'plot_data'/f.name)
for f in (H/'lung_chord_clarity_20261017_v23/plot_data').glob('trajectory_backbone*'):shutil.copy2(f,P/'plot_data'/f.name)
readme='v27 配套复核材料\n\nplot_data：报告图号对应的图值；F09/F10的Top5选择表另见GSEA_top5_selected.tsv。S_A/S_B为队列示意与补充图谱来源；S_C是概念交集示意，计数见F21及screening表。\ntables：完整统计及原检验家族，含阴性结果。统计原值不因本轮叙事调整重新筛选。\nevidence：SKIL完整50项FDR复算、H2-DM来源解释、接收/发送覆盖及背景指标源代码。H2-DM条目保留在原数据库和模型输出，但报告单列抗原呈递注释，不作为普通表面信号推荐；本轮未重新估计全部网络。\nF21/F22图值为从已核对筛选汇总整理的显示表；完整成员见tables/screening。\nfigure_sources.json记录图值来源，file_inventory.json为文件完整性清单。\n本包为统计和图值复核材料，不包含全部原始测序文件；不将其称为从FASTQ开始重跑的全流程包。\n'
(P/'README.txt').write_text(readme,encoding='utf-8');(P/'figure_sources.json').write_text(json.dumps(idx,ensure_ascii=False,indent=2),encoding='utf-8')
inv=[{'path':str(f.relative_to(P)),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in P.rglob('*') if f.is_file() and f.name!='file_inventory.json'];(P/'file_inventory.json').write_text(json.dumps(inv,ensure_ascii=False,indent=2),encoding='utf-8')
print('SUPPLEMENT',len(inv),'files',sum(x['bytes'] for x in inv))