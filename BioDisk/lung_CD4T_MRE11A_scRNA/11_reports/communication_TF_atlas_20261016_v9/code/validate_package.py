from pathlib import Path
import pandas as pd, numpy as np, json, hashlib, shutil, zipfile
import fitz
R=Path('C:/codex_test'); O=R/'output/silicosis_v9'; S=R/'tmp/silicosis_v9/v4'; T=R/'tmp/silicosis_v9'
checks={}
def same(plot,source,keys,values):
 a=pd.read_csv(O/'plot_data'/plot,sep='\t'); b=pd.read_csv(S/source,sep='\t')
 assert not a.duplicated(keys).any(),plot
 x=a.merge(b,on=keys,suffixes=('_plot','_source'),how='left',validate='1:1',indicator=True)
 assert x._merge.eq('both').all(),plot
 for v in values: assert np.allclose(x[v+'_plot'],x[v+'_source'],equal_nan=True), (plot,v)
 checks[plot]={'rows':len(a),'source_values_match':True,'checked_columns':values}
same('F01_LR_coverage.tsv','02_sender/all_prior_LR_expression_coverage.tsv.gz',['sender','scope','ligand','receptor'],['sender_case_detected','receiver_case_detected','sender_case_mean_pct','receiver_case_mean_pct'])
same('F02_ligand_prediction.tsv','02_sender/source_specific_activity_ranks.tsv',['sender','scope','direction','ligand'],['source_specific_rank','aupr_corrected'])
same('F02_receptor_joint_detection.tsv','03_receiver/receptor_complex_codetection.tsv',['lineage','scope','sample','ligand','rule','components'],['joint_pct','cells'])
same('F03_sample_AUC.tsv','03_receiver/TF_all_regulon_sample_activity.tsv.gz',['seed','scope','sample','regulon'],['mean_auc','cells'])
same('F03_disease_effects.tsv','03_receiver/all_TF_activity_group_effects.tsv',['seed','scope','regulon'],['delta_mean_auc','P','BH'])
same('F04_TF_program_associations.tsv','03_receiver/all_TF_program_associations.tsv.gz',['seed','scope','model','x','y'],['rho','P','BH'])
same('F06_sender_STAT1_bridge.tsv','04_chains/all_sender_receiver_associations.tsv.gz',['seed','scope','sender','model','x','y'],['rho','P','BH'])
a=pd.read_csv(O/'plot_data/F05_prior_module_edges.tsv',sep='\t'); b=pd.read_csv(S/'03_receiver/focal_TF_prior_targets_and_modules.tsv',sep='\t'); pd.testing.assert_frame_equal(a,b[b.module.notna()].reset_index(drop=True)); checks['F05_prior_module_edges']='exact subset'
a=pd.read_csv(O/'plot_data/F05_target_DE.tsv',sep='\t')
for scope,g in a.groupby('scope'):
 b=pd.read_csv(S/f'03_receiver/{scope}_DE_all_genes.tsv.gz',sep='\t').set_index('canonical')
 for _,r in g.iterrows(): assert np.allclose([r.logFC,r.P,r.BH],b.loc[r.gene,['logFC','PValue','FDR']].astype(float))
checks['F05_target_DE']='all values match original full-family DE tables'
d=fitz.open(O/'silicosis_communication_TF_v9.pdf'); old=fitz.open(O/'baseline/report_v8.pdf')
assert len(d)==len(old)+7
for i in range(len(old)):
 assert d[i+7].get_text()==old[i].get_text()
 assert d[i+7].get_pixmap().samples==old[i].get_pixmap().samples
checks['baseline_preservation']='all 19 pages text and rendered pixels identical'
assert all('\ufffd' not in p.get_text() for p in d)
for i in range(7):
 for block in d[i].get_text('blocks'): assert block[0]>=0 and block[1]>=0 and block[2]<=843 and block[3]<=596
checks['new_page_bounds']='all text inside page'; checks['replacement_glyphs']=0
for i in [2,6]: d[i].get_pixmap(matrix=fitz.Matrix(1.5,1.5)).save(T/f'final_{i+1:02d}.png')
fonts={f[1] for i in range(7) for f in d[i].get_fonts()}; checks['new_pages_font_extensions']=sorted(fonts)
src=json.loads((O/'provenance/source_sha256.json').read_text(encoding='utf-8'))
for rel,h in src.items():
 p=S/rel; assert hashlib.sha256(p.read_bytes()).hexdigest()==h
 dst=O/'source_tables'/rel; dst.parent.mkdir(parents=True,exist_ok=True); shutil.copy2(p,dst)
checks['source_files_verified']=len(src); checks['new_inferential_tests']=0
(O/'provenance/final_validation.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2),encoding='utf-8')
(O/'README.txt').write_text('矽肺通讯与TF图谱补充整合版 v9\n\n先打开 silicosis_communication_TF_v9.pdf：26页，前7页为新增图谱，后19页为v8完整报告。\nreport_v9.html：可放大查看新增图件，并链接v8原HTML。\nfigures：六组图的PNG、SVG及矢量PDF。\nplot_data：逐图取值，保留完整原始P、BH和模型字段。\nsource_tables：本轮使用的冻结输入表，与服务器v4结果逐文件SHA256一致。\nsupplement：v8全部16份补充表。\nprovenance：来源SHA256、绘图口径和验收记录。\ncode/build_v9.py：本轮构建脚本；顶部ROOT/SRC/OLD/OUT为本机执行路径，移机复现时需改为对应目录。\n\n本轮未新增推断性检验；TF展示子集不作为新候选定论。原v8各页文字及像素与合并PDF中的对应页一致。新增页已渲染检查；本轮不是对v8全部结论的独立重审。\n服务器项目目录：/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/\n',encoding='utf-8')
shutil.copy2(__file__,O/'code/validate_package.py')
rows=[]
for p in sorted(O.rglob('*')):
 if p.is_file() and p.name!='MANIFEST_SHA256.tsv': rows.append({'file':p.relative_to(O).as_posix(),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
pd.DataFrame(rows).to_csv(O/'MANIFEST_SHA256.tsv',sep='\t',index=False,encoding='utf-8')
zip_path=O.parent/'silicosis_communication_TF_v9.zip'
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED) as z:
 for p in sorted(O.rglob('*')):
  if p.is_file(): z.write(p,p.relative_to(O).as_posix())
with zipfile.ZipFile(zip_path) as z: assert z.testzip() is None
print(json.dumps({'status':'PASS','checks':checks,'archive_bytes':zip_path.stat().st_size,'archive_sha256':hashlib.sha256(zip_path.read_bytes()).hexdigest()},ensure_ascii=False,indent=2))
