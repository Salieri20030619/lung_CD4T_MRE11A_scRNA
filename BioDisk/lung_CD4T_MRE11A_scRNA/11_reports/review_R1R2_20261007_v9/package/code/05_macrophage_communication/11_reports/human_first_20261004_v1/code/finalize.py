from pathlib import Path
import hashlib,json,zipfile,pandas as pd,datetime,platform,sys
from PIL import Image
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/human_first_20261004_v1';A=P/'02_analysis';R=O/'results';D=O/'deliverables';T=P/'11_reports/technical_delivery_20261002_ikzf2_stat1_v1/revision_v2/final_release'
inputs=[A/'mechanism_step2_20260927_154919_EsN4gk/input'/('human_'+s) for s in ['cells.tsv','genes.tsv','counts.mtx']]
inputs += [A/'scenic3_20261001_123452_Kvp9h5/input'/('human_mac_'+s) for s in ['cells.tsv','genes.tsv','counts.mtx']]
inputs += [T/'prior_review_repaired/scenic3_review/results/reviewed_cells.tsv.gz',T/'source_tables/M2_human_sample_scores.tsv',T/'data/M2_all_module_group_models.tsv',A/'mechanism_step3_20260928_103057_lIRL38/tables/human_source_models_BH.tsv',P/'11_reports/front_half_20261004_v1/results/group_context_models.tsv',P/'11_reports/missing_links_20261003_v1/results/sender_all_samples.tsv']
inputs += list((A/'mechanism_step3_20260928_103057_lIRL38/tables').glob('DE_CD4_*.tsv'))
inputs += [A/'mechanism_step1_20260927_152218_s4V5eF/comm_input'/('CD4_'+s) for s in ['cells.tsv','genes.tsv','counts.mtx.gz']]
inputs += [A/'mechanism_step1_20260927_152218_s4V5eF/tables/group_effects.tsv',A/'mechanism_step5_20260928_153825_Q4sB1x/tables/geneset_members.tsv',A/'mechanism_step5_20260928_153825_Q4sB1x/tables/geneset_status.tsv',A/'mechanism_step5_20260928_153825_Q4sB1x/tables/ligand_activities.tsv']
inputs += [Path(v['path']) for v in json.loads((O/'provenance/human_prior_manifest.json').read_text())]
inputs += [P/'00_raw/nichenet_v2_mouse/ligand_target_matrix_nsga2r_final_mouse.rds',P/'00_raw/nichenet_v2_mouse/lr_network_mouse_21122021.rds']
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
manifest=[dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(set(inputs))];pd.DataFrame(manifest).to_csv(O/'provenance/source_manifest_complete.tsv',sep='\t',index=False)
checks=json.loads((O/'provenance/validation.json').read_text());images=[]
for f in sorted((O/'figures').glob('*.png')):
 with Image.open(f) as im:im.verify()
 with Image.open(f) as im:images.append(dict(file=f.name,width=im.width,height=im.height))
checks['figure_structural_QA']=images;checks['visual_QA_status']='PENDING_EXPLICIT_TRANSFER_APPROVAL; no local visual inspection claimed'
for f in sorted(R.glob('stage2_LT_*.tsv.gz')):
 d=pd.read_csv(f,sep='\t');assert not d.duplicated(['ligand','target']).any();assert d.weight.notna().all() and (d.weight>=0).all()
checks['all_exported_ligand_target_pairs_unique_finite_nonnegative']=True
html=(D/'human_first_evidence_report_zh.html').read_text();assert html.count('<figure>')==3 and html.count('</figure>')==3 and '231,168' in html
checks['report_embeds_three_figures']=True
checks['completed_time_UTC']=datetime.datetime.now(datetime.timezone.utc).isoformat();checks['python']=sys.version
(O/'provenance/final_QA.json').write_text(json.dumps(checks,indent=2))
(D/'REPRODUCTION_zh.txt').write_text('本包提供完整代码及派生统计结果。原始矩阵和RDS对象保留服务器，不包含于ZIP。\n使用同一服务器环境从原始输入复现：\n/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python code/reproduce.py /mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/human_first_repro_NEWNAME\n新目录必须不存在。代码核验资源SHA256，禁止覆盖本轮结果。\n输入路径和SHA256见provenance/source_manifest_complete.tsv。源计数与临床扩展资料不能由本ZIP自行替代。\n人类经典prior=Zenodo3260758；鼠prior=Zenodo7074291；不要比较两版本rank的数值等价性。\nstatistical_families_and_lineage.tsv严格区分旧未计算FDR、修复表和当前BH家族。\n',encoding='utf-8')
# Do not package raw single-cell rosters, count matrices, aggregate RDS, old source archive, or interrupted download.
files=[]
for sub in ['code','results','provenance','figures']:
 for f in sorted((O/sub).iterdir()):
  if not f.is_file():continue
  if f.name in ['stage1_cell_dictionary.tsv.gz','stage1_sender_cell_dictionary.tsv.gz','stage2_aggregates.rds','context.zip','DELIVERY_STATUS.json']:continue
  if f.suffix in ['.lock','.pid'] or f.name in ['fetch_human_prior.py','use_cached_prior.py']:continue
  files.append(f)
files+=sorted(f for f in D.iterdir() if f.suffix in ['.html','.txt'])
outmanifest=[dict(file=str(f.relative_to(O)),bytes=f.stat().st_size,sha256=sha(f)) for f in files];pd.DataFrame(outmanifest).to_csv(D/'OUTPUT_MANIFEST.tsv',sep='\t',index=False);files.append(D/'OUTPUT_MANIFEST.tsv')
zp=D/'human_first_reproducible_20261004.zip'
with zipfile.ZipFile(zp,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for f in files:z.write(f,str(f.relative_to(O)))
with zipfile.ZipFile(zp) as z:assert z.testzip() is None;assert not any(n.endswith('.rds') or 'cell_dictionary' in n or n.endswith('.mtx') for n in z.namelist())
receipt=dict(report=str(D/'human_first_evidence_report_zh.html'),report_sha256=sha(D/'human_first_evidence_report_zh.html'),zip=str(zp),zip_bytes=zp.stat().st_size,zip_sha256=sha(zp),entries=len(files),server_complete=True,library_status='PENDING_TRANSFER_APPROVAL')
(O/'provenance/DELIVERY_STATUS.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt,indent=2),flush=True)
