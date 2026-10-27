from pathlib import Path
import pandas as pd,json,hashlib,zipfile,subprocess
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/human_alignment_20261005_v1';R=O/'results';H=P/'11_reports/human_first_20261004_v1'
rd=lambda n:pd.read_csv(R/n,sep='\t',low_memory=False)
assert len(rd('all_fixed_candidates_full_context.tsv.gz'))==4032
assert len(rd('fixed_TF_bridge_associations.tsv.gz'))==2604
assert len(rd('previous_nine_associations_reproduced.tsv'))==9
assert len(rd('TF_metadata_mismatch_audit.tsv.gz'))==0
assert len(rd('fixed_candidates_mouse_direction_context.tsv.gz'))>0
assert rd('ligand_disease_target_TF_edge_alignment.tsv.gz').shape[0]>0
assert not (rd('fixed_TF_bridge_associations.tsv.gz').query("link == 'ligand_to_TF'").FDR<.05).any()
p=rd('stage1_all_finite_models.tsv');print('PHENOTYPE_TREG_ROWS',len(p[(p.scope.str.contains('GZMA'))&(p.metric=='score_TReg')&(p['filter']=='All')]))
old=O/'provenance/compute.log'
if old.exists():old.rename(O/'provenance/INVALID_initial_draft_compute.log')
for f in (H/'provenance').glob('*'):
 if f.is_file() and any(t in f.name.lower() for t in ['prior','resource','manifest','spec']):
  dest=O/'provenance'/('reused_'+f.name)
  if f.stat().st_size<1000000:dest.write_bytes(f.read_bytes())
readme='''本轮使用服务器原有环境与对象，不联网安装、不改变配置。\n执行顺序：\n1 Rscript code/align_export.R\n2 Rscript code/align_prior.R\n3 python code/align_compute.py\n4 python code/align_supplement.py\n5 python code/align_roster.py\n6 python code/align_report.py\n7 python code/align_package.py\nR=/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript\nPython=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python\nOPENBLAS_NUM_THREADS=1; OMP_NUM_THREADS=1。脚本使用绝对源路径；原始数据不随此包复制。\n所有结果以最终results及align_compute.log为准。INVALID_initial_draft_compute.log仅保留调试审计，来源列尚未合并的初轮AUC统计无效。\nZIP仅服务器保留。用户仅授权HTML报告与三图转移和Library上传。\n'''
(O/'README_reproduce_zh.txt').write_text(readme,encoding='utf-8')
inputs=[P/'02_analysis/mechanism_step2_20260927_154919_EsN4gk/input/human_cells.tsv',P/'02_analysis/scenic3_20261001_123452_Kvp9h5/input/human_mac_cells.tsv',H/'results/stage2_aggregates.rds',H/'results/stage2_evidence_all_candidates.tsv.gz',H/'results/stage2_targetset_members.tsv.gz',H/'results/stage2_background_members.tsv.gz',P/'11_reports/technical_delivery_20261002_ikzf2_stat1_v1/revision_v2/final_release/data/S2_TF_sample_activity.tsv',P/'11_reports/technical_delivery_20261002_ikzf2_stat1_v1/revision_v2/final_release/data/TF_network_all_edges.tsv']
source=[]
for f in inputs:
 h=hashlib.sha256()
 with f.open('rb') as fh:
  for chunk in iter(lambda:fh.read(8*1024*1024),b''):h.update(chunk)
 source.append(dict(path=str(f),bytes=f.stat().st_size,sha256=h.hexdigest()))
(O/'provenance/source_manifest.json').write_text(json.dumps(source,indent=2))
files=[f for f in O.rglob('*') if f.is_file() and f.suffix!='.zip' and f.name not in ['file_manifest.json','package_receipt.json']]
manifest=[dict(path=str(f.relative_to(O)),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in files];(O/'provenance/file_manifest.json').write_text(json.dumps(manifest,indent=2))
archive=O/'deliverables/human_alignment_20261005_reproducible_SERVER_ONLY.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for f in files+[O/'provenance/file_manifest.json']:z.write(f,str(f.relative_to(O)))
with zipfile.ZipFile(archive) as z:assert z.testzip() is None
receipt={'archive':str(archive),'bytes':archive.stat().st_size,'sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),'files':len(files)+1,'status':'SERVER_ONLY_NO_TRANSFER','report':str(O/'deliverables/human_alignment_20261005_report_zh.html'),'mouse_rows':len(rd('fixed_candidates_mouse_direction_context.tsv.gz'))}
(O/'provenance/package_receipt.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt,indent=2))
