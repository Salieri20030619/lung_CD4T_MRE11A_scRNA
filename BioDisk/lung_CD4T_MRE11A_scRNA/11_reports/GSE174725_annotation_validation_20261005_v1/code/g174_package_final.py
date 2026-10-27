from pathlib import Path
import json,hashlib,zipfile
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/GSE174725_annotation_validation_20261005_v1';D=O/'deliverables'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
refs=['11_reports/GSE174725_reaudit_20261005_v1/deliverables/SHA256_manifest.json','11_reports/human_first_20261004_v1/results/stage1_complete_program_members.tsv','11_reports/human_first_20261004_v1/results/stage2_targetset_members.tsv.gz','11_reports/human_first_20261004_v1/results/stage2_CD4_mapping.tsv','11_reports/human_gap_diagnostics_20261005_v1/results/independence_program_members.tsv.gz','11_reports/human_gap_diagnostics_20261005_v1/results/matched_control_members.tsv.gz','11_reports/human_alignment_20261005_v1/results/all_fixed_candidate_LR_pairs.tsv.gz','02_analysis/mechanism_step2_20260927_154919_EsN4gk/tables/human_score_controls.tsv.gz']
(O/'provenance/dependency_source_sha256.json').write_text(json.dumps([dict(path=str(P/p),bytes=(P/p).stat().st_size,sha256=sha(P/p)) for p in refs],indent=2))
(O/'provenance/annotation_scope_addendum.txt').write_text('Supplementary mature/non-AM marker review performed without additional clustering or outcome-based relabeling. Existing clusters0/6 support mature macrophage with alveolar-like state. No other coherent mature macrophage cluster supported by the reviewed C1Q/APOE/APOC1/MRC1/TREM2/GPNMB/LGMN program plus competing lineage evidence. Rare or transitional non-AM macrophages not biologically excluded. LYVE1 is not common-mapped, not an expression-negative finding. BALF_AM does not mean all macrophages in the organism. See mature_nonAM_macrophage_marker_audit.tsv.\n')
report=D/'GSE174725_annotation_validation_20261005_report_zh.html';files=[p for folder in ['code','results','provenance','figures'] for p in (O/folder).rglob('*') if p.is_file()]+[report,D/'README.txt'];man=[dict(path=str(p.relative_to(O)),bytes=p.stat().st_size,sha256=sha(p),in_archive=p.suffix!='.h5ad') for p in files];m=D/'SHA256_manifest.json';m.write_text(json.dumps(man,indent=2));z=D/'GSE174725_annotation_validation_reproducible_SERVER_ONLY.zip'
with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED) as out:
 for p in files+[m]:
  if p.suffix!='.h5ad':out.write(p,str(p.relative_to(O)))
with zipfile.ZipFile(z) as out:
 assert out.testzip() is None
 for row in man:
  if row['in_archive']:assert hashlib.sha256(out.read(row['path'])).hexdigest()==row['sha256']
x=json.loads((D/'completion.json').read_text());x.update(archive_sha256=sha(z),archive_bytes=z.stat().st_size,manifest_entries=len(man),archive_integrity='ALL_ARCHIVED_BYTES_MATCH_MANIFEST');(D/'completion.json').write_text(json.dumps(x,indent=2));print(json.dumps(x,indent=2))
