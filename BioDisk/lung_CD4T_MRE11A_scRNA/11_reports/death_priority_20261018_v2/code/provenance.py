from pathlib import Path
import json,hashlib,importlib.metadata as md,sys
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/death_priority_20261018_v2');P=O.parent.parent;B=P/'11_reports/three_cohort_four_axis_20261017_v16'
paths=list((P/'00_raw/GSE183682').glob('*matrix.mtx.gz'))+list((P/'00_raw/GSE183682').glob('*features.tsv.gz'))+list((P/'00_raw/GSE183682').glob('*barcodes.tsv.gz'))
paths += [B/'resources/CellChat_protein_LR_expanded.tsv',B/'inputs/edgeR_fits.rds',P/'00_raw/nichenet_v2_mouse/ligand_target_matrix_nsga2r_final_mouse.rds',P/'02_analysis/scenic1_20260929_103513_xRoEnr/Mouse/All/seed17/regulons.csv',P/'02_analysis/GSE314943/stage3_20260905_214228_VPtd0r/objects/whole_annotated.h5ad']
paths += list((P/'02_analysis/mechanism_step1_20260927_152218_s4V5eF/comm_input').glob('Mac_*'))
rows=[]
for p in paths:
 if not p.is_file():continue
 h=hashlib.sha256()
 with p.open('rb') as f:
  for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
 rows.append(dict(path=str(p),bytes=p.stat().st_size,sha256=h.hexdigest()))
(O/'inputs/additional_provenance.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
versions={k:md.version(k) for k in ['numpy','pandas','scipy','scanpy','anndata','harmonypy','matplotlib']}
(O/'logs/python_environment.json').write_text(json.dumps({'python':sys.version,'packages':versions},indent=2),encoding='utf-8')
print('PROVENANCE',len(rows))
