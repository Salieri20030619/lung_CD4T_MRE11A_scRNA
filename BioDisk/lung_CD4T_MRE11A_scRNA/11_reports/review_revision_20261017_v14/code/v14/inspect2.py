from pathlib import Path
import pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA'); V=P/'11_reports/integrated_macrophage_CD4_20261015_v4'; B=P/'11_reports/STAT1_pyro_bridge_20261015_v5'
for f in [V/'03_receiver/receptor_complex_codetection.tsv',V/'03_receiver/TF_regulon_coverage.tsv',V/'02_sender/lead_ligand_source_identity_audit.tsv',B/'tables/source_complex_identity_audit.tsv',B/'tables/strict_identity_sensitivity_roster.tsv',B/'tables/CD4_executor_identity_audit.tsv']:
 a=pd.read_csv(f,sep='\t');print('\nFILE',f.name,a.shape,list(a));print(a.head(3).to_string(index=False))
for f in (P/'02_analysis/scenic2_20261001_092655_ZKNWD6/objects').iterdir():print('SCENIC_OBJ',f.name)
import h5py
for rel in ['stage3_20260905_214228_VPtd0r/objects/whole_annotated.h5ad','stage4_20260905_215942_3XPis4/objects/TNK_harmony.h5ad']:
 with h5py.File(P/'02_analysis/GSE314943'/rel) as h:
  print('H5',rel,'obs',list(h['obs']),'var',list(h['var']),'layers',list(h.get('layers',{})))
print('PYTHON_OK')