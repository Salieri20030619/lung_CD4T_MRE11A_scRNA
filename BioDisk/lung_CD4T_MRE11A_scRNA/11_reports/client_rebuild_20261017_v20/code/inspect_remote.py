from pathlib import Path
import anndata as ad,json
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/client_rebuild_20261017_v20'
for sub in ['stage3_20260905_214228_VPtd0r/objects/whole_annotated.h5ad','stage5_20260905_222457_a2dApT/objects/CD4_core.h5ad','stage9_20260906_080538_wN3sjO/objects/K15_trajectory.h5ad']:
 a=ad.read_h5ad(P/'02_analysis/GSE314943'/sub,backed='r');print(sub,a.shape,list(a.obs.columns),list(a.obsm.keys()));print(a.obs.head(1).to_string());print(a.var.head(2).to_string());a.file.close()
for root in [P/'02_metadata',P/'02_analysis/GSE314943']:
 files=[str(p) for p in root.rglob('*') if p.is_file() and any(k in p.name.lower() for k in ['kegg','gmt','enrich'])];print('RESOURCES',files[:60])