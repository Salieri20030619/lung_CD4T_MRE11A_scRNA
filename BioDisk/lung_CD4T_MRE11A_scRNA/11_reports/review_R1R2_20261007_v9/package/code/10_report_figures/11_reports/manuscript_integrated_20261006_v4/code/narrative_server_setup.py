from pathlib import Path
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/manuscript_integrated_20261006_v4'
for sub in ['code','deliverables','figures','provenance']:(O/sub).mkdir(parents=True,exist_ok=True)
print(str(O))
