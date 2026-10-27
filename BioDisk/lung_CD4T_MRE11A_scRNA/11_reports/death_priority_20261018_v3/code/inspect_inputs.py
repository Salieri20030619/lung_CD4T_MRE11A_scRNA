from pathlib import Path
import pandas as pd,json
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
m=pd.read_csv(P/'02_analysis/mechanism_step2_20260927_154919_EsN4gk/input/human_cells.tsv',sep='\t')
print(m.groupby(['series','group','subtype'])['sample'].nunique().to_string())
print(m.groupby(['sample','series','group','subtype']).size().unstack(fill_value=0).to_string())
for p in (P/'02_analysis/GSE174725/suppC2_20260911_164417_l173TC').iterdir():print(p.name)
