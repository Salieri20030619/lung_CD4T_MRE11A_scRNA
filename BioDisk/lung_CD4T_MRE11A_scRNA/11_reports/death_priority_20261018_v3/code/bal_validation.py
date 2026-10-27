from pathlib import Path
import pandas as pd,numpy as np
from scipy.io import mmread
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');I=P/'02_analysis/GSE174725/suppC2_20260911_164417_l173TC/input';O=P/'11_reports/death_priority_20261018_v3';V=P/'11_reports/death_priority_20261018_v2'
m=pd.read_csv(I/'cells.tsv.gz',sep='\t');g=pd.read_csv(I/'genes.tsv',sep='\t');print('GENE_COL',g.columns.tolist(),g.head(2).to_dict('records'));print(m.groupby(['group','sample','subtype']).size().to_string())
x=mmread(I/'counts.mtx').tocsr();print('SHAPE',x.shape,len(m),len(g))
