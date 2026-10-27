from pathlib import Path
import anndata as ad,pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
a=ad.read_h5ad(P/'11_reports/GSE174725_annotation_validation_20261005_v1/results/broad_annotation_review_SERVER_ONLY.h5ad',backed='r')
print('OBJ',a.shape,list(a.layers.keys()),list(a.obs.columns))
for c in a.obs.columns:
 if any(v in c.lower() for v in ['anno','broad','type','label','sample','group']):print(c,a.obs[c].value_counts().head(25).to_dict())
print('OBS',a.obs.head(2).to_dict('records'))
I=P/'02_analysis/GSE174725/suppC2_20260911_164417_l173TC/input'
m=pd.read_csv(I/'cells.tsv.gz',sep='\t');print('T_COVERAGE',m.groupby(['sample','group','subtype']).size().to_string())
