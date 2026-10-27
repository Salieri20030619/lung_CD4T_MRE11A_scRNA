from pathlib import Path
import importlib.util,json,pandas as pd,numpy as np
from scipy import sparse
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');B=P/'11_reports/three_cohort_four_axis_20261017_v16';O=P/'11_reports/death_priority_20261018_v31';O.mkdir(exist_ok=True)
for n in ['code','logs','tables','objects','figures','inputs']:(O/n).mkdir(exist_ok=True)
for mod in ['scanpy','anndata','harmonypy','igraph','leidenalg','sklearn','umap','statsmodels']:print('PACKAGE',mod,bool(importlib.util.find_spec(mod)),flush=True)
for ds,sc in [('GSE314943','Core'),('GSE183682','Core'),('GSE274372','TLike')]:
 stem=ds+'__'+sc;m=pd.read_csv(B/'inputs'/(stem+'.cells.tsv.gz'),sep='\t');x=sparse.load_npz(B/'inputs'/(stem+'.cellcounts.npz'));g=pd.read_csv(B/'inputs'/(stem+'.genes.txt'),header=None)[0];print('INPUT',ds,x.shape,len(g),m.columns.tolist(),flush=True);print(m.groupby(['sample','group','day','subtype']).size().to_string(),flush=True)
go=pd.read_csv(B/'resources/GO_mouse_membership.tsv',sep='\t');print('GO_FIELDS',go.columns.tolist(),flush=True);print(go.head(1).to_string(index=False),flush=True)