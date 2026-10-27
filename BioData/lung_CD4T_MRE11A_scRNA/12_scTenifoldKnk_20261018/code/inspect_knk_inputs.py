from pathlib import Path
import pandas as pd
from scipy import sparse
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/three_cohort_four_axis_20261017_v16/inputs')
for ds in ['GSE314943','GSE183682']:
 m=pd.read_csv(P/(ds+'__Core.cells.tsv.gz'),sep='\t');g=pd.read_csv(P/(ds+'__Core.genes.txt'),header=None)[0];x=sparse.load_npz(P/(ds+'__Core.cellcounts.npz'))
 print(ds,x.shape,m.columns.tolist(),flush=True)
 print(m.groupby(['sample','subtype','group']).size().to_string(),flush=True)
 print(g.head().to_list(), x.data[:20],flush=True)