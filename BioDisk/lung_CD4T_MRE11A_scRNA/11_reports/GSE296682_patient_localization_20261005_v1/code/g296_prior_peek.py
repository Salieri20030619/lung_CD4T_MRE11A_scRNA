from pathlib import Path
import pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports')
for p in [P/'human_first_20261004_v1/results/stage2_LR_expression_all_donors.tsv.gz',P/'human_alignment_20261005_v1/results/IL15_matched_coverage_summary.tsv']:
 d=pd.read_csv(p,sep='\t');print(p.name,d.shape,d.columns.tolist());print(d.head(4).to_string(index=False))
