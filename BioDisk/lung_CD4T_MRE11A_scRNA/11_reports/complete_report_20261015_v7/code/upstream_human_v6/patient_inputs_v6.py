from pathlib import Path
import sys,json,pandas as pd
O=Path(sys.argv[1]);P=O.parents[1];J=P/'11_reports/GSE296682_patient_localization_20261005_v1';print('FILES')
for root in [P/'00_raw/GSE296682',J/'code',J/'provenance']:
 print(root)
 if root.exists():
  for p in root.rglob('*'):
   if p.is_file():print(str(p),p.stat().st_size)
print('GENES',pd.read_csv(J/'results/fixed_feature_mapping.tsv',sep='\t').gene.tolist())