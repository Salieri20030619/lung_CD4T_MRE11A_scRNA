from pathlib import Path
import pandas as pd
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/review_revision_20261017_v14')
for f in sorted((O/'supplement_v14').glob('*')): print(f.name,f.stat().st_size)
a=pd.read_csv(O/'supplement_v14/S21_BH_recalculation_audit.tsv',sep='\t');print('BH',len(a),a.status.value_counts().to_dict())
a=pd.read_csv(O/'supplement_v14/S25_regulon_module_control_overlap.tsv',sep='\t');print(a[a.seed==17][['tf','module','control_overlap']].to_string(index=False))