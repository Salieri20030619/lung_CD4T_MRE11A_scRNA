from pathlib import Path
import pandas as pd
H=Path(__file__).parent.parent;a=pd.read_csv(H/'lung_three_cohort_20261017_v16/data/tables/sample_inventory.tsv',sep='\t');z=a[a.scope=='Core'];print(z.to_string(index=False));print(a[(a.dataset=='GSE183682')&a.scope.isin(['Naive','ActCD4','Treg'])][['sample','scope','group','day','cells']].to_string(index=False))