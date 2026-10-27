from pathlib import Path
import zipfile,hashlib,pandas as pd
W=Path('C:/Users/31937/lung_four_axes_20261017_v15');f=W/'v15_computational_tables.zip';assert hashlib.sha256(f.read_bytes()).hexdigest()=='da918da250c838c8292632b8ec87c4d15804a2c003d6291234072830f5b72cc1'
with zipfile.ZipFile(f) as z:assert z.testzip() is None;z.extractall(W/'tables')
a=pd.read_csv(W/'tables/05_intersection/exploratory_data_candidates.tsv',sep='\t');print('EXPLORATORY',len(a),'UNIQUE',a.gene.nunique());print(a[['gene','scope','logFC','PValue','FDR','damage_primary','mac_primary','damage_pooled','mac_pooled']].to_string(index=False))