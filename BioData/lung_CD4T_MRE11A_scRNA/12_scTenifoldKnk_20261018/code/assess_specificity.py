from pathlib import Path
import json
import pandas as pd
from scipy.stats import spearmanr
O=Path('/mnt/BioData/bio/wangmingyuan/lung_CD4T_MRE11A_scRNA/12_scTenifoldKnk_20261018')
marks=list((O/'runs').glob('*/complete.json'));print('COMPLETE',len(marks));print('ERRORS',[str(f) for f in marks if json.loads(f.read_text())['status']!='success'])
for state in ['Saline','Silica']:
 d=O/'runs'/('All_Naive_'+state+'_s17');a=pd.read_csv(d/'HIF1A.tsv',sep='\t').set_index('gene');b=pd.read_csv(d/'ARL6IP1.tsv',sep='\t').set_index('gene');ii=a.index.intersection(b.index).difference(['HIF1A','ARL6IP1']);print(state,'HIF1A vs ARL6IP1 response rho',spearmanr(a.loc[ii,'distance'],b.loc[ii,'distance']).statistic)