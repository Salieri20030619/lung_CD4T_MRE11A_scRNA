from pathlib import Path
import sys,pandas as pd
O=Path(sys.argv[1]);a=pd.read_csv(O/'04_chains/all_sender_receiver_associations.tsv.gz',sep='\t')
b=a[(a.sender=='AMlike')&(a.x.isin(['IL15','IL18','TNF']))&(a.y.str.contains('IKZF2|STAT1|Ikzf2|Stat1',regex=True))&(a.model=='Pooled')]
print(b[['scope','x','y','rho','P','BH','loo_same_sign']].to_string(index=False))
f=O/'02_sender/ligand_activity_all_states.tsv'
if f.exists():
 a=pd.read_csv(f,sep='\t');print('TOP_LIGANDS');print(a[a['rank']<=3][['scope','direction','test_ligand','rank','aupr_corrected']].to_string(index=False))
for s in ['AMlike','AllMac']:
 a=pd.read_csv(O/'02_sender'/f'{s}_DE_all_genes.tsv.gz',sep='\t');print('SOURCE_DE',s,list(a.columns));print(a[a.canonical.isin(['IL15','IL18','TNF'])].to_string(index=False))