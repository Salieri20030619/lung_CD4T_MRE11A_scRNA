from pathlib import Path
import pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');V=P/'11_reports/integrated_macrophage_CD4_20261015_v4';M=P/'02_analysis/mechanism_step1_20260927_152218_s4V5eF';g=pd.read_csv(M/'comm_input/CD4_genes.tsv',sep='\t');print('GENE_MAP',g.head().to_dict('records'))
for sc in ['Naive','ActCD4','Treg','Core']:
 d=pd.read_csv(V/f'03_receiver/{sc}_DE_all_genes.tsv.gz',sep='\t');a=d[~d.canonical.isin(g.canonical)];print(sc,'missing',len(a));print(a[['gene_id','symbol','canonical']].head(20).to_string(index=False))