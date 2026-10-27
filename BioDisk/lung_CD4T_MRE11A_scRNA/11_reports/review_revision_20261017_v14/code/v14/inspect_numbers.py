from pathlib import Path
import pandas as pd
P=Path('C:/Users/31937/lung_review_20261017_v14/baseline_v13/supplement')
for name in ['S09_STAT1_gene_associations.tsv','S10_fixed_composition_STAT1_genes.tsv','S19_NicheNet_all_source_rankings.tsv']:
 a=pd.read_csv(P/name,sep='\t');print(name,list(a))
 if 'y' in a:a=a[a.y=='GSDMD']
 if 'seed' in a:a=a[a.seed==17]
 if 'scope' in a:a=a[a.scope.isin(['Core','ActCD4'])]
 if 'test_ligand' in a:a=a[a.test_ligand=='CXCL16']
 if 'ligand' in a:a=a[a.ligand=='CXCL16']
 print(a.head(12).to_string(index=False))