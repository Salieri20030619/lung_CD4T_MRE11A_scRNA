from pathlib import Path
import pandas as pd
O=Path(__file__).parent;U=O/'server_tables';T=O/'tables'
for n in ['ORA_all.tsv','GSEA_all.tsv']:
 d=pd.read_csv(U/n,sep='\t');print(n,len(d));p='P' if n.startswith('ORA') else 'pval';print(d.sort_values(p).groupby(['scope','database']).head(2)[[c for c in ['scope','database','description','direction','k','n','NES',p,'FDR','padj'] if c in d]].to_string(index=False))
sc=pd.read_csv(O.parent/'lung_append_analysis_20261017_v18/tables/programme_group_effects.tsv',sep='\t');print(sc[(sc.cohort.astype(str)=='314')&(sc.variant=='Raw')&(sc.scope.isin(['ActCD4','Treg']))&sc.program.isin(['PyroExec','ApoExec','Inflam','MitoISR'])][['scope','program','effect','P','FDR']].to_string(index=False))
print(pd.read_csv(T/'candidate_detection_summary.tsv',sep='\t').query("dataset=='GSE314943' and variant=='BG1' and ((gene=='IL1R2' and scope=='Treg') or (gene=='SKIL' and scope=='Treg') or (gene=='TTC39C' and scope=='ActCD4'))").to_string(index=False))