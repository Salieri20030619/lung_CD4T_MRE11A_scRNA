from pathlib import Path
import pandas as pd,json
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/four_axis_screen_20261017_v15');assert (O/'provenance/computation_complete.json').exists()
a=pd.read_csv(O/'05_intersection/all_gene_state_evidence.tsv.gz',sep='\t');b=a[a.data_primary].copy();tier='primary_data'
if b.empty:
 b=a[a.data_exploratory].copy();tier='exploratory_only'
 if b.empty:b=a[a.quality_pass&a.DE_pass&a.scope.ne('Core')].copy();tier='nearest_only'
 b['absFC']=b.logFC.abs();b=b.sort_values(['damage_pooled','mac_pooled','LOO_direction_agree','absFC','gene'],ascending=[False,False,False,False,True]).drop_duplicates('gene').head(10)
else:b=b.sort_values(['damage_primary','mac_primary','LOO_direction_agree','logFC','gene'],ascending=[False,False,False,False,True])
b['literature_tier']=tier;b.to_csv(O/'06_novelty/literature_worklist.tsv',sep='\t',index=False)
print('STAGE_COUNTS');print(pd.read_csv(O/'05_intersection/stage_counts.tsv',sep='\t').to_string(index=False));print('WORKLIST');print(b[['gene','scope','logFC','PValue','FDR','damage_primary','mac_primary','damage_pooled','mac_pooled','literature_tier']].to_string(index=False))
print('PRIOR_CANDIDATES');print(a[a.gene.isin(['IKZF2','STAT1','CXCL16','IL1R2','PENK','TNFRSF9','LAG3','TTC39C','FOXP1'])][['gene','scope','DE_pass','quality_pass','damage_primary','mac_primary','data_primary','missing_axes']].to_string(index=False))
for file in ['shortlist_damage_links.tsv','shortlist_macrophage_links.tsv']:
 f=O/'05_intersection'/file;d=pd.read_csv(f,sep='\t');d=d[d.gene.isin(b.gene)]
 if 'sender' in d:
  d=d[d.primary_pass];cols=['gene','scope','sender','ligand','P','BH','rho','adj_P','adj_BH','detected_receptors','DE_receptors','ligand_DE_P']
 else:d=d[d.primary_pass];cols=['gene','scope','module','rho','P','BH','adj_P','adj_BH']
 print(file);print(d[cols].head(80).to_string(index=False))