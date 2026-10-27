from pathlib import Path
import pandas as pd,sys
O=Path(sys.argv[1]);a=pd.read_csv(O/'inputs/TF_target_prior.tsv',sep='\t');print('PRIOR',a.shape,list(a.columns),a.head(1).to_dict('records'))
a=pd.read_csv(O/'04_chains/all_sender_receiver_associations.tsv.gz',sep='\t');tf=pd.read_csv(O/'03_receiver/IKZF2_STAT1_primary_and_sensitivity.tsv',sep='\t');p=tf[tf.P<.05]
for r in p.itertuples():
 b=a[(a.scope==r.scope)&(a.sender=='AMlike')&(a.y=='AUC::'+r.x)&(a.model=='Pooled')&(a.P<.05)]
 print('LINK',r.scope,r.x,r.y,'SOURCE_HITS',len(b));print(b[['x','rho','P','BH']].head(12).to_string(index=False))
co=pd.read_csv(O/'03_receiver/receptor_complex_codetection.tsv',sep='\t');print('COMPLEX CASE');print(co[(co.group=='Silica')&(co.scope.isin(['ActCD4','Treg','AMlike']))&(co.ligand.isin(['IL15','IL18','TNF','IL15_transpresentation']))][['ligand','scope','sample','cells','joint_pct']].to_string(index=False))