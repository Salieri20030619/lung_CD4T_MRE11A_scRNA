from pathlib import Path
import pandas as pd
W=Path('C:/Users/31937/lung_review_20261017_v14/revision_v14');Q=W/'supplement_v14'
for f in ['S22_fixed_composition_gene_inputs.tsv','S22_Core_fixed_composition_TF_inputs.tsv','S22_sample_scope_roster.tsv']:
 a=pd.read_csv(Q/f,sep='\t');print(f,list(a));print(a.head(1).to_dict('records'))
a=pd.read_csv(W/'supplement/S20_GO_autophagy_measured_members.tsv',sep='\t');print('GO',list(a))
a=pd.read_csv(Q/'S23_targeted_proximal_associations.tsv',sep='\t');print('n_valid',a.P.notna().sum(),'n_invalid',a.P.isna().sum())