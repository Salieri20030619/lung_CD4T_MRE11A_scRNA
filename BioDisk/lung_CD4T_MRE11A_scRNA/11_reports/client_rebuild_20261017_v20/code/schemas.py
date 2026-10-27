from pathlib import Path
import pandas as pd
H=Path(r'C:\Users\31937')
for rel in ['lung_append_analysis_20261017_v18/tables/programme_group_effects.tsv','lung_append_analysis_20261017_v18/tables/clean_TF_sample_AUC.tsv','lung_append_analysis_20261017_v18/tables/focal_TF_disease_effects.tsv','lung_three_cohort_20261017_v16/data/tables/candidate_complete_routes.tsv','lung_three_cohort_20261017_v16/data/tables/expressed_LR_with_decoy_literature.tsv']:
 d=pd.read_csv(H/rel,sep='\t');print(rel, d.shape, list(d.columns));print(d.head(1).to_string(index=False))