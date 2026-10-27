from pathlib import Path
import pandas as pd,sys
O=Path(sys.argv[1])
for n in ['01_CD4_death/disease_program_comparisons.tsv','02_mito_autophagy/disease_program_comparisons.tsv']:
 a=pd.read_csv(O/n,sep='\t');print(n);print(a[a['filter']=='All'][['scope','module','difference','P','BH','loo_same_sign']].to_string(index=False))
a=pd.read_csv(O/'03A_state_pseudotime/program_trajectory_tests.tsv',sep='\t');print('PROGRAM_TRAJECTORY');print(a[a.run=='K15_R1'][['state','feature','mean_slope','P','BH']].to_string(index=False))