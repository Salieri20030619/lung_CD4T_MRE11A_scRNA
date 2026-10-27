from pathlib import Path
import pandas as pd,sys
O=Path(sys.argv[1]);P=O.parents[1]
r=pd.read_csv(O/'03A_state_pseudotime/genomewide_trajectory_tests.tsv.gz',sep='\t')
for state,gs in [('Naive',['AHNAK']),('ActCD4',['MDFIC','TTC39C','ST6GALNAC3'])]:print(r[(r.state==state)&r.feature.isin(gs)][['state','run','feature','mean_slope','P','same_sign_fraction']].to_string(index=False))
H=P/'11_reports/GSE174725_annotation_validation_20261005_v1/results'
print((H/'frozen_scope_donor_counts.tsv').read_text(encoding='utf-8'))
a=pd.read_csv(H/'frozen_annotation_cells_SERVER_ONLY.tsv.gz',sep='\t',nrows=2);print('FROZEN_COLS',list(a.columns));print(a.to_string(index=False)[:1500])
for f in (P/'11_reports/GSE174725_annotation_validation_20261005_v1/code').glob('*.py'):
 lines=f.read_text(encoding='utf-8').splitlines()
 hits=[(i,l) for i,l in enumerate(lines) if 'CD4' in l and any(s in l for s in ['scope','frozen','old_CD4'])]
 if hits:print(f.name,'\n'.join(str(i+1)+':'+l for i,l in hits[:12]))