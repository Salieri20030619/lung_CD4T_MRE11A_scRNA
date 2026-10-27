from pathlib import Path
import pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');B=P/'11_reports';S=B/'independent_candidates_revision_20261005_v1/results/source_evidence'
d=pd.read_csv(next(S.glob('old_TF_program.tsv*')),sep='\t')
q=d[(d.tf=='STAT1')&(d['filter']=='All')&(d.variant=='DropOverlap')&(d.scope.str.contains('GZMA'))&(d.model=='SourceDep')&(d.feature=='score_Inflam')]
print('EXACT_GZMA',q.to_string(index=False))
T=B/'technical_delivery_20261002_ikzf2_stat1_v1/revision_v2/final_release'
for f in (T/'code').glob('*.py'):
 txt=f.read_text();ls=txt.splitlines()
 for i,l in enumerate(ls):
  if 'S2_TF_associations' in l or ('FDR' in l and ('groupby' in l or 'bh' in l.lower())): print(f.name,i+1,'\n'.join(ls[max(0,i-2):i+4]))
for root in [P/'02_analysis',P/'10_tables/CD4_final/pathway_v1']:
 for f in root.rglob('*'):
  if f.suffix.lower() in ['.py','.r'] and f.stat().st_size<500000:
   txt=f.read_text(errors='replace');ls=txt.splitlines()
   if 'source_adjusted_full_DEG' in txt or ('S2_TF_associations' in txt and 'groupby' in txt):
    print('SOURCE',f)
    for i,l in enumerate(ls):
     if any(w in l for w in ['source_adjusted_full_DEG','groupby','p.adjust','ncol','design','sample','FDR']):print(i+1,l[:250])
