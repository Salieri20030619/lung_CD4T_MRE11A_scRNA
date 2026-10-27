from pathlib import Path
import pandas as pd,json,zipfile,hashlib,shutil
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/four_axis_screen_20261017_v15';R=O/'report/review_tables';R.mkdir(parents=True,exist_ok=True)
for d in ['01_quality','02_DE','03_damage','04_macrophage','05_intersection','06_novelty','provenance']:
 for f in (O/d).glob('*'):
  if f.is_file() and f.suffix in ['.tsv','.gz','.json','.log','.txt']:
   # retain full inference tables; exclude single-cell diagnosis data in reviewer bundle
   if 'cell_' in f.name or f.name.startswith('decontx_'):continue
   dest=R/d/f.name;dest.parent.mkdir(exist_ok=True,parents=True);shutil.copy2(f,dest)
shutil.copy2(O/'PLAN_FROZEN.json',R/'PLAN_FROZEN.json')
a=pd.read_csv(O/'05_intersection/all_gene_state_evidence.tsv.gz',sep='\t');old=a[a.gene.isin(['IKZF2','STAT1','CXCL16','IL1R2','PENK','TNFRSF9','LAG3','TTC39C','FOXP1'])];old.to_csv(O/'05_intersection/prior_candidates_audit.tsv',sep='\t',index=False);shutil.copy2(O/'05_intersection/prior_candidates_audit.tsv',R/'05_intersection/prior_candidates_audit.tsv')
counts={'gene_state_tests_primary':int(a.scope.ne('Core').sum()),'DE_primary':int((a.scope.ne('Core')&a.DE_pass).sum()),'Q_DE_primary':int((a.scope.ne('Core')&a.DE_pass&a.quality_pass).sum()),'data_primary':int(a.data_primary.sum()),'data_exploratory':int(a.data_exploratory.sum()),'damage_test_rows':len(pd.read_csv(O/'03_damage/all_gene_program_tests.tsv.gz',sep='\t',usecols=['gene'])),'mac_test_rows':sum(len(pd.read_csv(f,sep='\t',usecols=['gene'])) for f in (O/'04_macrophage').glob('*_all_links.tsv.gz') if f.name!='consolidated_all_links.tsv.gz')};(O/'provenance/screen_counts.json').write_text(json.dumps(counts,indent=2),encoding='utf-8');shutil.copy2(O/'provenance/screen_counts.json',R/'provenance/screen_counts.json')
zpath=O/'report/v15_computational_tables.zip'
with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED) as z:
 for f in R.rglob('*'):
  if f.is_file():z.write(f,f.relative_to(R).as_posix())
print(json.dumps(counts));print('BUNDLE',zpath.stat().st_size,hashlib.sha256(zpath.read_bytes()).hexdigest())