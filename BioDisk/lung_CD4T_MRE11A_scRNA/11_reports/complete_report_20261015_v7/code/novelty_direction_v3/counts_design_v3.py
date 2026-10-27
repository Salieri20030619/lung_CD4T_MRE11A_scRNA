from pathlib import Path
import pandas as pd,sys,json
O=Path(sys.argv[1]);a=pd.read_csv(O/'literature/search_log.tsv',sep='\t');print(a[['gene','topic','hit_count','returned','status']].to_string(index=False))
for acc in ['GSE314943','GSM9529140','GSE183682','GSM5567502']:
 print('\nGEO',acc)
 for l in (O/'direction'/(acc+'_public_metadata.txt')).read_text(encoding='utf-8').splitlines():
  if any(t in l.lower() for t in ['!series_overall_design','!sample_characteristics','!sample_treatment_protocol','!series_pubmed_id']):print(l[:1800])
a=pd.read_csv(O/'literature/retrieved_article_metadata.tsv',sep='\t').fillna('');sel=a[a.title.str.contains('Standigm|AHNAK Loss|Cross-talk|IFNα8|Biomarker potential|Molecular cloning and expression of human ST6',case=False)]
print(sel[['gene','pmid','pmcid','doi','title']].to_string(index=False))