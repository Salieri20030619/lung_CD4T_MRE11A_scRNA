from pathlib import Path
import pandas as pd,json,re,sys
O=Path(sys.argv[1]);a=pd.read_csv(O/'literature/search_log.tsv',sep='\t');print(a[['gene','topic','hit_count','returned','status']].to_string(index=False))
for acc in ['GSE314943','GSM9529140','GSE183682','GSM5567502','GSE274372']:
 print('GEO',acc)
 s=(O/'direction'/(acc+'_public_metadata.txt')).read_text(encoding='utf-8')
 for line in s.splitlines():
  if any(t in line.lower() for t in ['!series_overall_design','!sample_title','!sample_characteristics','!sample_treatment_protocol','!sample_growth_protocol','!series_pubmed_id']):print(line[:2500])
ids=['38349059','41526926','35235341','34461439','31188487','36303158','10671520']
seen=set()
for f in (O/'literature/cache').glob('*.json'):
 for r in json.loads(f.read_text(encoding='utf-8')).get('resultList',{}).get('result',[]):
  if (r.get('pmid') in ids or ('IFN' in r.get('title','') and 'memory' in r.get('title',''))) and r.get('id') not in seen:
   seen.add(r.get('id'));print('\nPAPER',r.get('pmid'),r.get('doi'),r.get('title'));print(re.sub('<[^>]+>',' ',r.get('abstractText',''))[:4200])