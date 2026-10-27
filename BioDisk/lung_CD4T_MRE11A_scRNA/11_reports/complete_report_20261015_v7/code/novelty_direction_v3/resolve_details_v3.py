from pathlib import Path
import sys,json,re,xml.etree.ElementTree as ET,pandas as pd
O=Path(sys.argv[1])
seen=False
for f in (O/'literature/cache').glob('MDFIC*.json'):
 for r in json.loads(f.read_text(encoding='utf-8')).get('resultList',{}).get('result',[]):
  if 'IFN' in r.get('title','') and 'memory CD4' in r.get('title','') and not seen:
   print('PREPRINT',r.get('doi'),r.get('firstPublicationDate'),r.get('pubTypeList'));print(re.sub('<[^>]+>',' ',r.get('abstractText',''))[:3800]);seen=True
f=O/'literature/cache/PMC7290873.xml'
root=ET.fromstring(f.read_bytes())
for n in root.iter('table-wrap'):
 if 'ST6GALNAC3' in ''.join(n.itertext()):print('IPF_TABLE_HEADING',' '.join(n.itertext())[:700])
a=pd.read_csv(O/'direction/sample_expression_estimands.tsv',sep='\t');print('EXTERNAL_COUNTS');print(a[(a.dataset=='GSE183682')&(a.scope=='ActCD4')&(a['filter']=='All')&(a.gene.isin(['MDFIC','TTC39C']))][['day','sample','group','gene','cells','positive_cells','target_UMI','full_library_UMI']].to_string(index=False))
a=pd.read_csv(O/'direction/direction_by_metric_and_filter.tsv',sep='\t');print('BG1');print(a[(a.scope=='ActCD4')&(a['filter']=='BG1')&(a.gene.isin(['MDFIC','TTC39C']))&(a.metric=='mean_CP10k')][['dataset','day','gene','difference']].to_string(index=False))
a=pd.read_csv(O/'direction/human_Core_composition_decomposition.tsv',sep='\t');print('HUMAN AHNAK PARTS');print(a[(a['filter']=='ReviewedAll')&(a.gene=='AHNAK')].to_string(index=False))