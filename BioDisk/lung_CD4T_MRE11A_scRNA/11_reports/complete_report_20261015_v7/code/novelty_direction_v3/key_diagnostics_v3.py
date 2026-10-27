from pathlib import Path
import json,sys,pandas as pd
O=Path(sys.argv[1]);r=json.loads((O/'literature/cache/PMID_41862476.json').read_text(encoding='utf-8'));print('DISCOVERY_PAPER',r.get('pmcid'),r.get('doi'),r.get('title'));print(r.get('abstractText','')[:1800])
if r.get('pmcid'):
 f=O/'literature/cache'/(r['pmcid']+'_relevant_paragraphs.txt')
 if f.exists():
  for s in f.read_text(encoding='utf-8').split('\n\n'):
   if 'silica' in s.lower() and any(k in s.lower() for k in ['mg','day','week']):print('METHOD',s[:2500])
h=pd.read_csv(O/'direction/human_scope_group_directions.tsv',sep='\t');print('AHNAK_HUMAN');print(h[h.gene=='AHNAK'][['filter','scope','n_control','n_case','difference','P','BH']].to_string(index=False))
a=pd.read_csv(O/'direction/human_leave_one_donor_out.tsv',sep='\t');print('AHNAK_LOO');print(a[(a.gene=='AHNAK')&(a['filter']=='ReviewedAll')][['scope','omitted_donor','full_difference','loo_difference']].to_string(index=False))
c=pd.read_csv(O/'direction/detection_intensity_decomposition.tsv',sep='\t');print('DECOMP');print(c[(c.gene.isin(['MDFIC','TTC39C']))&(c.scope=='ActCD4')&(c['filter']=='All')][['dataset','day','gene','delta','detection_component','positive_intensity_component']].to_string(index=False))