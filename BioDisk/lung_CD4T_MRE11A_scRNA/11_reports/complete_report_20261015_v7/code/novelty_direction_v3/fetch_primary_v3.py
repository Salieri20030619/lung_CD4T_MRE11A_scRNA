from pathlib import Path
import urllib.request,urllib.parse,json,xml.etree.ElementTree as ET,re,sys,pandas as pd,concurrent.futures
O=Path(sys.argv[1]);cache=O/'literature/cache';ids=['41862476','38349059','28223352','34461439','18191595','33987830','36303158','31188487','35235341','41526926','32457453']
allrecords=[]
for pmid in ids:
 u='https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urllib.parse.urlencode({'query':'EXT_ID:'+pmid+' AND SRC:MED','format':'json','resultType':'core'})
 d=json.loads(urllib.request.urlopen(u,timeout=40).read());r=d['resultList']['result'][0];allrecords.append(r);(cache/('PMID_'+pmid+'.json')).write_text(json.dumps(r,ensure_ascii=False),encoding='utf-8')
 if pmid in ['41862476','38349059']:print('PRIMARY',pmid,r.get('pmcid'),r.get('doi'),r.get('title'))
for r in allrecords:
 if r.get('pmid') not in ['38349059','41862476','28223352','36303158']:continue
 pmcid=r.get('pmcid')
 if not pmcid:continue
 try:
  b=urllib.request.urlopen('https://www.ebi.ac.uk/europepmc/webservices/rest/'+pmcid+'/fullTextXML',timeout=45).read();(cache/(pmcid+'.xml')).write_bytes(b);root=ET.fromstring(b)
  hits=[]
  for node in root.iter('p'):
   s=' '.join(node.itertext())
   keys=['MDFIC','A549'] if r['pmid']=='38349059' else ['silica','SiO','28 days','56 days'] if r['pmid']=='41862476' else ['MDFIC','Ttc39c']
   if any(x.lower() in s.lower() for x in keys):hits.append(s)
  (cache/(pmcid+'_relevant_paragraphs.txt')).write_text('\n\n'.join(hits),encoding='utf-8')
  print('FULLTEXT',r['pmid'],len(hits))
  if r['pmid']=='38349059':
   for s in hits:
    if ('MDFIC' in s and any(x in s for x in ['TGF','A549','knockdown','EMT'])):print(s[:2600])
  if r['pmid']=='41862476':
   for s in hits:
    if any(x in s.lower() for x in ['mg/kg','mg/ml','silica-induced','silica suspension','after','day']):print(s[:2500])
 except Exception as e:print('FULLTEXT_ERROR',pmcid,str(e))