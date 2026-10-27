from pathlib import Path
import sys,urllib.request,urllib.parse,json,xml.etree.ElementTree as ET,re,csv
O=Path(sys.argv[1]);targets=[('MDFIC','28341862'),('MDFIC','25437532'),('ST6GALNAC3','32596346'),('TTC39C','19936710'),('AHNAK','21624382'),('AHNAK','29724814')]
for gene,pmid in targets:
 url='https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urllib.parse.urlencode({'query':'EXT_ID:'+pmid+' AND SRC:MED','format':'json','resultType':'core'})
 r=json.loads(urllib.request.urlopen(url,timeout=40).read())['resultList']['result'][0];(O/'literature/cache'/('PMID_'+pmid+'.json')).write_text(json.dumps(r,ensure_ascii=False),encoding='utf-8');print('PAPER',gene,pmid,r.get('pmcid'),r.get('doi'),r.get('title'))
 if not r.get('pmcid'):continue
 try:
  b=urllib.request.urlopen('https://www.ebi.ac.uk/europepmc/webservices/rest/'+r['pmcid']+'/fullTextXML',timeout=50).read();(O/'literature/cache'/(r['pmcid']+'.xml')).write_bytes(b);root=ET.fromstring(b);hits=[]
  for tag in ['p','table-wrap']:
   for node in root.iter(tag):
    s=' '.join(node.itertext())
    if re.search(r'(?<![A-Za-z0-9])'+gene+r'(?![A-Za-z0-9])',s,re.I):
     hits.append(s)
     m=re.search(gene,s,re.I);print('CONTEXT',s[max(0,m.start()-350):m.end()+1000])
  (O/'literature/cache'/(r['pmcid']+'_'+gene+'_context.txt')).write_text('\n\n'.join(hits),encoding='utf-8')
  if not hits:print('NO_EXACT_GENE_IN_MAIN_TEXT')
 except Exception as e:print('FETCH_ERROR',str(e))