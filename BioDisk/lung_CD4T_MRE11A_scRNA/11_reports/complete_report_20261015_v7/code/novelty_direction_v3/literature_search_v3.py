from pathlib import Path
import urllib.request,urllib.parse,json,csv,datetime,concurrent.futures,sys,re
O=Path(sys.argv[1]);cache=O/'literature/cache';cache.mkdir(exist_ok=True)
aliases={'AHNAK':['AHNAK','AHNAK1','desmoyokin'],'MDFIC':['MDFIC','MyoD family inhibitor domain-containing','human I-mfa domain-containing'],'TTC39C':['TTC39C'],'ST6GALNAC3':['ST6GALNAC3','ST6GALNACIII','ST6GalNAc III','SIAT7C']}
topics={'all':'','silica':'(TITLE_ABS:silicosis OR TITLE_ABS:silica OR TITLE_ABS:pneumoconiosis)','fibrosis':'(TITLE_ABS:fibrosis OR TITLE_ABS:fibrotic)','T_cells':'(TITLE_ABS:"T cell" OR TITLE_ABS:"T cells" OR TITLE_ABS:CD4 OR TITLE_ABS:Treg OR TITLE_ABS:lymphocyte)','death_autophagy':'(TITLE_ABS:autophag* OR TITLE_ABS:mitophag* OR TITLE_ABS:mitochondri* OR TITLE_ABS:pyropto* OR TITLE_ABS:apoptos*)'}
jobs=[]
for g,als in aliases.items():
 base='('+ ' OR '.join('TITLE_ABS:"'+x+'"' for x in als)+')'
 for name,t in topics.items():jobs.append((g,name,base+(' AND '+t if t else '')+' AND FIRST_PDATE:[1900-01-01 TO 2026-10-15]'))
 jobs.append((g,'fulltext_silica','('+' OR '.join('"'+x+'"' for x in als)+') AND (silicosis OR silica OR pneumoconiosis) AND FIRST_PDATE:[1900-01-01 TO 2026-10-15]'))
def query(job):
 g,t,q=job;url='https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urllib.parse.urlencode({'query':q,'format':'json','pageSize':100,'resultType':'core','sort':'FIRST_PDATE_D desc'})
 try:
  req=urllib.request.Request(url,headers={'User-Agent':'CD4-research-evidence-review/1.0'});data=json.loads(urllib.request.urlopen(req,timeout=60).read().decode('utf-8'));(cache/(g+'_'+t+'.json')).write_text(json.dumps(data,ensure_ascii=False),encoding='utf-8')
  return {'gene':g,'topic':t,'query':q,'source':'Europe PMC','url':url,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'hit_count':data['hitCount'],'returned':len(data.get('resultList',{}).get('result',[])),'status':'OK'},data.get('resultList',{}).get('result',[])
 except Exception as e:return {'gene':g,'topic':t,'query':q,'status':'ERROR','error':str(e)},[]
logs=[];articles={}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
 for meta,result in ex.map(query,jobs):
  logs.append(meta)
  for r in result:
   key=(meta['gene'],r.get('source'),r.get('id'))
   if key not in articles:articles[key]={'gene':meta['gene'],'source':r.get('source'),'id':r.get('id'),'pmid':r.get('pmid',''),'pmcid':r.get('pmcid',''),'doi':r.get('doi',''),'title':r.get('title',''),'year':r.get('pubYear',''),'first_publication_date':r.get('firstPublicationDate',''),'publication_types':';'.join(r.get('pubTypeList',{}).get('pubType',[])),'topics':set()}
   articles[key]['topics'].add(meta['topic'])
for r in articles.values():r['topics']=';'.join(sorted(r['topics']))
def write(rows,name):
 fields=list(dict.fromkeys(k for r in rows for k in r))
 with (O/'literature'/name).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
write(logs,'search_log.tsv');write(list(articles.values()),'retrieved_article_metadata.tsv')
print('SEARCH_COUNTS')
for r in logs:print(r['gene'],r['topic'],r.get('hit_count'),r['status'])
print('DIRECT_RELEVANCE_TITLES')
for r in articles.values():
 if any(x in r['topics'].split(';') for x in ['silica','T_cells','death_autophagy','fibrosis']):print(r['gene'],r['pmid'],r['year'],r['title'])
# Fetch accession records as text, never raw sequencing data.
for acc in ['GSE314943','GSM9529140','GSE183682','GSM5567502','GSE274372']:
 url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?'+urllib.parse.urlencode({'acc':acc,'targ':'self','form':'text','view':'full'})
 try:
  s=urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0'}),timeout=45).read().decode('utf-8');(O/'direction'/(acc+'_public_metadata.txt')).write_text(s,encoding='utf-8')
  print('GEO',acc,'VALID',('!Sample_' in s or '!Series_' in s))
 except Exception as e:print('GEO_ERROR',acc,str(e))