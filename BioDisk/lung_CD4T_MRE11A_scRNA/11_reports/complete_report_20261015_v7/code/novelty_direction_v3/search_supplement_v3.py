from pathlib import Path
import urllib.request,urllib.parse,json,csv,sys,datetime,concurrent.futures
O=Path(sys.argv[1]);jobs=[]
aliases={'AHNAK':'(AHNAK OR AHNAK1 OR desmoyokin)','MDFIC':'MDFIC','TTC39C':'(TTC39C OR C18orf17 OR HsT2697)','ST6GALNAC3':'(ST6GALNAC3 OR ST6GALNACIII OR "ST6GalNAc III" OR SIAT7C)'}
for g,b in aliases.items():
 for t,k in [('strict_silicosis_fulltext','(silicosis OR pneumoconiosis OR "silica-induced pulmonary")'),('CD4_fulltext','(CD4 OR "T cell" OR "T cells")')]:jobs.append((g,t,b+' AND '+k+' AND FIRST_PDATE:[1900-01-01 TO 2026-10-15]'))
jobs.append(('TTC39C','legacy_alias_titleabs','(TITLE_ABS:C18orf17 OR TITLE_ABS:HsT2697) AND FIRST_PDATE:[1900-01-01 TO 2026-10-15]'))
def query(job):
 g,t,q=job;u='https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urllib.parse.urlencode({'query':q,'format':'json','pageSize':100,'resultType':'core'})
 d=json.loads(urllib.request.urlopen(u,timeout=60).read());(O/'literature/cache'/(g+'_'+t+'.json')).write_text(json.dumps(d,ensure_ascii=False),encoding='utf-8');return dict(gene=g,topic=t,query=q,url=u,hit_count=d['hitCount'],returned=len(d['resultList']['result']),retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),d['resultList']['result']
meta=[];papers=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
 for m,rs in ex.map(query,jobs):
  meta.append(m)
  for r in rs:papers.append(dict(gene=m['gene'],topic=m['topic'],pmid=r.get('pmid',''),pmcid=r.get('pmcid',''),doi=r.get('doi',''),title=r.get('title',''),year=r.get('pubYear','')))
for rows,name in [(meta,'supplemental_search_log.tsv'),(papers,'supplemental_article_metadata.tsv')]:
 with (O/'literature'/name).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
for m in meta:print(m['gene'],m['topic'],m['hit_count'],m['returned'])
print('STRICT_DISEASE_HITS')
for r in papers:
 if r['topic']=='strict_silicosis_fulltext':print(r['gene'],r['pmid'],r['title'])