from pathlib import Path
import json,requests,concurrent.futures,datetime,time
O=Path(r'C:\Users\31937\lung_three_cohort_20261017_v16');(O/'literature_success').mkdir(exist_ok=True)
s=(O/'inspect_candidates.run.log').read_text(encoding='utf-8');line=next(x for x in s.splitlines() if x.startswith('DE_FUNCTION '));genes=set(line.split(' ',2)[2].split(';'))|{'IL1R2','DUSP5','STAT1','TNFRSF9','GSN','FOXP1'}
import pandas as pd
genes=set(pd.read_csv(O/'data/tables/three_axis_candidates_all_contrasts.tsv',sep='\t').gene)
alias={'TNFRSF1A':['TNFR1'],'TNFRSF1B':['TNFR2'],'TNFRSF9':['4-1BB','CD137'],'GSN':['gelsolin'],'PMAIP1':['NOXA'],'BCL2L11':['BIM'],'BCL2L1':['BCL-XL'],'IL1R2':['IL-1R2'],'HIF1A':['HIF-1alpha'],'DDIT4':['REDD1'],'BCL2A1A':['BCL2A1','BFL1'],'BCL2A1B':['BCL2A1','BFL1'],'BCL2A1D':['BCL2A1','BFL1'],'CASP4':['caspase-11','CASP11'],'CDKN1A':['p21'],'PARK7':['DJ-1'],'TNFSF10':['TRAIL']}
alias.update({'CASP1':['caspase-1'],'SKIL':['SnoN'],'ERN1':['IRE1alpha','IRE1-alpha'],'CTNNB1':['beta-catenin'],'GZMB':['granzyme B'],'PPM1K':['PP2Cm'],'FAM162A':['HGTD-P'],'BCL3':['Bcl-3'],'MAP2K3':['MKK3'],'MAP3K1':['MEKK1'],'MAPK9':['JNK2'],'PTPN1':['PTP1B'],'TPT1':['TCTP'],'IER3':['IEX-1'],'ZEB2':['SIP1']})
queries=[]
for g in sorted(genes):
 term='('+' OR '.join('TITLE_ABS:"'+v+'"' for v in [g]+alias.get(g,[]))+')'
 for kind,context in [('disease','(TITLE_ABS:silicosis OR TITLE_ABS:"silica-induced" OR TITLE_ABS:"silica induced")'),('mechanism','(TITLE_ABS:pyroptosis OR TITLE_ABS:"mitochondrial apoptosis" OR TITLE_ABS:"intrinsic apoptosis" OR (TITLE_ABS:mitochondria AND TITLE_ABS:apoptosis))'),('CD4_death','(TITLE_ABS:CD4 OR TITLE_ABS:Treg OR TITLE_ABS:"regulatory T") AND (TITLE_ABS:apoptosis OR TITLE_ABS:pyroptosis OR TITLE_ABS:mitochondrial)')]:
  queries.append((g,kind,term+' AND '+context+' AND FIRST_PDATE:[1900-01-01 TO 2026-10-17] NOT PUB_TYPE:review'))
def work(item):
 g,k,q=item;p=O/'literature_success'/f'{g}__{k}.json'
 if p.exists():
  old=json.loads(p.read_text(encoding='utf-8'))
  if old['query']==q and old['status']=='success':return old
  p.with_suffix('.previous.json').write_bytes(p.read_bytes())
 result={'gene':g,'kind':k,'query':q,'searched_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'endpoint':'https://www.ebi.ac.uk/europepmc/webservices/rest/search'}
 for attempt in range(3):
  try:
   r=requests.get(result['endpoint'],params={'query':q,'format':'json','pageSize':10,'resultType':'core'},timeout=(10,35));r.raise_for_status();j=r.json();result.update(hitCount=j['hitCount'],results=j.get('resultList',{}).get('result',[]),status='success');break
  except Exception as e:
   result.update(status='failed',error=str(e));time.sleep(1)
 p.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8');return result
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as ex:
 out=list(ex.map(work,queries))
summary=[{k:r.get(k) for k in ['gene','kind','hitCount','status','query','searched_at']} for r in out];(O/'literature/search_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8');print('GENES',len(genes),'QUERIES',len(out),'FAILED',sum(r['status']!='success' for r in out))