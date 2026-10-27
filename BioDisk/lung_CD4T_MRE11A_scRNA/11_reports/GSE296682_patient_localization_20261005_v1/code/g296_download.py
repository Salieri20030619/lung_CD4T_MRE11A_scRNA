from pathlib import Path
import urllib.request,json,hashlib,time,gzip,shutil
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/GSE296682_patient_localization_20261005_v1';C=O/'public_cache';V=O/'provenance'
def get(url):
 with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'ResearchDataAudit/1.0'}),timeout=60) as r:return r.read()
api=V/'zenodo_record_18626573.json'
if not api.exists():api.write_bytes(get('https://zenodo.org/api/records/18626573'))
record=json.loads(api.read_text());entries={x['key']:x for x in record['files']};print('OFFICIAL_FILES',json.dumps([{k:x[k] for k in ['key','size','checksum']} for x in entries.values()]),flush=True)
results=[]
for name in ['BAL_FINAL_META.csv','BAL_FINAL.h5ad']:
 e=entries[name];target=C/name;url='https://zenodo.org/records/18626573/files/'+name+'?download=1';expected=e['checksum'].split(':')[-1]
 if not target.exists():
  assert shutil.disk_usage(C).free>e['size']+2*1024**3
  part=C/(name+'.part');assert not part.exists(),'Partial exists; inspect before resume'
  print('DOWNLOAD_BEGIN',name,e['size'],url,flush=True);start=time.time();total=0;last=start
  with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'ResearchDataAudit/1.0','Accept-Encoding':'identity'}),timeout=60) as r,part.open('wb') as f:
   headers=dict(r.headers);print('HEADERS',json.dumps(headers),flush=True)
   for b in iter(lambda:r.read(1024*1024),b''):
    f.write(b);total+=len(b)
    if time.time()-last>20:print('PROGRESS',name,total,'elapsed',round(time.time()-start,1),flush=True);last=time.time()
  if headers.get('Content-Encoding',headers.get('content-encoding','')).lower()=='gzip':
   with gzip.open(part,'rb') as src,target.open('wb') as dst:shutil.copyfileobj(src,dst)
  else:part.rename(target)
  (V/(name+'.download_headers.json')).write_text(json.dumps(headers,indent=2))
 md=hashlib.md5();sh=hashlib.sha256()
 with target.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):md.update(b);sh.update(b)
 result=dict(file=name,path=str(target),bytes=target.stat().st_size,md5=md.hexdigest(),sha256=sh.hexdigest(),official_md5=expected,url=url);assert result['md5']==expected,result;assert result['bytes']==e['size'],result;results.append(result);print('VERIFIED',json.dumps(result),flush=True)
(V/'official_download_verified.json').write_text(json.dumps(results,indent=2))
