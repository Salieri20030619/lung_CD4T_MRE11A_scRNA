from pathlib import Path
import subprocess,json,hashlib,shutil,time
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/GSE296682_patient_localization_20261005_v1';C=O/'public_cache';V=O/'provenance'
entries=[('BAL_FINAL_META.csv','89a6e220d3ea1da5418027c6fd71f8b0',28044314),('BAL_FINAL.h5ad','86e4165b84f965e807e4cf361d3d6bb6',None)];results=[]
(V/'official_resource_snapshot.json').write_text(json.dumps(dict(source='https://zenodo.org/records/18626573',version='v2',retrieved_UTC='2026-10-05',files=[dict(file=n,md5=m,bytes=z) for n,m,z in entries],note='Official web record checksums; API timed out and was stopped. Official file HEAD succeeded. No alternate host or access-control bypass.'),indent=2))
for name,expected,size in entries:
 target=C/name;url='https://zenodo.org/records/18626573/files/'+name+'?download=1'
 if not target.exists():
  assert shutil.disk_usage(C).free>3*1024**3
  part=C/(name+'.part');assert not part.exists(),'Inspect existing partial before resume'
  print('DOWNLOAD_BEGIN',name,url,flush=True)
  proc=subprocess.run(['curl','--fail','-L','--connect-timeout','15','--max-time','900','--silent','--show-error','--output',str(part),'--write-out','HTTP=%{http_code} BYTES=%{size_download} SPEED=%{speed_download}\n',url],capture_output=True,text=True,timeout=920)
  (V/(name+'.curl_status.txt')).write_text(proc.stdout+'\n'+proc.stderr);print(proc.stdout,proc.stderr,flush=True);assert proc.returncode==0,proc.returncode;part.rename(target)
 md=hashlib.md5();sh=hashlib.sha256()
 with target.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):md.update(b);sh.update(b)
 r=dict(file=name,path=str(target),bytes=target.stat().st_size,md5=md.hexdigest(),sha256=sh.hexdigest(),official_md5=expected,url=url);assert r['md5']==expected,r
 if size:assert r['bytes']==size
 results.append(r);(V/'official_download_verified.json').write_text(json.dumps(results,indent=2));print('VERIFIED',json.dumps(r),flush=True)
