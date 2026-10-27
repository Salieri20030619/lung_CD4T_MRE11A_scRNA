#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA
mkdir -p "$ROOT/logs" "$ROOT/00_raw" "$ROOT/tmp"
export ROOT
exec 9>"$ROOT/tmp/download_lung_screen6.lock"
flock -n 9 || { echo DOWNLOAD_ALREADY_RUNNING; exit 0; }
TASK=$(mktemp -d "$ROOT/tmp/download_screen6_XXXXXX")
export TASK
trap 'rc=$?; echo "SCREEN6_DOWNLOAD_FAIL EXIT_CODE=$rc"; exit "$rc"' ERR
cat > "$TASK/download.py" <<'PY'
import os,subprocess,gzip,tarfile,concurrent.futures,csv,time,re
from pathlib import Path
from urllib.parse import urlparse
root=Path(os.environ['ROOT']);task=Path(os.environ['TASK'])
env=os.environ.copy()
proxy=os.environ.get('SCREEN_PROXY','http://127.0.0.1:17897')
for k in ['http_proxy','https_proxy','HTTP_PROXY','HTTPS_PROXY']:env[k]=proxy
# Explicit proxy selection prevents inherited no_proxy from bypassing this transfer route.
base=['curl','--proxy',proxy,'--noproxy','','-L','--fail','--show-error','--silent','--connect-timeout','30']
url='https://ftp.ncbi.nlm.nih.gov/geo/README.txt'
subprocess.run(base+['--max-time','90','-o',str(task/'proxy_test.txt'),url],env=env,check=True)
print('PROXY_TEST_PASS',flush=True)
series=['GSE184854','GSE183682','GSE296682','GSE295463','GSE174725','GSE159354']
manifest=[]
def fetch_small(url,dest):
 subprocess.run(base+['--retry','4','--retry-delay','5','--max-time','180','-o',str(dest)+'.part',url],env=env,check=True)
 os.replace(str(dest)+'.part',dest)
def soft_parse(path):
 groups={};current=None;seriesfiles=[]
 with gzip.open(path,'rt',errors='replace') as f:
  for line in f:
   line=line.rstrip('\r\n')
   if line.startswith('^SAMPLE = '):current=line.split(' = ',1)[1];groups[current]=[]
   elif line.startswith('^'):current=None
   if current:groups[current].append(line)
   elif line.startswith('!Series_supplementary_file = '):seriesfiles.append(line.split(' = ',1)[1])
 return groups,seriesfiles
def urls(lines):
 return [x.split(' = ',1)[1] for x in lines if x.startswith('!Sample_supplementary_file') and ' = ' in x and x.split(' = ',1)[1].startswith(('ftp://','https://'))]
for g in series:
 dest=root/'00_raw'/g;dest.mkdir(exist_ok=True)
 soft=dest/(g+'_family.soft.gz')
 fetch_small(f'https://ftp.ncbi.nlm.nih.gov/geo/series/{g[:-3]}nnn/{g}/soft/{g}_family.soft.gz',soft)
 groups,sfiles=soft_parse(soft)
 selected=[];scope='all public sample supplementary files'
 if g=='GSE296682':
  selected=[u for u in sfiles if u.endswith('_BAL_FINAL.h5ad') or u.endswith('_BAL_FINAL_metadata.txt.gz')]
  scope='processed h5ad and cell metadata'
 elif g=='GSE174725':
  selected=[u for u in sfiles if u.endswith(('_RAW.tar','_barcodes.tsv.gz','_features.tsv.gz','_matrix.mtx.gz'))]
  scope='series processed matrix and original supplemental archive'
 elif g=='GSE159354':
  chosen={s:l for s,l in groups.items() if any('CD45' in x.upper() and not re.search(r'CD45\s*(?:[-−]|NEG)',x.upper()) for x in l if x.startswith(('!Sample_title','!Sample_characteristics_ch1')))}
  if not chosen:raise RuntimeError(g+': no CD45 samples identified; inspect family.soft.gz')
  for s,l in chosen.items():selected.extend(urls(l))
  scope='CD45 samples selected from title/characteristics: '+','.join(chosen)
 else:
  for l in groups.values():selected.extend(urls(l))
  selected.extend(u for u in sfiles if u.endswith('.xlsx'))
 if not selected:raise RuntimeError(g+': no matching public files; inspect family.soft.gz')
 for u in sorted(set(selected)):
  u=u.replace('ftp://','https://',1)
  name=Path(urlparse(u).path).name
  manifest.append({'dataset':g,'url':u,'path':str(dest/name),'scope':scope})
 print(g,'FILES='+str(len(set(selected))),scope,flush=True)
with (root/'logs/screen6_download_manifest.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=['dataset','url','path','scope'],delimiter='\t');w.writeheader();w.writerows(manifest)
def valid(p):
 if not p.exists() or p.stat().st_size==0:return False
 try:
  if p.name.endswith('.gz'):
   with gzip.open(p,'rb') as f:
    while f.read(4*1024*1024):pass
  elif p.name.endswith('.tar'):
   with tarfile.open(p,'r:') as t:
    if not any(True for _ in t):return False
  elif p.name.endswith(('.h5ad','.h5')):
   with p.open('rb') as f:
    if f.read(8)!=b'\x89HDF\r\n\x1a\n':return False
  elif p.name.endswith('.xlsx'):
   import zipfile
   with zipfile.ZipFile(p) as z:
    if z.testzip() is not None:return False
  return True
 except Exception:return False
def job(row):
 p=Path(row['path']);log=root/'logs'/(row['dataset']+'_'+p.name+'.download.log')
 if valid(p):print('EXISTS_PASS='+p.name,flush=True);return row['dataset'],p.name,'PASS'
 if p.exists():p.rename(p.with_name(p.name+'.invalid_'+str(time.time_ns())))
 part=p.with_name(p.name+'.part')
 print('DOWNLOAD_START='+p.name,flush=True)
 # Each attempt starts a fresh curl process and resumes from the current byte offset.
 for attempt in range(1,21):
  with log.open('a') as f:
   f.write(f'ATTEMPT={attempt}\n');f.flush()
   rc=subprocess.run(base+['--continue-at','-','--speed-limit','1024','--speed-time','180','--max-time','7200','-o',str(part),row['url']],env=env,stdout=f,stderr=f).returncode
  if rc==0:
   # Validate with the real extension; retain corrupt transfers for diagnosis.
   os.replace(part,p)
   if valid(p):print('DOWNLOAD_PASS='+p.name,flush=True);return row['dataset'],p.name,'PASS'
   p.rename(p.with_name(p.name+'.invalid_'+str(time.time_ns())))
   print('VALIDATION_RETRY='+p.name,flush=True)
  if rc==33:
   print('RANGE_UNSUPPORTED='+p.name,flush=True)
   return row['dataset'],p.name,'FAIL_RANGE'
  time.sleep(10)
 return row['dataset'],p.name,'FAIL'
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:results=list(pool.map(job,manifest))
with (root/'logs/screen6_download_audit.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['dataset','file','status']);w.writerows(results)
for g in series:
 r=[x for x in results if x[0]==g]
 print(g+'_DOWNLOAD_'+('PASS' if r and all(x[2]=='PASS' for x in r) else 'FAIL'),flush=True)
if not all(x[2]=='PASS' for x in results):raise RuntimeError('Some files failed; see screen6_download_audit.tsv')
print('SCREEN6_DOWNLOAD_PASS',flush=True)
PY
bash -n "$0"
conda run -n batch_bench python -m py_compile "$TASK/download.py"
echo SYNTAX_CHECK=PASS
conda run --no-capture-output -n batch_bench python -u "$TASK/download.py"
