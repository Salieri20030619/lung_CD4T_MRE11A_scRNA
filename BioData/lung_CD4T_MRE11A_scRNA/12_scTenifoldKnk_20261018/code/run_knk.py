import os,json,subprocess,time,concurrent.futures
from pathlib import Path
O=Path('/mnt/BioData/bio/wangmingyuan/lung_CD4T_MRE11A_scRNA/12_scTenifoldKnk_20261018')
env=os.environ.copy();env.update(TMPDIR=str(O/'tmp'),R_LIBS_USER=str(O/'library'),OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2',MKL_NUM_THREADS='2')
jobs=json.loads((O/'jobs.json').read_text());r='/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript'
def fit(i):
 j=jobs[i-1];mark=O/'runs'/j['id']/'complete.json'
 if mark.exists() and json.loads(mark.read_text())['status']=='success':return dict(id=j['id'],returncode=0,cached=True)
 with (O/'logs'/(j['id']+'.log')).open('wb') as f:p=subprocess.run([r,'--vanilla',str(O/'code/fit_knk.R'),str(i)],env=env,stdout=f,stderr=subprocess.STDOUT)
 result=dict(id=j['id'],returncode=p.returncode)
 (O/'logs'/(j['id']+'.exit.json')).write_text(json.dumps(result),encoding='utf-8');print(result,flush=True);return result
if __name__=='__main__':
 (O/'supervisor_pid.json').write_text(json.dumps(dict(pid=os.getpid(),started=time.time())),encoding='utf-8')
 pilot=fit(1)
 if pilot['returncode']!=0:raise RuntimeError('Pilot failed; see log')
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:results=[pilot]+list(ex.map(fit,range(2,len(jobs)+1)))
 (O/'compute_complete.json').write_text(json.dumps(dict(status='success' if all(r['returncode']==0 for r in results) else 'error',runs=results),indent=2),encoding='utf-8')