from pathlib import Path
import os,subprocess,json,fcntl,shutil,hashlib,datetime
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/human_first_20261004_v1')
lock=open(O/'provenance/launch.lock','w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
pidfile=O/'provenance/stage2.pid'
if pidfile.exists():
 pid=int(pidfile.read_text());p=Path('/proc')/str(pid)/'cmdline'
 if p.exists() and 'stage2.R' in p.read_text():raise RuntimeError('Existing stage2 process; will not duplicate')
assert not (O/'provenance/stage2_complete.txt').exists()
assert (O/'provenance/human_prior_manifest.json').exists(),'Human prior not ready'
shutil.copy2(O/'code/STAGE2_SPEC.txt',O/'provenance/STAGE2_SPEC.txt')
(O/'provenance/stage2_frozen.json').write_text(json.dumps({'time':datetime.datetime.now().isoformat(),'code_sha256':hashlib.sha256((O/'code/stage2.R').read_bytes()).hexdigest(),'spec_sha256':hashlib.sha256((O/'provenance/STAGE2_SPEC.txt').read_bytes()).hexdigest()},indent=2))
env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
r='/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript'
subprocess.run([r,'--vanilla','-e',f'parse("{O}/code/stage2.R");cat("PARSE_OK")'],check=True,stdout=subprocess.DEVNULL,env=env)
log=open(O/'provenance/stage2_run.log','w');p=subprocess.Popen([r,'--vanilla',str(O/'code/stage2.R')],stdout=log,stderr=subprocess.STDOUT,start_new_session=True,env=env);pidfile.write_text(str(p.pid));print('STARTED_STAGE2',p.pid)
