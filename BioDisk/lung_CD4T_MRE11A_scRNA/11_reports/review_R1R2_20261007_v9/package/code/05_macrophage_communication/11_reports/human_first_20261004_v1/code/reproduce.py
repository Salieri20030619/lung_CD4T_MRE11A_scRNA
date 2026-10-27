"""Re-run in a NEW server directory using verified existing inputs and installed packages."""
from pathlib import Path
import sys,os,shutil,subprocess,json,hashlib,datetime
src=Path(__file__).resolve().parent.parent
if len(sys.argv)!=2:raise SystemExit('Usage: python code/reproduce.py /mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/human_first_repro_NAME')
dst=Path(sys.argv[1]).resolve();root=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports').resolve()
assert dst.parent==root and dst.name.startswith('human_first_repro_') and not dst.exists(),'Choose a new human_first_repro_* directory directly below project/11_reports'
for r in json.loads((src/'provenance/human_prior_manifest.json').read_text()):assert hashlib.sha256(Path(r['path']).read_bytes()).hexdigest()==r['sha256']
for d in ['code','results','provenance','figures','deliverables','resources']:(dst/d).mkdir(parents=True,exist_ok=True)
for f in (src/'code').iterdir():
 if f.is_file() and f.suffix in ['.py','.R','.txt']:(dst/'code'/f.name).write_text(f.read_text().replace(src.name,dst.name))
shutil.copy2(src/'provenance/human_prior_manifest.json',dst/'provenance/human_prior_manifest.json');shutil.copy2(dst/'code/STAGE2_SPEC.txt',dst/'provenance/STAGE2_SPEC.txt')
(dst/'provenance/reproduction_origin.json').write_text(json.dumps({'source_run':str(src),'time':datetime.datetime.now().isoformat()},indent=2))
py='/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python';r='/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript'
env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
for name in ['stage1.py','stage2.R','stage2_followup.py','stage3_mapping.R','mouse_background.R','stage3.py','sender_context.R','enrich_and_validate.py','build_report.py']:
 cmd=[r,'--vanilla',str(dst/'code'/name)] if name.endswith('.R') else [py,str(dst/'code'/name)]
 print('RUN',name,flush=True)
 with (dst/'provenance'/('reproduce_'+name+'.log')).open('w') as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,env=env,check=True)
print('REPRODUCTION_COMPLETE',dst,flush=True)
