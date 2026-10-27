from pathlib import Path
import subprocess,os,sys,json
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/manuscript_graphics_20261006_v5')
name=sys.argv[1];assert name in ['graphics_inputs.R','graphics_gsea.R','graphics_umap.R','graphics_mouse_gsea.R']
script=str(O/'code'/name);ps=subprocess.check_output(['ps','-eo','args'],text=True)
if any(script in s and 'graphics_run.py' not in s for s in ps.splitlines()):raise SystemExit('ALREADY_RUNNING')
env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
with (O/'provenance'/f'{name}.log').open('w') as f:
 p=subprocess.Popen(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript',script],stdout=f,stderr=subprocess.STDOUT,start_new_session=True,env=env)
print(json.dumps({'pid':p.pid,'script':name,'log':str(O/'provenance'/f'{name}.log')}))
