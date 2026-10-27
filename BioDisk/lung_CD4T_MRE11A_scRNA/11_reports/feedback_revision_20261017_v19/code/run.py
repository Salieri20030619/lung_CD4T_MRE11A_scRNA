import subprocess,sys,pathlib
base='/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/feedback_revision_20261017_v19/code/'
f=pathlib.Path(sys.argv[1]);dest=base+f.name
flags={'creationflags':subprocess.CREATE_NO_WINDOW}
subprocess.run(['scp','-q',str(f),'sgs10:'+dest],check=True,**flags)
exe='/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript --vanilla' if f.suffix=='.R' else '/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python -u'
if f.name=='aucell_clean.py':exe='/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/software/pyscenic_0.12.1_py310/bin/python -u'
r=subprocess.run(['ssh','-T','-o','BatchMode=yes','sgs10',exe+' '+dest],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,**flags); (f.parent/(f.stem+'.run.log')).write_bytes(r.stdout); sys.stdout.buffer.write(r.stdout); sys.exit(r.returncode)