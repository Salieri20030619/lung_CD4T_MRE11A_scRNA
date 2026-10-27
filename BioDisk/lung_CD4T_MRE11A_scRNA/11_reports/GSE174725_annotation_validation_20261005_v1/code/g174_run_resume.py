from pathlib import Path
import subprocess
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/GSE174725_annotation_validation_20261005_v1')
(O/'provenance/resume_reason.txt').write_text('Initial build completed common count export, fixed per-donor Scrublet and graph but marker export failed because scipy sparse indexing does not accept pandas Series. Fixed conversion to numpy boolean. Resume from saved common counts without re-reading CSV; identical doublet and graph parameters. Graph/cell-score checkpoints added. No biological/model specification change.\n')
with (O/'provenance/resume.log').open('w') as f:p=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python',str(O/'code/g174_resume_annotation.py')],stdout=f,stderr=subprocess.STDOUT)
print('EXIT',p.returncode);print((O/'provenance/resume.log').read_text()[-4000:]);raise SystemExit(p.returncode)
