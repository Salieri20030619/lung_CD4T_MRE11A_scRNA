from pathlib import Path
import subprocess
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/GSE174725_annotation_validation_20261005_v1')
with (O/'provenance/build.log').open('w') as f:p=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python',str(O/'code/g174_build_annotation.py')],stdout=f,stderr=subprocess.STDOUT)
print('EXIT',p.returncode);print((O/'provenance/build.log').read_text()[-6000:]);raise SystemExit(p.returncode)
