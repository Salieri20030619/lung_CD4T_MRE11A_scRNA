from pathlib import Path
import hashlib,json
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
for stem in ['mechanism_step1','mechanism_step2','mechanism_step5','mechanism_step6']:
 for d in P.glob('02_analysis/'+stem+'*'):
  print('TABLES',d.name,[f.name for f in (d/'tables').glob('*')])
for f in P.glob('02_analysis/GSE314943/stage*/objects/*.h5ad'):print('OBJ',f,f.stat().st_size)
print('ENVS',[str(p) for p in Path('/mnt/BioDisk/wangmingyuan/miniconda3/envs').glob('*/bin/python')])