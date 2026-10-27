from pathlib import Path
p=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/02_analysis/mechanism_step1_20260927_152218_s4V5eF/analyze.py')
print('\n'.join(p.read_text().splitlines()[:70]))