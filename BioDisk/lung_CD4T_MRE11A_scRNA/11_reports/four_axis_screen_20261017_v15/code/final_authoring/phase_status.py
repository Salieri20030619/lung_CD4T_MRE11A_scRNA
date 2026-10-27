from pathlib import Path
import pandas as pd
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/four_axis_screen_20261017_v15')
print((O/'provenance/screen_run.log').read_text()[-1800:]);f=O/'provenance/computation_complete.json'
if f.exists():print(f.read_text())