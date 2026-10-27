from pathlib import Path
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/four_axis_screen_20261017_v15')
p=O/'provenance/decontx_diagnostic.log';print(p.read_text()[-2000:] if p.exists() else 'exporting input');print('DONE',(O/'provenance/decontx_exit.json').exists())