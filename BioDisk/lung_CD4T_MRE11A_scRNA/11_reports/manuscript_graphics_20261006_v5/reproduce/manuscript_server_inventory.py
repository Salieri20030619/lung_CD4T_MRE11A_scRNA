from pathlib import Path
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');H=P/'11_reports/human_first_20261004_v1/results';V=P/'11_reports/phase_delivery_IKZF2_STAT1_20261006_v2'
print('H',[(f.name,f.stat().st_size) for f in H.iterdir() if f.is_file() and ('member' in f.name or 'background' in f.name or 'gene' in f.name or 'DE_' in f.name)])
print('OLD',[(str(f.relative_to(V)),f.stat().st_size) for f in V.rglob('*') if f.is_file() and ('provenance' in f.parts or f.parent.name=='code')])
print('MOUSEGENES',[(str(f),f.stat().st_size) for f in (P/'02_analysis/mechanism_step5_20260928_153825_Q4sB1x/tables').glob('*members*')])
