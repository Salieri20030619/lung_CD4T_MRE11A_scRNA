from pathlib import Path
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');M=P/'02_analysis/mechanism_step1_20260927_152218_s4V5eF'
for f in M.glob('*.py'):
 s=f.read_text()
 for line in s.splitlines():
  if any(w in line for w in ['BACKGROUND','BG=','BG =','bggenes','BGMARK','SFTPC','bg_frac','background_fraction']):print(f.name,line[:500])
print('CODE_FILES',[f.name for f in M.glob('*.py')])