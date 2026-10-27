from pathlib import Path
import re,pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
for filename,segments in [('scripts/14_GSE288887_CD4_damage_pyroptosis_modules.R',[(245,274),(128,158)]),('02_analysis/CD4_DEG_collect_20260914_100836_XJWbue/source/scripts/16_GSE288887_CD4_pathway_analysis.R',[(476,501),(1198,1213)]),('scripts/15_GSE288887_CD4_subtype_module_audit.R',[(210,270),(290,320)])]:
 f=P/filename
 if not f.exists():continue
 ls=f.read_text().splitlines();print(filename)
 for a,b in segments:print('\n'.join(f'{i+1}: {ls[i]}' for i in range(a-1,min(b,len(ls)))))
for root in [P/'logs',P/'07_DEG/CD4_pseudobulk_v1/pathway_v1']:
 print('FILES',str(root),[f.name for f in root.iterdir() if f.is_file() and re.search('pathway|16_|17_|session',f.name,re.I)])
for f in (P/'logs').glob('*'):
 if f.is_file() and re.search(r'16_|pathway|17_',f.name,re.I) and f.stat().st_size<1000000:
  text=f.read_text(errors='replace');lines=[l for l in text.splitlines() if re.search('msigdbr_|clusterProfiler_|GSVA_|UCell_|202[0-9]\.\d\.Hs|Human.*MSig|DB_VERSION|db_version|subjects|Control|CWP',l,re.I)]
  if lines:print('LOG',f.name,'\n'.join(lines[-20:]))
