from pathlib import Path
import re,json
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/manuscript_integrated_20261006_v3'
print('ROOTS',[x.name for x in P.iterdir()])
found=[]
for f in P.rglob('*'):
 if f.is_file() and any(x in str(f).lower() for x in ['gsea','pathway_v1','pathway_sensitivity','damage_pyroptosis']) and f.suffix.lower() in ['.tsv','.rds','.r','.json','.txt']:
  found.append({'path':str(f.relative_to(P)),'bytes':f.stat().st_size})
(O/'provenance/enrichment_inventory.json').write_text(json.dumps(found,indent=2))
print('NONFIG_RESULTS',json.dumps([f for f in found if not any(x in f['path'] for x in ['_labels.','source/scripts','context/','figures/'])],ensure_ascii=False))
for name in ['14_GSE288887_CD4_damage_pyroptosis_modules.R','16_GSE288887_CD4_pathway_analysis.R','17_GSE288887_CD4_pathway_sensitivity.R']:
 f=P/'02_analysis/CD4_DEG_collect_20260914_100836_XJWbue/source/scripts'/name;ls=f.read_text().splitlines();(O/'provenance'/name).write_text('\n'.join(ls))
 print('SCRIPT_HEAD',name,'\n'+'\n'.join(f'{i+1}: {s}' for i,s in enumerate(ls[:100])))
 if name.startswith('14'):
  print('APOPT_LINES','\n'.join(f'{i+1}: {s}' for i,s in enumerate(ls) if re.search(r'apopt|ApoExec|module_tests|write_tsv|mito|pyro',s,re.I)))
for f in [P/'02_analysis/GSE314943/stage7_20260905_234200_aiQYL0/tables/methods.txt',P/'02_analysis/GSE314943/stage8_20260906_003153_6mFKmg/tables/methods.txt']:
 if f.exists():print('METHODS',str(f.relative_to(P)),f.read_text())
