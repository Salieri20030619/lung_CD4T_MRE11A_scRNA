from pathlib import Path
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
f=P/'02_analysis/CD4_DEG_collect_20260914_100836_XJWbue/source/scripts/16_GSE288887_CD4_pathway_analysis.R';ls=f.read_text().splitlines()
print('\n'.join(f'{i+1}: {ls[i]}' for i in list(range(350,430))+list(range(485,520))))
for root in [P/'10_tables/CD4_final/pathway_v1',P/'04_objects']:
 if root.exists():print('FILES',root,[str(x.relative_to(root)) for x in root.rglob('*') if any(w in x.name.lower() for w in ['sample','design','pseudobulk','status','summary'])][:80])
for f in (P/'11_reports').rglob('*.py'):
 if f.stat().st_size>200000:continue
 txt=f.read_text(errors='replace')
 if ('S2_TF_associations' in txt and 'groupby' in txt and ('fdr_bh' in txt or 'bh(' in txt)) or ('tf_associations.tsv' in txt.lower() and 'groupby' in txt):
  print('FDR_CODE',f)
  for i,l in enumerate(txt.splitlines()):
   if any(w in l for w in ['groupby','FDR','associations.tsv']):print(i+1,l[:230])
