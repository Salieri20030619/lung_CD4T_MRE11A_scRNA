from pathlib import Path
import pandas as pd,re
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
for root in [P/'10_tables/CD4_final/module_audit_v1',P/'10_tables/CD4_final/subtype_module_audit_v1']:
 for f in root.glob('*.tsv'):
  if any(k in f.name for k in ['cell','subject_scores','detection','correlations','counts']):continue
  d=pd.read_csv(f,sep='\t');print('TABLE_SCHEMA',str(f.relative_to(P)),d.shape,list(d));print(d.head(2).to_string(index=False)[:1200])
for name,segments in [('14_GSE288887_CD4_damage_pyroptosis_modules.R',[(159,242),(295,330),(337,380),(390,435),(470,510)]),('16_GSE288887_CD4_pathway_analysis.R',[(350,469),(632,670),(1130,1202),(1239,1295)]),('17_GSE288887_CD4_pathway_sensitivity.R',[(446,468),(590,607)])]:
 f=P/'scripts'/name
 if not f.exists():f=P/'02_analysis/CD4_DEG_collect_20260914_100836_XJWbue/source/scripts'/name
 ls=f.read_text().splitlines();print('CODE',name)
 for a,b in segments:print('\n'.join(f'{i+1}: {ls[i]}' for i in range(a-1,b)))
for root in [P/'10_tables/CD4_final/pathway_v1',P/'02_analysis/GSE314943/stage7_20260905_234200_aiQYL0/tables',P/'02_analysis/GSE314943/stage8_20260906_003153_6mFKmg/tables']:
 for f in root.glob('*session*'):
  print('SESSION',f);print('\n'.join(x for x in f.read_text().splitlines() if any(k in x for k in ['clusterProfiler','fgsea','msigdbr','GSVA','UCell','org.Mm','GO.db','R version'])))
