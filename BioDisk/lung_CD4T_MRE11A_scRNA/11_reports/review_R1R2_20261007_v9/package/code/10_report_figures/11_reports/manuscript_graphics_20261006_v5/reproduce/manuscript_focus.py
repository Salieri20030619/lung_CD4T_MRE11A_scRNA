from pathlib import Path
import pandas as pd,json,re,hashlib
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');B=P/'11_reports';T=B/'technical_delivery_20261002_ikzf2_stat1_v1/revision_v2/final_release';O=B/'manuscript_integrated_20261006_v3'
for sub in ['code','results','provenance','figures','deliverables']:(O/sub).mkdir(parents=True,exist_ok=True)
roots=[P/'02_analysis/CD4_DEG_collect_20260914_100836_XJWbue/source/10_tables/CD4_final',P/'02_analysis/GSE314943/stage8_20260906_003153_6mFKmg',T/'data',T/'source_tables']
for root in roots:
 print('DIRECTORY',str(root.relative_to(P)))
 print('\n'.join(str(f.relative_to(root))+' '+str(f.stat().st_size) for f in root.rglob('*') if f.is_file() and f.suffix in ['.tsv','.R','.py','.md','.json','.txt'] and 'sample' not in f.name.lower() and 'donor' not in f.name.lower()))
for f in [T/'source_tables/M1_module_definition.tsv',T/'data/IKZF2_Core_death_module_direction_audit.tsv',P/'02_analysis/GSE314943/stage8_20260906_003153_6mFKmg/tables/module_definitions.tsv',P/'02_analysis/GSE314943/stage7_20260905_234200_aiQYL0/tables/GSEA_summary.tsv']:
 print('TABLE',str(f.relative_to(P)));print(f.read_text()[:10000])
for f in [T/'data/M2_all_module_group_models.tsv',T/'source_tables/S2_TF_module_associations.tsv',T/'source_tables/M1_group_models.tsv']:
 if f.exists():
  d=pd.read_csv(f,sep='\t');print('SCHEMA',str(f.relative_to(P)),d.shape,list(d));print(d.head(2).to_string(index=False))
