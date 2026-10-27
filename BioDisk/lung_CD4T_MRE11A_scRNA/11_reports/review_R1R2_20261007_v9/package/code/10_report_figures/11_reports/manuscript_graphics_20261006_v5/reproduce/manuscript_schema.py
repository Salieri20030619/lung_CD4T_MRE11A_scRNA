from pathlib import Path
import pandas as pd,re
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');T=P/'11_reports/technical_delivery_20261002_ikzf2_stat1_v1/revision_v2/final_release'
for root in [P/'10_tables/CD4_final/module_audit_v1',P/'10_tables/CD4_final/subtype_module_audit_v1',P/'11_reports/stepwise_diagnostic_20261002_v1/results']:
 print('DIR',root);print('\n'.join(f.name for f in root.iterdir() if f.is_file()))
files=[P/'10_tables/CD4_final/pathway_v1/GO_KEGG_GSEA.tsv',P/'10_tables/CD4_final/pathway_v1/GO_KEGG_GSVA_tests.tsv',P/'10_tables/CD4_final/pathway_sensitivity_v1/mitochondrial_pathway_sensitivity.tsv',P/'10_tables/CD4_final/pathway_sensitivity_v1/mitochondrial_pathway_stability.tsv',P/'02_analysis/GSE314943/stage7_20260905_234200_aiQYL0/tables/Core_Full_Standard_GO_BP_GSEA.tsv',T/'source_tables/M1_group_effects.tsv',T/'data/S2_TF_associations.tsv',T/'source_tables/M2_human_source_models.tsv']
for f in files:
 d=pd.read_csv(f,sep='\t');print('SCHEMA',str(f.relative_to(P)),d.shape,list(d));print(d.head(2).to_string(index=False)[:4500])
 for c in ['dataset','scope','model','filter','feature','metric','module','species','variant']:
  if c in d:print(c,d[c].drop_duplicates().tolist()[:30])
for f in [P/'10_tables/CD4_final/pathway_sensitivity_v1/sensitivity_GSEA_summary.tsv',P/'10_tables/CD4_final/pathway_v1/pathway_DEG_summary.tsv',P/'10_tables/CD4_final/pathway_sensitivity_v1/rank_input_manifest.tsv',P/'10_tables/CD4_final/pathway_sensitivity_v1/mitochondrial_leading_edge_summary.tsv']:
 print('TABLE',str(f.relative_to(P)));print(f.read_text()[:5000])
