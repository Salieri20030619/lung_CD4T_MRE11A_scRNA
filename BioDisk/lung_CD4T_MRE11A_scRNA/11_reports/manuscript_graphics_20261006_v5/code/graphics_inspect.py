from pathlib import Path
import json,re
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/manuscript_graphics_20261006_v5'
inv=json.loads((O/'provenance/object_inventory.json').read_text())
print('RELEVANT_OBJECTS',json.dumps([x for x in inv if x['path'].startswith(('04_pre','06_CD4','07_DEG','10_tables','02_analysis/GSE314943')) and not re.search('geneset_CD|zoneinfo',x['path'])]))
print('ENVS',[x.name for x in Path('/mnt/BioDisk/wangmingyuan/miniconda3/envs').iterdir()])
for root in ['10_tables/broad_clustering_v1','10_tables/CD4_final','02_analysis/GSE314943/stage7_20260905_234200_aiQYL0','11_reports/technical_delivery_20261002_ikzf2_stat1_v1/revision_v2/final_release/provenance']:
 d=P/root
 print('DIRECTORY',root,[x.name for x in d.iterdir()] if d.exists() else 'missing')
for name in ['scripts/16_GSE288887_CD4_pathway_analysis.R','scripts/17_GSE288887_CD4_pathway_sensitivity.R','02_analysis/GSE314943/stage7_20260905_234200_aiQYL0/analysis.R']:
 f=P/name;ls=f.read_text().splitlines();print('CODE',name)
 print('\n'.join(f'{i+1}: {l}' for i,l in enumerate(ls) if re.search('readRDS|read_tsv|read.csv|saveRDS|geneList|ranking|ranked|rank_vec|gsea|GSEA|TERM2|gmt|gene_sets|pathway|stats|sqrt|FDR',l)))
for f in (P/'scripts').glob('*.R'):
 if re.search('broad|CD4_final',f.name,re.I):
  print('ANNOTATION_CODE',f.name);ls=f.read_text().splitlines();print('\n'.join(f'{i+1}: {l}' for i,l in enumerate(ls) if re.search('readRDS|saveRDS|lineage|celltype|CD4_sub|write_tsv|UMAP|umap|Subject|source',l)))
