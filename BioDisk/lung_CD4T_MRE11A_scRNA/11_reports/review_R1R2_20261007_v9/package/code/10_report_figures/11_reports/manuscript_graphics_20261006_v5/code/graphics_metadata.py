from pathlib import Path
import pandas as pd,json
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/manuscript_graphics_20261006_v5'
for n in ['Human_CD4_UMAP_SERVER_ONLY.tsv','Human_lung_UMAP_SERVER_ONLY.tsv']:
 d=pd.read_csv(O/'tables'/n,sep='\t',low_memory=False)
 print('META',n,'shape',d.shape)
 for c in ['dataset','sample_id','subject_id','group','source_type','cd4_label','cd4_label_all','cd4_decision','lineage_prelim']:
  if c in d:print(c,d[c].value_counts(dropna=False).to_dict())
for f in [P/'02_analysis/mechanism_step2_20260927_154919_EsN4gk/input/human_cells.tsv',P/'10_tables/broad_clustering_v1/GSE288887_preliminary_lineage_summary.tsv']:
 d=pd.read_csv(f,sep='\t',nrows=3);print('REFERENCE',str(f),list(d),d.to_dict('records'))
