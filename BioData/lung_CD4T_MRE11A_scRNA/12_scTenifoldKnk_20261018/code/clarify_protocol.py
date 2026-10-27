from pathlib import Path
import json
O=Path('/mnt/BioData/bio/wangmingyuan/lung_CD4T_MRE11A_scRNA/12_scTenifoldKnk_20261018')
p=json.loads((O/'protocol.json').read_text());p['sampling']='each seed caps donor cells at150 without replacement; package internal sampling WITH replacement, min(500,floor(0.8*available_cells)) draws per network,10 networks';p['gene_FDR']='package BH across full modeled gene set including knockout target; target excluded from enrichment and interpretation';p['module_members']='retained==True in frozen v2 module_members.tsv; target excluded for each knockout';p['secondary_QC']='cross-target response similarity and surfactant/secretory epithelial transcript audit'
(O/'protocol.json').write_text(json.dumps(p,indent=2),encoding='utf-8')