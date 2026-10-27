from pathlib import Path
import pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');V=P/'11_reports/integrated_macrophage_CD4_20261015_v4';S=P/'11_reports/methods_complete_20261016_v13/supplement'
for f in [P/'02_analysis/scenic1_20260929_103513_xRoEnr/Mouse/All/regulon_edges.tsv',S/'S20_all_fixed_module_genes.tsv',S/'S20_all_fixed_module_controls.tsv.gz',S/'S20_GO_autophagy_control_members.tsv',V/'02_sender/all_ligand_sample_expression.tsv.gz']:
 a=pd.read_csv(f,sep='\t',nrows=2);print(f.name,list(a));print(a.to_dict('records'))
print((P/'02_analysis/scenic2_20261001_092655_ZKNWD6/analyze.py').read_text()[:6500])