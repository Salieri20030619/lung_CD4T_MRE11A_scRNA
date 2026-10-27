from pathlib import Path
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
s=(P/'02_analysis/scenic2_20261001_092655_ZKNWD6/stats.py').read_text();print(s[:1900])
import pandas as pd
old=pd.read_csv(P/'02_analysis/mechanism_step1_20260927_152218_s4V5eF/tables/module_definition.tsv',sep='\t');new=pd.read_csv(P/'11_reports/methods_complete_20261016_v13/supplement/S20_all_fixed_module_genes.tsv',sep='\t');print('MOD_DIFF',set(old.gene)-set(new.gene),set(new.gene)-set(old.gene))