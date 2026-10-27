from pathlib import Path
import os,json,ast
os.environ['OPENBLAS_NUM_THREADS']='1'
import pandas as pd,numpy as np
from scipy import stats
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/human_gap_diagnostics_20261005_v1';R=O/'results';V=P/'11_reports/human_alignment_20261005_v1/results';H=P/'11_reports/human_first_20261004_v1/results';A=P/'02_analysis'
s=(O/'code/gap_compute.py').read_text();tree=ast.parse(s);exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef)],type_ignores=[]),'<frozen functions>','exec'),globals())
up=rd(R/'all_candidate_TF_stepwise_models.tsv.gz');pr=rd(R/'TF_program_stepwise_models.tsv');chains=rd(R/'nominal_upstream_same_roster_program_checks.tsv');flow=rd(R/'IL15_fixed_condition_attrition.tsv');tf=rd(V/'TF_readout_donors.tsv.gz');tf=tf[tf['filter']=='All'];meas=rd(R/'raw_donor_gene_measurements.tsv.gz')
exec(compile(s[s.index('# Independent program recomputations'):],'<resume corrected original gene space>','exec'),globals())
