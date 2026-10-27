from pathlib import Path
import os,ast,json
os.environ['OPENBLAS_NUM_THREADS']='1'
import pandas as pd,numpy as np
from scipy import stats
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/human_gap_diagnostics_20261005_v1';R=O/'results';V=P/'11_reports/human_alignment_20261005_v1/results'
tree=ast.parse((O/'code/gap_compute.py').read_text());exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef)],type_ignores=[]),'<frozen diagnostic functions>','exec'),globals())
(O/'provenance/contrast_influence_spec.txt').write_text('Final finite influence check: six existing IL7R IKZF2 AUC versions versus the already fixed original-minus-mean3controls score; leave each of the16 donors out, same model. No new endpoint/threshold/seed. Descriptive signs, ranges and model estimability, no LOO discovery claims.\n')
tf=rd(V/'TF_readout_donors.tsv.gz');tf=tf[(tf.scope=='CD4-IL7R')&(tf.tf=='IKZF2')&(tf.readout=='AUC')&(tf['filter']=='All')];pg=rd(R/'IL7R_fixed_reference_contrast_donors.tsv');rows=[]
for (variant,seed),v in tf.groupby(['variant','seed']):
 d=v[['sample','group','series','cells','mean_depth','value']].rename(columns={'cells':'receiver_cells','mean_depth':'receiver_depth','value':'x'}).merge(pg[['sample','excess']],on='sample').rename(columns={'excess':'y'});d=complete(d)
 for sample in d['sample']:rows.append(dict(variant=variant,seed=seed,omitted=sample,**cor(d[d['sample']!=sample],'SourceReceiverDepth')))
wr(rows,'excess_program_donor_influence.tsv');q=pd.DataFrame(rows).groupby(['variant','seed']).agg(deletions=('omitted','size'),computed=('P','count'),rho_min=('rho','min'),rho_max=('rho','max'),sign_flips=('rho',lambda x:sum(x>=0)),nominal=('P',lambda x:sum(x<.05)),P_max=('P','max')).reset_index();wr(q,'excess_program_influence_summary.tsv');print(q.to_string(index=False));(O/'provenance/influence_complete.json').write_text(json.dumps({'rows':len(rows),'all_six_versions_checked':True},indent=2))
