from pathlib import Path
import os,ast,json
os.environ['OPENBLAS_NUM_THREADS']='1'
import pandas as pd,numpy as np
from scipy import stats
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/human_gap_diagnostics_20261005_v1';R=O/'results';V=P/'11_reports/human_alignment_20261005_v1/results'
tree=ast.parse((O/'code/gap_compute.py').read_text());exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef)],type_ignores=[]),'<frozen functions>','exec'),globals());(O/'provenance/CONTRAST_FOLLOWUP_SPEC.txt').write_bytes((O/'code/CONTRAST_FOLLOWUP_SPEC.txt').read_bytes())
control=rd(R/'matched_control_donor_scores.tsv.gz');control=control[control.endpoint.str.startswith('MatchedControl')].groupby('sample').score.mean().rename('reference');ind=rd(R/'independence_program_donors.tsv.gz');orig=ind[ind.program=='IL7R_nominal_original'][['sample','score']].rename(columns={'score':'original'});scores=orig.merge(control,on='sample');scores['excess']=scores.original-scores.reference;wr(scores,'IL7R_fixed_reference_contrast_donors.tsv')
tf=rd(V/'TF_readout_donors.tsv.gz');tf=tf[(tf.scope=='CD4-IL7R')&(tf.tf=='IKZF2')&(tf.readout=='AUC')&(tf['filter']=='All')];rr=[]
for (variant,seed),d in tf.groupby(['variant','seed']):
 base=d[['sample','group','series','cells','mean_depth','value']].rename(columns={'cells':'receiver_cells','mean_depth':'receiver_depth','value':'x'}).merge(scores,on='sample')
 for endpoint in ['original','reference','excess']:
  q=complete(base.rename(columns={endpoint:'y'}));rr.append(dict(variant=variant,seed=seed,endpoint=endpoint,adjust='SourceReceiverDepth',**cor(q,'SourceReceiverDepth')))
 q=complete(base.rename(columns={'original':'y'}),extra='reference');rr.append(dict(variant=variant,seed=seed,endpoint='original',adjust='SourceReceiverDepthPlusReference',**cor(q,'SourceReceiverDepth',extra='reference')))
wr(rr,'IKZF2_matched_reference_contrast_models.tsv')
gene=rd(V/'alignment_gene_donors.tsv.gz');des=rd(V/'stage1_donor_design.tsv');mm=gene[(gene.key=='CD4-IL7R_All')&(gene.gene=='IKZF2')][['sample','cells','mean_depth']].rename(columns={'cells':'receiver_cells','mean_depth':'receiver_depth'});base=scores.merge(des[['sample','group','series']],on='sample').merge(mm,on='sample');chains=[]
for lig,sender in [('ICAM2','MonoMac'),('ICAM2','MacWide'),('VEGFB','MacWide')]:
 sx=gene[(gene.key==sender)&(gene.gene==lig)][['sample','cells','mean_depth','mean_CP10k']].rename(columns={'cells':'sender_cells','mean_depth':'sender_depth','mean_CP10k':'x'})
 for endpoint in ['original','reference','excess']:
  d=complete(base.merge(sx,on='sample').rename(columns={endpoint:'y'}),True);chains.append(dict(ligand=lig,sender=sender,endpoint=endpoint,**cor(d,'SourceBothDepth')))
wr(chains,'selected_nominal_chains_reference_contrast.tsv');print('CONTRAST',pd.DataFrame(rr).drop(columns=['samples']).to_string(index=False));print('CHAIN_CONTRAST',pd.DataFrame(chains).drop(columns=['samples']).to_string(index=False))
(O/'provenance/contrast_complete.json').write_text(json.dumps({'TF_diagnostics':len(rr),'selected_chain_diagnostics':len(chains),'interpretation':'Diagnostic P/CI only, not additional discovery family'},indent=2))
