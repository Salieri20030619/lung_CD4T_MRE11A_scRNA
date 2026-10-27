import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
from pathlib import Path
import sys,json,subprocess
import pandas as pd,numpy as np
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');V=P/'11_reports/death_priority_20261018_v2';O=P/'11_reports/death_priority_20261018_v3'
if not (O/'network_loo/scoring_complete.json').exists():
 pids=[]
 for pp in Path('/proc').iterdir():
  if pp.name.isdigit():
   try:
    args=(pp/'cmdline').read_bytes().split(b'\0')
    if any(a.endswith(b'/code/network_score.py') for a in args):pids.append(int(pp.name))
   except OSError:pass
 assert len(pids)==1,pids
 subprocess.run(['tail','--pid='+str(pids[0]),'-f','/dev/null'],check=True)
assert json.loads((O/'network_loo/scoring_complete.json').read_text())['status']=='success'
sys.path.insert(0,str(V/'code'));import statslib as S
S.O=O
scores=pd.read_csv(V/'tables/sample_scores.tsv',sep='\t');a=pd.read_csv(O/'tables/network_LOO_sample_AUC.tsv',sep='\t');routes=pd.read_csv(V/'tables/joint_network_candidate_TF_routes.tsv',sep='\t')
for (arm,left,scope),aa in a.groupby(['arm','left_out','scope']):
 var='Raw' if arm=='All' else 'BG1'
 z=scores[(scores.variant==var)&(scores.scope==scope)&scores.dataset.ne('GSE274372')&scores['sample'].ne(left)].merge(aa.pivot(index=['dataset','sample'],columns='tf',values='AUC').add_prefix('TF_'),on=['dataset','sample']).sort_values(['dataset','sample']).reset_index(drop=True)
 assert len(z)==9
 for model in ['JointStratum','JointStratumCondition']:
  for r in routes[(routes.arm==arm)&(routes.scope==scope)].itertuples():
   S.assoc(z,'TF_'+r.tf,'G_'+r.gene,'LOO_'+left,var,'candidate_'+r.branch,scope,model,extra={'left_out':left,'arm':arm,'gene':r.gene,'tf':r.tf,'branch':r.branch})
out=S.finalize(S.results,'network_LOO_candidate_associations.tsv')
base=pd.read_csv(V/'tables/joint_network_TF_associations.tsv',sep='\t')
rows=[]
for keys,z in out.groupby(['arm','variant','scope','model','gene','tf','branch']):
 arm,var,scope,model,gene,tf,branch=keys
 b=base[(base.variant==var)&(base.scope==scope)&(base.model==model)&(base.gene==gene)&(base.tf==tf)&base.kind.str.contains('candidate_')]
 assert len(b)==1
 r=b.iloc[0];valid=z[z.status=='OK']
 rows.append(dict(zip(['arm','variant','scope','model','gene','tf','branch'],keys),folds_total=10,folds_evaluable=len(valid),same_direction=int((np.sign(valid.rho)==np.sign(r.rho)).sum()),nominal=int((valid.P<.05).sum()),FDR05=int((valid.FDR<.05).sum()),rho_min=valid.rho.min(),rho_max=valid.rho.max(),baseline_rho=r.rho,baseline_P=r.P,baseline_FDR=r.FDR))
pd.DataFrame(rows).to_csv(O/'tables/network_LOO_association_summary.tsv',sep='\t',index=False)
# Transparent suitability decision for the optional in-silico intervention.
cov=pd.read_csv(V/'tables/joint_network_TF_coverage.tsv',sep='\t');e=pd.read_csv(O/'tables/network_LOO_edge_support.tsv.gz',sep='\t')
elig=[]
for gene in ['HIF1A','IL2RA','ARL6IP1']:
 z=cov[cov.tf==gene];ev=e[(e.tf==gene)&e.stable_2of3.eq(True)]
 elig.append(dict(gene=gene,All_clean_targets=int(z.loc[z.arm=='All','retained_targets'].sum()),BG1_clean_targets=int(z.loc[z.arm=='BG1','retained_targets'].sum()),All_outgoing_edges_LOO_ge8=int(((ev.arm=='All')&(ev.folds_present>=8)).sum()),decision='NOT_RUN_MODEL_NOT_SUPPORTED',reason='TF-only SCENIC is not a receptor/ER protein perturbation model' if gene!='HIF1A' else 'Only 6 cleaned All targets and no BG1 scored regulon; no cross-scope stable intervention model. Actual mouse and human HIF1A perturbation evidence analyzed instead.'))
pd.DataFrame(elig).to_csv(O/'tables/computational_perturbation_suitability.tsv',sep='\t',index=False)
(O/'network_loo/statistics_complete.json').write_text(json.dumps({'status':'success','tests':len(out),'summaries':len(rows)},indent=2),encoding='utf-8')
print('NETWORK_LOO_STATISTICS_COMPLETE',flush=True)
