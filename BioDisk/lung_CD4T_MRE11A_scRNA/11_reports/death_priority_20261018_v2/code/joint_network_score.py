import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',NUMBA_NUM_THREADS='1')
from pathlib import Path
import numpy as np,pandas as pd,json,loompy
from pyscenic.cli.utils import load_signatures
from pyscenic.aucell import aucell
from ctxcore.genesig import GeneSignature
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/death_priority_20261018_v2';N=O/'joint_network'
assert json.loads((N/'network_inference_complete.json').read_text())['status']=='success'
alias={'MRE11A':'MRE11','ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','PARK2':'PRKN'}
def canon(g):return alias.get(str(g).upper(),str(g).upper())
cand=pd.read_csv(O/'tables/frozen_candidate_pairs.tsv',sep='\t');excluded=set(cand.gene)|set(pd.read_csv(O/'tables/module_members.tsv',sep='\t').gene)
all_edges=[];covers=[];sample_rows=[];runs=[];candidate_edges=[]
for arm,variant in [('All','Raw'),('BG1','BG1')]:
 d=N/'Mouse'/arm;m=pd.read_csv(d/'cells.tsv',sep='\t')
 edges=[]
 for seed in [17,29,43]:
  regs=load_signatures(str(d/f'seed{seed}/regulons.csv'));runs.append(dict(arm=arm,seed=seed,regulons=len(regs)))
  for reg in regs:
   assert '(-)' not in reg.name,'negative regulon requires separate scoring'
   for g,w in reg.gene2weight.items():edges.append(dict(arm=arm,seed=seed,tf=canon(reg.transcription_factor),target=canon(g),motif_symbol=g,weight=w))
 e=pd.DataFrame(edges);stable=e.groupby(['tf','target','motif_symbol']).agg(seed_support=('seed','nunique'),weight=('weight','mean')).reset_index();stable['arm']=arm;stable['stable_2of3']=stable.seed_support>=2;all_edges.append(stable)
 candidate_edges.append(stable[stable.target.isin(set(cand.gene))])
 with loompy.connect(str(d/'score.loom'),mode='r') as ds:
  cells=np.asarray(ds.ca.CellID).astype(str);genes=np.asarray(ds.ra.Gene).astype(str);values=ds[:,:].T
 expr=pd.DataFrame(values,index=cells,columns=genes).loc[m.cell]
 sigs=[]
 for tf,z in stable[stable.stable_2of3].groupby('tf'):
  keep=z[~z.target.isin(excluded|{tf})&z.motif_symbol.isin(genes)]
  covers.append(dict(arm=arm,variant=variant,tf=tf,stable_targets=len(z),retained_targets=len(keep),status='SCORED' if len(keep)>=5 else 'INSUFFICIENT'))
  if len(keep)>=5:sigs.append(GeneSignature(name=tf,gene2weight=dict(zip(keep.motif_symbol,keep.weight))))
 print('JOINT_AUC',arm,len(m),expr.shape[1],len(sigs),flush=True)
 a=aucell(expr,sigs,auc_threshold=.05,seed=17,num_workers=4,noweights=True,normalize=False).reindex(m.cell)
 assert np.isfinite(a.to_numpy()).all();a.to_csv(d/'clean_stable_cell_AUC.tsv.gz',sep='\t')
 for scope in ['Core','Naive','ActCD4','Treg']:
  mm=m if scope=='Core' else m[m.subtype==scope]
  for sam,z in mm.groupby('sample'):
   assert len(z)>=1
   for tf,value in a.loc[z.cell].mean().items():sample_rows.append(dict(dataset=z.dataset.iloc[0],sample=sam,scope=scope,variant=variant,group=z.group.iloc[0],day=z.day.iloc[0],cells=len(z),tf=tf,AUC=value))
pd.concat(all_edges,ignore_index=True).to_csv(O/'tables/joint_network_edge_stability.tsv.gz',sep='\t',index=False)
pd.concat(candidate_edges,ignore_index=True).to_csv(O/'tables/joint_network_candidate_edges.tsv',sep='\t',index=False)
pd.DataFrame(covers).to_csv(O/'tables/joint_network_TF_coverage.tsv',sep='\t',index=False);pd.DataFrame(runs).to_csv(O/'tables/joint_network_run_summary.tsv',sep='\t',index=False)
s=pd.DataFrame(sample_rows);fixed=[]
for (ds,var,sam,tf),z in s[s.scope!='Core'].groupby(['dataset','variant','sample','tf']):
 if len(z)==3:
  row=z.iloc[0].to_dict();row.update(scope='Core',variant='FixedComp' if var=='Raw' else 'BG1FixedComp',AUC=float(z.AUC.mean()),cells=int(z.cells.sum()));fixed.append(row)
pd.concat([s,pd.DataFrame(fixed)],ignore_index=True).to_csv(O/'tables/joint_network_TF_sample_AUC.tsv',sep='\t',index=False)
(N/'scoring_complete.json').write_text(json.dumps({'status':'success','arms':2,'runs':6,'cells_All':3724,'scored_regulons':pd.DataFrame(covers).query("status=='SCORED'").groupby('arm').size().to_dict()},indent=2),encoding='utf-8')
print('JOINT_NETWORK_SCORING_COMPLETE',flush=True)
