import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',NUMBA_NUM_THREADS='1')
from pathlib import Path
import json,sys,subprocess
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/death_priority_20261018_v3';V=P/'11_reports/death_priority_20261018_v2';N=O/'network_loo'
os.environ['TMPDIR']=str(O/'tmp')
import pandas as pd,numpy as np,loompy
from pyscenic.cli.utils import load_signatures
from pyscenic.aucell import aucell
from ctxcore.genesig import GeneSignature
if not (N/'complete.json').exists():
 pids=[]
 for p in Path('/proc').iterdir():
  if not p.name.isdigit():continue
  try:
   if any(a.endswith(b'/code/network_loo.py') for a in (p/'cmdline').read_bytes().split(b'\0')):pids.append(int(p.name))
  except OSError:pass
 assert len(pids)==1,pids
 subprocess.run(['tail','--pid='+str(pids[0]),'-f','/dev/null'],check=True)
assert json.loads((N/'complete.json').read_text())['status']=='success'
aliases={'MRE11A':'MRE11','ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','PARK2':'PRKN'}
def canon(g):return aliases.get(g.upper(),g.upper())
cand=pd.read_csv(V/'tables/frozen_candidate_pairs.tsv',sep='\t')
excluded=set(cand.gene)|set(pd.read_csv(V/'tables/module_members.tsv',sep='\t').gene)
routes=pd.read_csv(V/'tables/joint_network_candidate_TF_routes.tsv',sep='\t')
allrows=[];aud=[];sample=[];cov=[]
for arm,var in [('All','Raw'),('BG1','BG1')]:
 bd=V/'joint_network/Mouse'/arm
 m=pd.read_csv(bd/'cells.tsv',sep='\t')
 with loompy.connect(str(bd/'score.loom'),'r') as ds:genes=np.asarray(ds.ra.Gene).astype(str);cells=np.asarray(ds.ca.CellID).astype(str);expr=pd.DataFrame(ds[:,:].T,index=cells,columns=genes).loc[m.cell]
 base=set((canon(r.transcription_factor),canon(g)) for r in load_signatures(str(bd/'seed17/regulons.csv')) for g in r.gene2weight)
 want=set(routes.loc[routes.arm==arm,'tf'])|{'NFKB1','NFKB2','HIF1A'}
 for d in sorted((N/arm).iterdir()):
  if not d.is_dir():continue
  sam=d.name;regs=load_signatures(str(d/'regulons.csv'));edges=set();sigs=[]
  for r in regs:
   tf=canon(r.transcription_factor)
   edges.update((tf,canon(g)) for g in r.gene2weight)
   if tf in want:
    keep={g:w for g,w in r.gene2weight.items() if canon(g) not in excluded|{tf} and g in expr.columns}
    cov.append(dict(arm=arm,left_out=sam,tf=tf,raw_targets=len(r.gene2weight),clean_targets=len(keep),scored=len(keep)>=5))
    if len(keep)>=5:sigs.append(GeneSignature(tf,keep))
  aud.append(dict(arm=arm,left_out=sam,regulons=len(regs),edges=len(edges),baseline_seed17_edges=len(base),shared=len(base&edges),jaccard=len(base&edges)/len(base|edges)))
  allrows.extend(dict(arm=arm,left_out=sam,tf=tf,target=g) for tf,g in edges)
  af=d/'focal_clean_AUC.tsv.gz'
  if af.exists():a=pd.read_csv(af,sep='\t',index_col=0)
  else:a=aucell(expr,sigs,auc_threshold=.05,seed=17,num_workers=4,noweights=True,normalize=False);a.to_csv(af,sep='\t')
  assert np.isfinite(a.to_numpy()).all()
  for (sub,ss),ii in m[m['sample']!=sam].groupby(['subtype','sample']).groups.items():
   zz=m.loc[ii]
   for tf,value in a.loc[zz.cell].mean().items():sample.append(dict(arm=arm,variant=var,left_out=sam,scope=sub,sample=ss,dataset=zz.dataset.iloc[0],group=zz.group.iloc[0],day=zz.day.iloc[0],tf=tf,AUC=value,cells=len(ii)))
  print('SCORED',arm,sam,len(sigs),flush=True)
e=pd.DataFrame(allrows);e.to_csv(O/'tables/network_LOO_all_edges.tsv.gz',sep='\t',index=False)
summary=e.groupby(['arm','tf','target']).left_out.nunique().rename('folds_present').reset_index();summary['folds_total']=10
base=pd.read_csv(V/'tables/joint_network_edge_stability.tsv.gz',sep='\t')
z=base.merge(summary,on=['arm','tf','target'],how='outer');z['folds_present']=z.folds_present.fillna(0).astype(int);z['folds_total']=10;z.to_csv(O/'tables/network_LOO_edge_support.tsv.gz',sep='\t',index=False)
z[z.target.isin(set(cand.gene))].to_csv(O/'tables/network_LOO_candidate_edges.tsv',sep='\t',index=False)
pd.DataFrame(aud).to_csv(O/'tables/network_LOO_fit_audit.tsv',sep='\t',index=False);pd.DataFrame(cov).to_csv(O/'tables/network_LOO_TF_coverage.tsv',sep='\t',index=False);pd.DataFrame(sample).to_csv(O/'tables/network_LOO_sample_AUC.tsv',sep='\t',index=False)
(N/'scoring_complete.json').write_text(json.dumps({'status':'success','folds':len(aud)},indent=2),encoding='utf-8')
print('NETWORK_LOO_SCORE_COMPLETE',flush=True)
