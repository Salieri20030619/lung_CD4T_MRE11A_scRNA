import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',NUMBA_NUM_THREADS='1')
from pathlib import Path
import pandas as pd,numpy as np,loompy,json
from pyscenic.cli.utils import load_signatures
from pyscenic.aucell import aucell
from ctxcore.genesig import GeneSignature
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/death_priority_20261018_v3';V=P/'11_reports/death_priority_20261018_v2'
alias={'MRE11A':'MRE11','ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','PARK2':'PRKN'}
def canon(g):return alias.get(g.upper(),g.upper())
cand=pd.read_csv(V/'tables/frozen_candidate_pairs.tsv',sep='\t');excluded=set(cand.gene)|set(pd.read_csv(V/'tables/module_members.tsv',sep='\t').gene)
routes=pd.read_csv(V/'tables/joint_network_candidate_TF_routes.tsv',sep='\t')
rows=[];cov=[]
for arm,var in [('All','Raw'),('BG1','BG1')]:
 d=V/'joint_network/Mouse'/arm;m=pd.read_csv(d/'cells.tsv',sep='\t')
 with loompy.connect(str(d/'score.loom'),'r') as ds:expr=pd.DataFrame(ds[:,:].T,index=np.asarray(ds.ca.CellID).astype(str),columns=np.asarray(ds.ra.Gene).astype(str)).loc[m.cell]
 sigs=[];want=set(routes.loc[routes.arm==arm,'tf'])|{'NFKB1','NFKB2','HIF1A'}
 for r in load_signatures(str(d/'seed17/regulons.csv')):
  tf=canon(r.transcription_factor)
  if tf not in want:continue
  keep={g:w for g,w in r.gene2weight.items() if canon(g) not in excluded|{tf} and g in expr.columns}
  cov.append(dict(arm=arm,tf=tf,clean_targets=len(keep),scored=len(keep)>=5))
  if len(keep)>=5:sigs.append(GeneSignature(tf,keep))
 a=aucell(expr,sigs,auc_threshold=.05,seed=17,num_workers=4,noweights=True,normalize=False)
 assert np.isfinite(a.to_numpy()).all()
 a.to_csv(O/'network_loo'/('full_seed17_'+arm+'_cell_AUC.tsv.gz'),sep='\t')
 for (sc,sam),ii in m.groupby(['subtype','sample']).groups.items():
  z=m.loc[ii]
  for tf,val in a.loc[z.cell].mean().items():rows.append(dict(arm=arm,variant=var,scope=sc,sample=sam,dataset=z.dataset.iloc[0],day=z.day.iloc[0],group=z.group.iloc[0],tf=tf,AUC=val,cells=len(z)))
pd.DataFrame(rows).to_csv(O/'tables/network_seed17_reference_AUC.tsv',sep='\t',index=False);pd.DataFrame(cov).to_csv(O/'tables/network_seed17_reference_TF_coverage.tsv',sep='\t',index=False)
print('FULL_SEED17_REFERENCE_SCORED',flush=True)
