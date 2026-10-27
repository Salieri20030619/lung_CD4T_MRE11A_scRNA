import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2',NUMBA_NUM_THREADS='2')
from pathlib import Path
import pandas as pd,numpy as np,json
from scipy import sparse
from pyscenic.cli.utils import load_signatures
from pyscenic.aucell import aucell
from ctxcore.genesig import GeneSignature
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/death_priority_20261018_v2';B=P/'11_reports/three_cohort_four_axis_20261017_v16'
alias={'MRE11A':'MRE11','ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','PARK2':'PRKN'}
def canon(g):return alias.get(str(g).upper(),str(g).upper())
inputs=[]
for ds in ['GSE314943','GSE183682']:
 stem=ds+'__Core';x=sparse.load_npz(B/'inputs'/(stem+'.cellcounts.npz'));m=pd.read_csv(B/'inputs'/(stem+'.cells.tsv.gz'),sep='\t');m['dataset']=ds;m['uid']=ds+'|'+m.cell;g=pd.read_csv(B/'inputs'/(stem+'.genes.txt'),header=None)[0].to_numpy();inputs.append((x,m,g))
common=sorted(set(inputs[0][2])&set(inputs[1][2]));mat=np.vstack([x[:,pd.Index(g).get_indexer(common)].toarray().astype('float32') for x,m,g in inputs]);m=pd.concat([m for x,m,g in inputs],ignore_index=True)
excl=set(pd.read_csv(O/'tables/module_members.tsv',sep='\t').gene)|set(pd.read_csv(O/'tables/frozen_candidate_pairs.tsv',sep='\t').gene)
regs=load_signatures(str(P/'02_analysis/scenic1_20260929_103513_xRoEnr/Mouse/All/seed17/regulons.csv'))
sigs=[];members=[];cov=[]
for r in regs:
 tf=canon(r.transcription_factor);mapped={canon(g):w for g,w in r.gene2weight.items() if canon(g) in common};kept={g:w for g,w in mapped.items() if g not in excl|{tf}}
 cov.append(dict(tf=tf,regulon=r.name,original=len(r.gene2weight),mapped=len(mapped),retained=len(kept),status='SCORED' if len(kept)>=5 else 'INSUFFICIENT'))
 for g,w in mapped.items():members.append(dict(tf=tf,regulon=r.name,gene=g,weight=w,retained=g in kept))
 if len(kept)>=5:sigs.append(GeneSignature(name=r.name,gene2weight=kept))
pd.DataFrame(members).to_csv(O/'tables/TF_regulon_members.tsv',sep='\t',index=False);pd.DataFrame(cov).to_csv(O/'tables/TF_coverage.tsv',sep='\t',index=False)
print('AUC',mat.shape,len(sigs),flush=True)
a=aucell(pd.DataFrame(mat,index=m.uid,columns=common),sigs,auc_threshold=.05,seed=17,num_workers=2,noweights=True,normalize=False)
a.columns=[canon(str(k).split('(')[0]) for k in a.columns];assert not a.columns.duplicated().any();a=a.reindex(m.uid);a.to_csv(O/'objects/TF_cell_AUC.tsv.gz',sep='\t')
rows=[]
for var in ['Raw','BG1']:
 for scope in ['Core','Naive','ActCD4','Treg']:
  take=np.ones(len(m),bool) if scope=='Core' else m.subtype.eq(scope).to_numpy()
  if var=='BG1':take&=m.background_fraction.to_numpy()<=.01
  for (ds,sam),mm in m.loc[take].groupby(['dataset','sample']):
   for tf,v in a.loc[mm.uid].mean().items():rows.append(dict(dataset=ds,sample=sam,scope=scope,variant=var,group=mm.group.iloc[0],day=mm.day.iloc[0],cells=len(mm),tf=tf,AUC=v))
d=pd.DataFrame(rows);fixed=[]
for (ds,var,sam,tf),z in d[d.scope!='Core'].groupby(['dataset','variant','sample','tf']):
 if len(z)==3:
  r=z.iloc[0].to_dict();r.update(scope='Core',variant='FixedComp' if var=='Raw' else 'BG1FixedComp',AUC=float(z.AUC.mean()),cells=int(z.cells.sum()));fixed.append(r)
pd.concat([d,pd.DataFrame(fixed)],ignore_index=True).to_csv(O/'tables/TF_sample_AUC.tsv',sep='\t',index=False)
(O/'logs/TF_scoring_complete.json').write_text(json.dumps({'status':'success','cells':len(m),'n':10,'regulons':len(sigs)}),encoding='utf-8')
print('TF_SCORED',len(sigs),flush=True)
