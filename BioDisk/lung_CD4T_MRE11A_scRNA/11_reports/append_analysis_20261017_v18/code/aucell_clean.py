import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMBA_NUM_THREADS='1')
from pathlib import Path
import json,pandas as pd,numpy as np,loompy
from pyscenic.cli.utils import load_signatures
from pyscenic.aucell import aucell
from ctxcore.genesig import GeneSignature
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/append_analysis_20261017_v18';S=P/'02_analysis/scenic1_20260929_103513_xRoEnr/Mouse/All';B=P/'11_reports/three_cohort_four_axis_20261017_v16';plan=json.loads((O/'protocol_frozen.json').read_text());cand=set(g for g,sc in plan['focus_pairs']);mods=pd.read_csv(O/'inputs/fixed_modules.tsv',sep='\t');gly=pd.read_csv(O/'inputs/glycolysis_GO_members.tsv',sep='\t');excluded=cand|set(mods.gene)|set(gly.gene)|{'PENK','TNFRSF9','LAG3','MRE11'}
tfs={'FOXP1','TRP53','NFKB1','NFKB2','RELA','RELB','JUN','JUNB','FOS','FOSB','NFATC1','NFATC2'}
alias={'MRE11A':'MRE11','ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','PARK2':'PRKN'}
def canon(g):g=str(g).upper();return alias.get(g,g)
regs=load_signatures(str(S/'seed17/regulons.csv'));sigs=[];members=[];cov=[]
for r in regs:
 tf=canon(r.transcription_factor)
 if tf not in tfs:continue
 kept={g:w for g,w in r.gene2weight.items() if canon(g) not in excluded|{tf}}
 cov.append(dict(tf=tf,regulon=r.name,original_genes=len(r.gene2weight),retained_genes=len(kept),status='SCORED' if kept else 'EMPTY'))
 for g,w in r.gene2weight.items():members.append(dict(tf=tf,target=g,canonical=canon(g),weight=w,retained=g in kept))
 if kept:sigs.append(GeneSignature(name=r.name,gene2weight=kept))
for tf in sorted(tfs-set(r['tf'] for r in cov)):cov.append(dict(tf=tf,status='REGULON_ABSENT'))
pd.DataFrame(members).to_csv(O/'tables/clean_regulon_members.tsv',sep='\t',index=False);pd.DataFrame(cov).to_csv(O/'tables/clean_regulon_coverage.tsv',sep='\t',index=False)
with loompy.connect(str(S/'score.loom'),mode='r') as ds:
 cells=np.asarray(ds.ca.CellID).astype(str);genes=np.asarray(ds.ra.Gene).astype(str);values=ds[:,:].T
expression=pd.DataFrame(values,index=cells,columns=genes);m=pd.read_csv(B/'inputs/GSE314943__Core.cells.tsv.gz',sep='\t');assert set(cells)==set(m.cell);expression=expression.loc[m.cell]
print('AUC_START',expression.shape,len(sigs),flush=True);a=aucell(expression,sigs,auc_threshold=.05,seed=17,num_workers=2,noweights=True,normalize=False);assert np.isfinite(a.to_numpy()).all();a.columns=[canon(str(c).split('(')[0]) for c in a.columns];a=a.reindex(m.cell);(O/'objects').mkdir(exist_ok=True);a.to_csv(O/'objects/clean_AUC_cells.tsv.gz',sep='\t');rows=[]
for var in ['Raw','BG1']:
 for scope in ['Core','Naive','ActCD4','Treg']:
  take=np.ones(len(m),bool) if scope=='Core' else m.subtype.eq(scope).to_numpy()
  if var=='BG1':take&=m.background_fraction.to_numpy()<=.01
  for sample,mm in m.loc[take].groupby('sample'):
   vals=a.loc[mm.cell].mean()
   for tf,v in vals.items():rows.append(dict(variant=var,scope=scope,sample=sample,group=mm.group.iloc[0],cells=len(mm),tf=tf,AUC=v))
d=pd.DataFrame(rows);fix=[]
for (var,sample,tf),z in d[d.scope.isin(['Naive','ActCD4','Treg'])].groupby(['variant','sample','tf']):
 if len(z)!=3:continue
 r=z.iloc[0].to_dict();r.update(variant='FixedComp' if var=='Raw' else 'BG1FixedComp',scope='Core',AUC=float(z.AUC.mean()),cells=int(z.cells.sum()));fix.append(r)
pd.concat([d,pd.DataFrame(fix)]).to_csv(O/'tables/clean_TF_sample_AUC.tsv',sep='\t',index=False);(O/'logs/AUC_complete.txt').write_text('success',encoding='utf-8');print('AUC_COMPLETE',len(d)+len(fix),flush=True)