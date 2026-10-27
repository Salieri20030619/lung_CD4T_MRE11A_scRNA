import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
from pathlib import Path
import pandas as pd,numpy as np,json,itertools
from scipy import sparse,stats
from scipy.io import mmread
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/death_priority_20261018_v3';V=P/'11_reports/death_priority_20261018_v2';H=P/'02_analysis/mechanism_step2_20260927_154919_EsN4gk/input'
def wr(d,n):pd.DataFrame(d).to_csv(O/'tables'/n,sep='\t',index=False)
mods=pd.read_csv(V/'tables/module_members.tsv',sep='\t');mods=mods[mods.retained.astype(str).str.lower().eq('true')]
m=pd.read_csv(H/'human_cells.tsv',sep='\t');g=pd.read_csv(H/'human_genes.tsv',sep='\t');print('GENES',g.columns.tolist(),flush=True)
x=mmread(H/'human_counts.mtx').tocsr().T.tocsr();genes=g.symbol.str.upper().to_numpy();assert x.shape==(len(m),len(g))
aliases={'MRE11A':'MRE11','ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','PARK2':'PRKN'}
genes=np.array([aliases.get(a,a) for a in genes]);u,inv=np.unique(genes,return_inverse=True);x=x@sparse.csr_matrix((np.ones(len(g)),(np.arange(len(g)),inv)),shape=(len(g),len(u)));genes=u
core=m.core.astype(str).str.lower().eq('true').to_numpy();m=m[core].reset_index(drop=True);x=x[core];assert np.all(np.asarray(x.sum(1)).ravel()<=m.library.to_numpy()+1e-6)
v=x.multiply(1e4/m.library.to_numpy()[:,None]).tocsr();v.data=np.log1p(v.data)
members={a:sorted(set(z.gene)&set(genes)) for a,z in mods.groupby('module')};wr([dict(module=a,gene=b,origin='frozen_mouse_program_conserved_symbol_projection') for a,bg in members.items() for b in bg],'human_module_members.tsv')
cand=pd.read_csv(V/'tables/frozen_candidate_pairs.tsv',sep='\t');cs=sorted(set(cand.gene)|{'PYCARD','GCLM','NFKB2'})
vals={a:np.asarray(v[:,np.isin(genes,b)].mean(1)).ravel() for a,b in members.items() if len(b)>=3}
vals.update({'G_'+a:v[:,np.flatnonzero(genes==a)[0]].toarray().ravel() for a in cs if a in genes})
rows=[];cover=[];counts=[];meta=[]
for (sam,sub),idx in m.groupby(['sample','subtype']).groups.items():
 z=m.loc[idx];row=dict(sample=sam,subtype=sub,group=z.group.iloc[0],series=z.series.iloc[0],cells=len(z));cover.append(row)
 if len(z)<10:continue
 r=row.copy();r.update({a:float(b[idx].mean()) for a,b in vals.items()});rows.append(r)
 if r['series']=='G227' and sub in ['Treg','CD4-CCR7','CD4-GZMA']:
  key=sam+'__'+sub;counts.append(np.asarray(x[idx].sum(0)).ravel());meta.append(dict(row,column=key,library=float(z.library.sum())))
wr(cover,'human_coverage.tsv');s=pd.DataFrame(rows);wr(s,'human_sample_scores.tsv')
pd.DataFrame(np.array(counts).T,index=genes,columns=[r['column'] for r in meta]).to_csv(O/'tables/human_G227_counts.tsv.gz',sep='\t');wr(meta,'human_G227_samples.tsv')
def corr(z,a,b,condition):
 design=pd.get_dummies(z['series'],drop_first=False,dtype=float).to_numpy()
 if condition:design=np.column_stack([design,z.group.eq('CWP').to_numpy(float)])
 M=np.eye(len(z))-design@np.linalg.pinv(design);xx=M@stats.rankdata(z[a]);yy=M@stats.rankdata(z[b])
 if np.linalg.norm(xx)<1e-10 or np.linalg.norm(yy)<1e-10:return np.nan,np.nan,0
 r=xx@yy/(np.linalg.norm(xx)*np.linalg.norm(yy))
 keys=z.series.astype(str)+(('|'+z.group.astype(str)) if condition else '')
 blocks=[np.flatnonzero(keys.to_numpy()==v) for v in sorted(keys.unique())]
 nperm=int(np.prod([float(__import__('math').factorial(len(b))) for b in blocks]))
 rng=np.random.default_rng(20261018)
 if nperm<=50000:perms=itertools.product(*[list(itertools.permutations(b)) for b in blocks]);exact=True
 else:perms=([rng.permutation(b) for b in blocks] for _ in range(19999));exact=False
 ge=0;nn=0
 for pp in perms:
  ix=np.arange(len(z))
  for bi,pi in zip(blocks,pp):ix[bi]=pi
  yp=M@yy[ix];den=np.linalg.norm(xx)*np.linalg.norm(yp);rp=0 if den<1e-10 else xx@yp/den
  ge+=abs(rp)>=abs(r)-1e-12;nn+=1
 return r,(ge/nn if exact else (ge+1)/(nn+1)),nn
out=[]
for cohort in ['AllHuman','G227']:
 for gene,sub in [('IL2RA','Treg'),('HIF1A','CD4-CCR7'),('ARL6IP1','CD4-CCR7')]:
  for threshold in [10,20]:
   z=s[(s.subtype==sub)&(s.cells>=threshold)&((s.series=='G227') if cohort=='G227' else True)].reset_index(drop=True)
   if len(z)<5:continue
   for adjust in [False,True]:
    for mod in members:
     if mod not in z:continue
     r,p,np_=corr(z,'G_'+gene,mod,adjust);out.append(dict(cohort=cohort,gene=gene,scope=sub,min_cells=threshold,model='SourceCondition' if adjust else 'Source',program=mod,n=len(z),rho=r,P=p,permutations=np_,family=f'{cohort}|{gene}|{threshold}|{adjust}|all14modules'))
o=pd.DataFrame(out);o['FDR']=np.nan
for _,idx in o.groupby('family').groups.items():
 valid=o.loc[idx,'P'].notna();ii=np.array(list(idx))[valid];o.loc[ii,'FDR']=multipletests(o.loc[ii,'P'],method='fdr_bh')[1]
wr(o,'human_candidate_associations.tsv')
# Original joint cohort directions remain internal consistency, not independent validation.
ss=pd.read_csv(V/'tables/sample_scores.tsv',sep='\t');out=[]
for var in ['Raw','BG1']:
 for sc in ['Naive','ActCD4','Treg']:
  z=ss[(ss.variant==var)&(ss.scope==sc)]
  for (ds,day),q in z.groupby(['dataset','day']):
   if q.group.nunique()!=2:continue
   for col in [c for c in q if c.startswith(('G_','M_'))]:
    eff=q.loc[q.group=='Silica',col].mean()-q.loc[q.group=='Saline',col].mean()
    out.append(dict(dataset=ds,day=day,scope=sc,variant=var,feature=col,n_control=int(q.group.eq('Saline').sum()),n_case=int(q.group.eq('Silica').sum()),effect=eff,evidence='descriptive_direction_no_new_P'))
wr(out,'cross_cohort_direction.tsv')
(O/'logs/human_validation_complete.txt').write_text('success',encoding='utf-8')
print('HUMAN_EXTERNAL_COMPLETE',len(o),flush=True)
