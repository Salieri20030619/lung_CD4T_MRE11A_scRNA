from pathlib import Path
import pandas as pd,numpy as np,itertools,json
from scipy.io import mmread
from scipy import sparse,stats
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/death_priority_20261018_v3';V=P/'11_reports/death_priority_20261018_v2'
mods=pd.read_csv(V/'tables/module_members.tsv',sep='\t');mods=mods[mods.retained.astype(str).str.lower().eq('true')]
I=P/'02_analysis/GSE174725/suppC2_20260911_164417_l173TC/input'
m=pd.read_csv(I/'cells.tsv.gz',sep='\t');genes=pd.read_csv(I/'genes.tsv',sep='\t').symbol.str.upper().to_numpy();x=mmread(I/'counts.mtx').tocsr().T.tocsr()
v=x.multiply(1e4/m.library_UMI.to_numpy()[:,None]).tocsr();v.data=np.log1p(v.data)
vals={a:np.asarray(v[:,np.isin(genes,z.gene)].mean(1)).ravel() for a,z in mods.groupby('module') if np.isin(genes,z.gene).sum()>=3}
vals.update({'G_'+a:v[:,np.flatnonzero(genes==a)[0]].toarray().ravel() for a in ['IL2RA','HIF1A','ARL6IP1','PYCARD','GCLM']})
rows=[]
for (sam,sub),idx in m[m.core.astype(str).str.lower().eq('true')].groupby(['sample','subtype']).groups.items():
 z=m.loc[idx];r=dict(sample=sam,subtype=sub,group=z.group.iloc[0],cells=len(z));r.update({a:float(b[idx].mean()) for a,b in vals.items()});rows.append(r)
s=pd.DataFrame(rows);s.to_csv(O/'tables/BAL174725_sample_scores.tsv',sep='\t',index=False);out=[]
for threshold in [1,10,20]:
 for sub,z in s[s.cells>=threshold].groupby('subtype'):
  for col in vals:
   a=z[z.group=='Silicosis'][col];b=z[z.group=='Exposed'][col]
   out.append(dict(min_cells=threshold,subtype=sub,feature=col,n_case=len(a),n_exposed=len(b),effect=float(a.mean()-b.mean()) if len(a)*len(b) else np.nan,P=np.nan,evidence='descriptive_only_sparse_subtype_coverage'))
pd.DataFrame(out).to_csv(O/'tables/BAL174725_descriptive_effects.tsv',sep='\t',index=False)
# Targeted member associations in the independent same-source human lung subset.
I=P/'02_analysis/mechanism_step2_20260927_154919_EsN4gk/input'
m=pd.read_csv(I/'human_cells.tsv',sep='\t');genes=pd.read_csv(I/'human_genes.tsv',sep='\t').symbol.str.upper().to_numpy();x=mmread(I/'human_counts.mtx').tocsr().T.tocsr();v=x.multiply(1e4/m.library.to_numpy()[:,None]).tocsr();v.data=np.log1p(v.data)
wanted=['IL2RA','HIF1A','ARL6IP1','PYCARD','GCLM','NCOR1','CYC1'];values={g:np.asarray(v[:,genes==g].sum(1)).ravel() for g in wanted};rows=[]
for (sam,sub),idx in m[(m.series=='G227')&m.core.astype(str).str.lower().eq('true')].groupby(['sample','subtype']).groups.items():
 if len(idx)<10:continue
 z=m.loc[idx];rows.append(dict(sample=sam,scope=sub,group=z.group.iloc[0],cells=len(idx),**{g:float(vv[idx].mean()) for g,vv in values.items()}))
s=pd.DataFrame(rows);s.to_csv(O/'tables/human_key_member_expression.tsv',sep='\t',index=False)
out=[]
for gene,member,sc in [('HIF1A','PYCARD','CD4-CCR7'),('ARL6IP1','GCLM','CD4-CCR7'),('IL2RA','NCOR1','Treg'),('IL2RA','CYC1','Treg')]:
 z=s[s.scope==sc].reset_index(drop=True)
 for adj in [False,True]:
  design=np.ones((len(z),1))
  if adj:design=np.column_stack([design,z.group.eq('CWP').to_numpy(float)])
  M=np.eye(len(z))-design@np.linalg.pinv(design);rx=M@stats.rankdata(z[gene]);ry=M@stats.rankdata(z[member]);den=np.linalg.norm(rx)*np.linalg.norm(ry)
  if den<1e-12:continue
  rho=rx@ry/den
  groups=[np.flatnonzero(z.group.to_numpy()==g) for g in sorted(z.group.unique())] if adj else [np.arange(len(z))]
  total=__import__('math').prod(__import__('math').factorial(len(g)) for g in groups);rng=np.random.default_rng(20261018)
  exact=total<=50000;choices=itertools.product(*[list(itertools.permutations(g)) for g in groups]) if exact else ([rng.permutation(g) for g in groups] for _ in range(19999))
  n=ge=0
  for part in choices:
   ix=np.arange(len(z))
   for g,pp in zip(groups,part):ix[g]=pp
   yp=M@stats.rankdata(z[member])[ix];d=np.linalg.norm(rx)*np.linalg.norm(yp)
   if d<1e-12:continue
   rp=rx@yp/d;ge+=abs(rp)>=abs(rho)-1e-12;n+=1
  out.append(dict(gene=gene,member=member,scope=sc,model='Condition' if adj else 'Unadjusted',n=len(z),rho=rho,P=ge/n if exact else (ge+1)/(n+1),permutations=n,exact=exact,family=f'{sc}|{adj}|targeted_members'))
o=pd.DataFrame(out);o['FDR']=o.groupby('family').P.transform(lambda x:multipletests(x,method='fdr_bh')[1]);o.to_csv(O/'tables/human_targeted_member_associations.tsv',sep='\t',index=False)
(O/'logs/BAL_member_validation_complete.txt').write_text('success',encoding='utf-8')
print('BAL_AND_MEMBER_COMPLETE',flush=True)
