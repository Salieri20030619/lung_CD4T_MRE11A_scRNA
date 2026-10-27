from pathlib import Path
import json,itertools
import numpy as np,pandas as pd
from scipy.stats import chi2,hypergeom,spearmanr
from statsmodels.stats.multitest import multipletests
O=Path('/mnt/BioData/bio/wangmingyuan/lung_CD4T_MRE11A_scRNA/12_scTenifoldKnk_20261018');T=O/'tables';jobs=json.loads((O/'jobs.json').read_text());checks=[]
for j in jobs:
 d=O/'runs'/j['id'];m=json.loads((d/'complete.json').read_text());assert m['status']=='success';checks.append(j['id']+' completion')
 for target in j['targets']:
  z=pd.read_csv(d/(target+'.tsv'),sep='\t')
  noise=pd.read_csv(T/'numerical_noise_audit.tsv',sep='\t');thr=float(noise.loc[(noise.id==j['id'])&(noise.target==target),'noiseLevel'].iloc[0]);expected=chi2.sf(z.FC.to_numpy(),df=1);expected[z.distance.to_numpy()<=thr]=1;assert np.allclose(expected,z['p.value'],rtol=1e-7,atol=1e-14);checks.append(j['id']+target+' model chi-square')
  q=multipletests(z['p.value'],method='fdr_bh')[1];assert np.allclose(q,z['p.adj'],rtol=1e-7,atol=1e-14);checks.append(j['id']+target+' BH')
p=pd.read_csv(T/'module_ORA.tsv',sep='\t')
for fam,z in p.groupby('family'):
 z=z[z.p.notna()];pp=hypergeom.sf(z.hits-1,z.universe,z.tested_members,z.selected);assert np.allclose(pp,z.p);assert np.allclose(multipletests(pp,method='fdr_bh')[1],z.FDR);checks.append(fam+' ORA and BH')
g=pd.read_csv(T/'all_gene_model_results.tsv.gz',sep='\t');c=pd.read_csv(T/'gene_seed_consensus.tsv.gz',sep='\t');assert len(g)==sum(j['genes']-1 for j in jobs for t in j['targets']);checks.append('row cardinality')
shared=[]
for kk,z in c[c.subtype=='Naive'].groupby(['arm','state']):
 a=z[z.target=='HIF1A'].set_index('gene');b=z[z.target=='ARL6IP1'].set_index('gene');ii=a.index.intersection(b.index);ta=set(a.nlargest(100,'median_rank').index);tb=set(b.nlargest(100,'median_rank').index)
 shared.append(dict(arm=kk[0],state=kk[1],targets='HIF1A vs ARL6IP1',shared_genes=len(ii),unsigned_rank_rho=float(spearmanr(a.loc[ii,'median_rank'],b.loc[ii,'median_rank']).statistic),top100_jaccard=len(ta&tb)/len(ta|tb)))
pd.DataFrame(shared).to_csv(T/'cross_target_specificity.tsv',sep='\t',index=False)
flags={'SFTPC','SFTPB','SFTPA1','SFTPD','SCGB1A1','SCGB3A1','SCGB3A2'};rows=[]
for kk,z in g.groupby(['arm','subtype','state','target','seed']):
 top=z.nlargest(100,'distance');hit=top[top.gene.isin(flags)];rows.append(dict(zip(['arm','subtype','state','target','seed'],kk))|dict(top100_epithelial_markers=len(hit),marker_genes=';'.join(hit.gene),top10=';'.join(top.head(10).gene)))
pd.DataFrame(rows).to_csv(T/'epithelial_marker_audit.tsv',sep='\t',index=False)
(O/'QA_complete.json').write_text(json.dumps(dict(status='success',checks=len(checks),items=checks),indent=2),encoding='utf-8');print('QA',len(checks),flush=True)