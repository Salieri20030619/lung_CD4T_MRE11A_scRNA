from pathlib import Path
import json
import pandas as pd,numpy as np,anndata as ad
from scipy.io import mmread
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
O=Path('/mnt/BioData/bio/wangmingyuan/lung_CD4T_MRE11A_scRNA/18_macrophage_supplement_20261019')
V=P/'11_reports/death_priority_20261018_v2';A=P/'11_reports/GSE174725_annotation_validation_20261005_v1/results'
spec=json.loads((O/'protocol.json').read_text());panels=spec['human_panels']
a=ad.read_h5ad(A/'broad_annotation_review_SERVER_ONLY.h5ad')
m=a.obs.copy();m['cell']=a.obs_names;keep=m.broad06.astype(str).isin(['0','6']);m=m.loc[keep];x=a.layers['counts'][keep.to_numpy()].tocsr();g=np.array(a.var_names)
I=P/'02_analysis/GSE174725/suppC2_20260911_164417_l173TC/input'
tm=pd.read_csv(I/'cells.tsv.gz',sep='\t');tx=mmread(I/'counts.mtx').T.tocsr();tg=pd.read_csv(I/'genes.tsv',sep='\t').symbol.to_numpy()
mods=pd.read_csv(V/'tables/module_members.tsv',sep='\t');mods=mods[mods.retained.astype(str).str.lower().eq('true')]
samples=['GSM5324866','GSM5324867','GSM5324868','GSM5324869','GSM5324870'];rows=[];values=[]
for sam in samples:
 mi=np.flatnonzero(m['sample'].eq(sam));row=dict(sample=sam,group='Exposed' if sam in samples[:2] else 'Silicosis',Mac_cells=len(mi))
 if len(mi):
  cp=np.asarray(x[mi].sum(0)).ravel()/float(m.iloc[mi].total_counts.sum())*1e4
  v=dict(sample=sam,group=row['group'],scope='BALF_AM',cells=len(mi))
  for k,genes in panels.items():v[k]=float(np.log1p(cp[np.isin(g,genes)]).mean()) if sum(np.isin(g,genes))>=3 else np.nan
  for lig in spec['ligands']:
   ii=np.flatnonzero(g==lig);v['L_'+lig]=float(np.log1p(cp[ii[0]])) if len(ii) else np.nan
  values.append(v)
 for scope in ['Core','NaiveCD4','MemCD4','Treg']:
  mask=tm['sample'].eq(sam)&(tm.core.astype(str).str.lower().eq('true') if scope=='Core' else tm.subtype.eq(scope));ii=np.flatnonzero(mask)
  row[scope+'_cells']=len(ii);row[scope+'_paired_ge10']=bool(len(mi)>=10 and len(ii)>=10)
  if not len(ii):continue
  cp=np.asarray(tx[ii].sum(0)).ravel()/float(tm.loc[mask,'library_UMI'].sum())*1e4
  v=dict(sample=sam,group=row['group'],scope=scope,cells=len(ii))
  for k,z in mods.groupby('module'):
   ix=np.isin(tg,z.gene);v[k]=float(np.log1p(cp[ix]).mean()) if ix.sum()>=3 else np.nan
  values.append(v)
 rows.append(row)
pd.DataFrame(rows).to_csv(O/'tables/silicosis_BAL_pairing.tsv',sep='\t',index=False)
pd.DataFrame(values).to_csv(O/'tables/silicosis_BAL_descriptive_scores.tsv',sep='\t',index=False)
summary=dict(status='completed_coverage_and_descriptive_only',cohort='GSE174725',disease='Silicosis',material='BALF',
 paired_ge10={sc:int(sum(r[sc+'_paired_ge10'] for r in rows)) for sc in ['Core','NaiveCD4','MemCD4','Treg']},
 no_inferential_statistics='No subtype has a matched donor with >=10 macrophages and >=10 receiver cells. Missing populations are not zero expression.',
 excluded=['CWP GSE288887/GSE227136 and pooled CWP controls'],P=None,FDR=None)
(O/'logs/silicosis_BAL_complete.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary));print(pd.DataFrame(rows).to_string(index=False))
