from pathlib import Path
import pandas as pd,numpy as np,scipy.sparse as sp,json
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports');B=P/'three_cohort_four_axis_20261017_v16';O=P/'dual_benchmark_20261017_v17'
genes=json.loads((O/'genes.json').read_text());rows=[]
for ds,core in [('GSE314943','Core'),('GSE183682','Core'),('GSE274372','TLike')]:
 stem=ds+'__'+core;m=pd.read_csv(B/'inputs'/f'{stem}.cells.tsv.gz',sep='\t');g=pd.read_csv(B/'inputs'/f'{stem}.genes.txt',header=None)[0].tolist();idx={v:i for i,v in enumerate(g)};x=sp.load_npz(B/'inputs'/f'{stem}.cellcounts.npz');tpos=np.asarray(x[:,[idx[v] for v in ['CD3D','CD3E','CD3G','TRAC','CD247'] if v in idx]].sum(1)).ravel()>0
 for scope in ([core,'Naive','ActCD4','Treg'] if core=='Core' else [core]):
  take=np.ones(len(m),bool) if scope==core else m.subtype.eq(scope).to_numpy()
  for sample,z in m.loc[take].groupby('sample'):
   ids=z.index.to_numpy();lib=z.library.to_numpy();bg=z.background_fraction.to_numpy();clean=bg<=.01
   for gene in genes:
    present=gene in idx;v=x[ids,idx[gene]].toarray().ravel() if present else np.zeros(len(ids));pos=v>0
    rows.append(dict(gene=gene,dataset=ds,scope=scope,sample=sample,group=z.group.iloc[0],cells=len(ids),feature_present=present,UMI=int(v.sum()),positive=int(pos.sum()),multiUMI=int((v>=2).sum()),Tmarker_positive=int((pos&tpos[ids]).sum()),CP10k=float(v.sum()/lib.sum()*1e4),BG1_cells=int(clean.sum()),BG1_CP10k=float(v[clean].sum()/lib[clean].sum()*1e4) if clean.any() else None))
pd.DataFrame(rows).to_csv(O/'tables/90_gene_sample_expression.tsv',sep='\t',index=False)
print('EXPRESSION',len(rows))