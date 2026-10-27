import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
from pathlib import Path
import anndata as ad,pandas as pd,numpy as np,json
from scipy import sparse
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/death_priority_20261018_v2');a=ad.read_h5ad(O/'objects/GSE183682_whole_clustered.h5ad');m=a.obs.copy();m['Broad']=m.cluster.astype(str).map(json.loads((O/'inputs/GSE183682_annotation.json').read_text()));x=sparse.csr_matrix(a.layers['counts']);gs=a.var_names
rows=[];count=[];frac=[];meta=[]
for sam,mm in m[m.Broad=='Mac'].groupby('sample'):
 ids=a.obs_names.get_indexer(mm.index);z=x[ids];nt=np.asarray((z[:,gs.isin(['CD3D','CD3E','TRAC'])]>0).sum(1)).ravel();nb=np.asarray((z[:,gs.isin(['MS4A1','CD79A','CD79B'])]>0).sum(1)).ravel()
 rows.append(dict(sample=sam,cells=len(ids),doublet_score_gt025=int((mm.doublet_score>.25).sum()),Tmarkers_ge2=int((nt>=2).sum()),Bmarkers_ge2=int((nb>=2).sum()),AMlike_cluster2=int((mm.cluster.astype(str)=='2').sum()),Mac_cluster4=int((mm.cluster.astype(str)=='4').sum())))
 count.append(np.asarray(z.sum(0)).ravel());frac.append(np.asarray((z>0).sum(0)).ravel()/len(ids));meta.append(dict(sample=sam,dataset='GSE183682',scope='AllMac',group=mm.group.iloc[0],day=mm.day.iloc[0],cells=len(ids),library=float(mm.library.sum())))
d=pd.DataFrame(meta);pd.DataFrame(np.array(count).T,index=pd.Index([{'MRE11A':'MRE11','ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','PARK2':'PRKN'}.get(g,g) for g in gs],name='gene'),columns=d['sample']).groupby(level=0).sum().to_csv(O/'inputs/GSE183682__AllMac.counts.tsv.gz',sep='\t');d.to_csv(O/'inputs/GSE183682__AllMac.samples.tsv',sep='\t',index=False);pd.DataFrame(rows).to_csv(O/'tables/macrophage_annotation_QC.tsv',sep='\t',index=False)
print(pd.DataFrame(rows).to_string(index=False))
