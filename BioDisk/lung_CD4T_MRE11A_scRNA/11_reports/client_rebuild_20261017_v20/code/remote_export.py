from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
import h5py,numpy as np,pandas as pd,scipy.sparse as sp,json
from anndata.io import read_elem
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/client_rebuild_20261017_v20';T=O/'tables'
paths={'atlas':'stage3_20260905_214228_VPtd0r/objects/whole_annotated.h5ad','cd4':'stage5_20260905_222457_a2dApT/objects/CD4_core.h5ad','trajectory':'stage9_20260906_080538_wN3sjO/objects/K15_trajectory.h5ad'}
markers=['CD3D','CD3E','TRAC','CD4','CD8A','MS4A1','CD79A','NKG7','LYZ2','CSF1R','S100A8','FCER1A','CST3','SFTPC','PECAM1','COL1A1','FOXP3','IL2RA','CTLA4','TCF7','CCR7','IL1R2','TTC39C','SKIL']
for label,rel in paths.items():
 with h5py.File(P/'02_analysis/GSE314943'/rel,'r') as h:
  m=read_elem(h['obs']);v=read_elem(h['var']);m['cell']=m.index;em='X_umap_Harmony' if label!='atlas' else 'X_umap';u=h['obsm'][em][:];m['UMAP1']=u[:,0];m['UMAP2']=u[:,1]
  m.to_csv(T/(label+'_cells.tsv.gz'),sep='\t',index=False)
  if label=='trajectory':continue
  sy=v.symbol.astype(str).str.upper();cols=np.flatnonzero(sy.isin(markers));ss=sy.iloc[cols].to_numpy();mat=h['layers']['counts'];ptr=mat['indptr'][:];out=[]
  assert mat.attrs['encoding-type']=='csr_matrix'
  for start in range(0,len(m),2000):
   end=min(start+2000,len(m));a,b=int(ptr[start]),int(ptr[end]);xx=sp.csr_matrix((mat['data'][a:b],mat['indices'][a:b],ptr[start:end+1]-a),shape=(end-start,len(v)));out.append(xx[:,cols].toarray())
  xx=np.vstack(out);lib=m.total_counts.to_numpy();scope='Broad' if label=='atlas' else 'Subtype';rows=[]
  for (sample,state),ids in m.groupby(['sample',scope],observed=True).indices.items():
   for gene in sorted(set(ss)):
    val=xx[ids][:,ss==gene].sum(1);rows.append(dict(sample=sample,group=m.iloc[ids[0]]['group'],state=state,gene=gene,cells=len(ids),positive_cells=int((val>0).sum()),fraction=float((val>0).mean()),mean_logCP10k=float(np.log1p(val/lib[ids]*1e4).mean()),CP10k=float(val.sum()/lib[ids].sum()*1e4)))
  pd.DataFrame(rows).to_csv(T/(label+'_markers.tsv'),sep='\t',index=False)
  print(label,len(m),m[scope].value_counts().to_dict(),flush=True)
# Within T-cell identity for focus positives, using frozen canonical count matrix.
B=P/'11_reports/three_cohort_four_axis_20261017_v16';x=sp.load_npz(B/'inputs/GSE314943__Core.cellcounts.npz');m=pd.read_csv(B/'inputs/GSE314943__Core.cells.tsv.gz',sep='\t');g=pd.read_csv(B/'inputs/GSE314943__Core.genes.txt',header=None)[0];ix={s:i for i,s in enumerate(g)}
tmark=np.asarray(x[:,[ix[s] for s in ['CD3D','CD3E','TRAC'] if s in ix]].sum(1)).ravel()>0;rows=[]
for gene,sc in [('IL1R2','Treg'),('SKIL','Treg'),('TTC39C','ActCD4')]:
 val=x[:,ix[gene]].toarray().ravel()
 for sample,mm in m[m.subtype==sc].groupby('sample'):
  ids=mm.index.to_numpy();pos=ids[val[ids]>0];rows.append(dict(gene=gene,scope=sc,sample=sample,group=mm.group.iloc[0],cells=len(ids),positive=len(pos),positive_Tmarker=int(tmark[pos].sum()),positive_BG1=int((m.loc[pos,'background_fraction']<=.01).sum()),positive_median_bg=float(m.loc[pos,'background_fraction'].median()) if len(pos) else np.nan))
pd.DataFrame(rows).to_csv(T/'candidate_positive_identity.tsv',sep='\t',index=False)
for n in ['KEGG_links.tsv','KEGG_names.tsv']:
 src=P/'02_analysis/GSE314943/stage7_20260905_234200_aiQYL0/database'/n;print(n,src.read_text().splitlines()[:3])
(O/'export_complete.json').write_text(json.dumps({'status':'success','sources':paths}),encoding='utf-8');print('EXPORT_SUCCESS',flush=True)