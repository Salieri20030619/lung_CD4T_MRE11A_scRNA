from pathlib import Path
import os,sys,csv,gzip,io,tarfile,json,traceback
import numpy as np
import pandas as pd
from scipy import sparse,io as sio
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(os.environ['PROJECT_ROOT']);OUT=Path(sys.argv[1]);OUT.mkdir(parents=True,exist_ok=True)
PANEL='CD3D CD3E TRAC CD4 CD8A CD8B CD8B1 IL7R TCF7 LEF1 CCR7 FOXP3 IL2RA CTLA4 ICOS NKG7 GNLY NCR1 KLRD1 GZMB CCL5 LYZ LYZ2 LST1 TYROBP EPCAM SFTPC MS4A1 CD79A MRE11 MRE11A'.split()
BROAD='CD3D CD3E CD3G TRAC TRBC1 TRBC2 NKG7 GNLY NCR1 KLRD1'.split()
status=[]
def collapse(X,genes,cells,sample):
 names=np.array([str(g).strip().upper() for g in genes]);unique,idx=np.unique(names,return_inverse=True)
 mapping=sparse.csr_matrix((np.ones(len(idx)),(np.arange(len(idx)),idx)),shape=(len(idx),len(unique)))
 X=(sparse.csr_matrix(X)@mapping).tocsr().astype('float32')
 assert X.shape[0]==len(cells) and np.isfinite(X.data).all() and (X.data>=0).all()
 a=ad.AnnData(X,obs=pd.DataFrame({'sample':sample,'barcode':list(cells)},index=[sample+':'+str(c) for c in cells]),var=pd.DataFrame(index=unique))
 assert a.obs_names.is_unique
 return a

def load174():
 candidates=[p for p in (ROOT/'02_analysis').glob('screen6_G174_recover_*') if (p/'recovery_audit.json').exists() and len(list(p.glob('GSM*_cells.tsv.gz')))==5]
 assert candidates,'Completed G174 recovery not found'
 source=max(candidates,key=lambda p:p.stat().st_mtime);print('G174_METADATA='+str(source),flush=True)
 objects=[]
 with tarfile.open(ROOT/'00_raw/GSE174725/GSE174725_RAW.tar') as t:
  members=[m for m in t if m.isfile() and m.name.endswith('.csv.gz')];assert len(members)==5
  for member in members:
   sample=Path(member.name).name.split('_')[0];meta=pd.read_csv(source/(sample+'_cells.tsv.gz'),sep='\t')
   cols=[k for k in BROAD if k in meta];keep=meta[cols].gt(0).any(axis=1);chosen=set(meta.loc[keep,'cell'])
   with t.extractfile(member) as b,gzip.GzipFile(fileobj=b) as z,io.TextIOWrapper(z) as h:
    reader=csv.reader(h);header=next(reader);ii=[j for j,c in enumerate(header[1:],1) if c in chosen];cells=[header[j] for j in ii]
    assert len(cells)==len(chosen) and len(cells)>0
    data=[];indices=[];indptr=[0];genes=[]
    for row in reader:
     if not row or not any(v.strip() for v in row):continue
     assert len(row)==len(header)
     v=np.array([row[j] for j in ii],dtype='float32');assert np.isfinite(v).all()
     nz=np.flatnonzero(v);indices.extend(nz);data.extend(v[nz]);indptr.append(len(data));genes.append(row[0])
    X=sparse.csr_matrix((np.asarray(data,dtype='float32'),np.asarray(indices),np.asarray(indptr)),shape=(len(genes),len(cells))).T
   objects.append(collapse(X,genes,cells,sample));print('LOADED',sample,len(cells),flush=True)
 return objects

def load183():
 objects=[];raw=ROOT/'00_raw/GSE183682';files=sorted(raw.glob('*matrix.mtx.gz'));assert len(files)==4
 tcrmap={'GSM5567501':'GSM5567509','GSM5567502':'GSM5567510','GSM5567503':'GSM5567511','GSM5567504':'GSM5567512'}
 for f in files:
  prefix=f.name[:-len('matrix.mtx.gz')];sample=f.name.split('_')[0]
  genes=pd.read_csv(raw/(prefix+'features.tsv.gz'),header=None,sep='\t')[1].astype(str)
  cells=pd.read_csv(raw/(prefix+'barcodes.tsv.gz'),header=None,sep='\t')[0].astype(str).values
  with gzip.open(f,'rb') as h:X=sio.mmread(h).T.tocsr()
  a=collapse(X,genes,cells,sample);present=a.var_names.isin(BROAD);keep=np.asarray((a.X[:,present]>0).sum(1)).ravel()>0
  tc=list(raw.glob(tcrmap[sample]+'*TCR*filtered_contig_annotations.csv.gz'))
  if tc:
   v=pd.read_csv(tc[0]);good=v.chain.isin(['TRA','TRB'])
   for k in ['productive','high_confidence']:good &= v[k].astype(str).str.lower().eq('true')
   a.obs['TCRpos']=a.obs.barcode.isin(v.loc[good,'barcode']).values;keep|=a.obs.TCRpos.values
  objects.append(a[keep].copy());print('LOADED',sample,int(keep.sum()),flush=True)
 return objects

def analyze(g,objects):
 d=OUT/g;d.mkdir();fig=d/'figures';fig.mkdir()
 features={str(a.obs['sample'].iloc[0]):set(a.var_names) for a in objects};common=set.intersection(*features.values())
 (d/'feature_audit.json').write_text(json.dumps({'features_per_sample':{k:len(v) for k,v in features.items()},'shared_features':len(common),'strategy':'intersection; no zero filling of missing genes'},indent=2))
 a=ad.concat(objects,join='inner',merge='same');a.X=sparse.csr_matrix(a.X)
 # QC uses only measured features shared across all libraries in this cohort.
 total=np.asarray(a.X.sum(1)).ravel();nf=np.asarray((a.X>0).sum(1)).ravel()
 mt=a.var_names.str.startswith(('MT-','MT.'));mt_count=np.asarray(a.X[:,mt].sum(1)).ravel()
 a.obs['nCount']=total;a.obs['nGene']=nf;a.obs['pctMT']=100*mt_count/np.maximum(total,1)
 a.obs['QCpass']=(nf>=200)&(total>=300)&(a.obs.pctMT<20)
 a.obs.to_csv(d/'candidate_QC.tsv.gz',sep='\t')
 a=a[a.obs.QCpass].copy();assert a.n_obs>=100
 a.layers['counts']=a.X.copy();sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a)
 sc.pp.highly_variable_genes(a,n_top_genes=min(3000,a.n_vars-1),flavor='seurat')
 a.var.loc[a.var_names.str.startswith(('MT-','MT.','RPL','RPS'))|a.var_names.isin(['MRE11','MRE11A']),'highly_variable']=False
 n_hvg=int(a.var.highly_variable.sum());assert n_hvg>=50
 b=a[:,a.var.highly_variable].copy();sc.pp.scale(b,max_value=10)
 sc.tl.pca(b,n_comps=min(30,b.n_obs-1,b.n_vars-1),random_state=42)
 a.obsm['X_pca']=b.obsm['X_pca'];sc.pp.neighbors(a,n_neighbors=min(20,a.n_obs-1),use_rep='X_pca',random_state=42)
 sc.tl.umap(a,random_state=42)
 for res,key in [(0.4,'R04'),(0.8,'R08'),(1.2,'R12')]:sc.tl.leiden(a,resolution=res,key_added=key,random_state=42)
 panel=[k for k in PANEL if k in a.var_names];C=a.layers['counts'];X=a.X
 values={k:np.asarray(C[:,a.var_names==k].sum(1)).ravel() for k in panel}
 for k in ['MRE11','MRE11A']:
  if k in values:a.obs[k]=values[k]
 a.obs.to_csv(d/'cell_metadata.tsv.gz',sep='\t')
 pd.DataFrame(a.obsm['X_umap'],index=a.obs_names,columns=['UMAP1','UMAP2']).to_csv(d/'umap.tsv.gz',sep='\t')
 rows=[];tops=[]
 for key in ['R04','R08','R12']:
  pd.crosstab(a.obs[key],a.obs['sample']).to_csv(d/(key+'_sample_counts.tsv'),sep='\t')
  for cluster in a.obs[key].cat.categories:
   mask=np.asarray(a.obs[key]==cluster);other=~mask
   for gene in panel:
    rows.append(dict(resolution=key,cluster=cluster,gene=gene,n_cells=int(mask.sum()),pct_detected=100*np.mean(values[gene][mask]>0),mean_count=np.mean(values[gene][mask])))
   if other.sum():
    mean=np.asarray(X[mask].mean(0)).ravel();rest=np.asarray(X[other].mean(0)).ravel();pct=np.asarray((C[mask]>0).mean(0)).ravel()
    score=mean-rest;order=np.argsort(-np.where(pct>=0.1,score,-np.inf))[:30]
    for j in order:
     if pct[j]>=0.1:tops.append(dict(resolution=key,cluster=cluster,gene=a.var_names[j],mean_log=mean[j],rest_mean_log=rest[j],difference=score[j],pct_detected=100*pct[j]))
 pd.DataFrame(rows).to_csv(d/'cluster_marker_detection.tsv',sep='\t',index=False)
 pd.DataFrame(tops).to_csv(d/'cluster_top30.tsv',sep='\t',index=False)
 xy=a.obsm['X_umap']
 for key in ['R04','R08','R12','sample']:
  f,ax=plt.subplots(figsize=(7,5));cats=pd.Categorical(a.obs[key]);colors=plt.get_cmap('tab20')
  for i,cat in enumerate(cats.categories):
   z=np.asarray(cats==cat);label=str(cat).replace('GSM53248','H').replace('GSM55675','M') if key=='sample' else str(cat)
   ax.scatter(xy[z,0],xy[z,1],s=3,c=[colors(i%20)],label=label,rasterized=True)
  ax.set(xlabel='UMAP1',ylabel='UMAP2',title=key);ax.legend(markerscale=3,fontsize=8,bbox_to_anchor=(1.01,1));f.tight_layout();f.savefig(fig/(key+'_UMAP.png'),dpi=300);plt.close(f)
 for gene in panel:
  f,ax=plt.subplots(figsize=(5,4));v=np.asarray(X[:,a.var_names==gene].sum(1)).ravel();order=np.argsort(v)
  p=ax.scatter(xy[order,0],xy[order,1],c=v[order],s=3,cmap='viridis',rasterized=True);f.colorbar(p,ax=ax,label='LogNorm');ax.set(xlabel='UMAP1',ylabel='UMAP2',title=gene);f.tight_layout();f.savefig(fig/(gene+'_UMAP.png'),dpi=250);plt.close(f)
 labels={s:s.replace('GSM53248','H').replace('GSM55675','M') for s in a.obs['sample'].unique()}
 (d/'figure_captions.txt').write_text('UMAP1/2: first/second UMAP coordinates. R04/R08/R12: Leiden resolution 0.4/0.8/1.2. Each point is a QC-passing broad T/NK candidate, not an assigned CD4 cell. Gene panels: log(1+counts/shared-feature library total*10000). LogNorm: log-normalized expression. Sample labels: '+json.dumps(labels)+'. Genes are uppercase symbols. Sample and cluster colors are categorical. Embedding excludes MRE11/MRE11A, mitochondrial and ribosomal HVGs. No batch correction or final doublet exclusion applied.\n')
 a.write_h5ad(d/'TNK_reclustered.h5ad',compression='gzip')
 print('DATASET_PASS',g,a.n_obs,flush=True)
for g,loader in [('GSE174725',load174),('GSE183682',load183)]:
 try:analyze(g,loader());status.append({'dataset':g,'status':'PASS'})
 except Exception as e:traceback.print_exc();status.append({'dataset':g,'status':'FAIL','error':str(e)})
(OUT/'status.json').write_text(json.dumps(status,indent=2))
with tarfile.open(ROOT/'tmp/screen2_TNK_return.tar.gz','w:gz') as t:
 for p in OUT.rglob('*'):
  if p.is_file() and p.suffix!='.h5ad':t.add(p,arcname=str(p.relative_to(OUT)))
passed=all(x['status']=='PASS' for x in status);print('SCREEN2_TNK_PASS' if passed else 'SCREEN2_TNK_PARTIAL',flush=True);sys.exit(0 if passed else 1)
