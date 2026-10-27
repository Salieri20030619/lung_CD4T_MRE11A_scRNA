import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2',NUMBA_NUM_THREADS='2')
from pathlib import Path
import numpy as np,pandas as pd,scanpy as sc,anndata as ad,json,gzip
from scipy import sparse,io
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/death_priority_20261018_v2';raw=P/'00_raw/GSE183682'
objs=[];qcs=[]
mp={'GSM5567501':('Saline',7),'GSM5567502':('Silica',7),'GSM5567503':('Saline',56),'GSM5567504':('Silica',56)}
for f in sorted(raw.glob('*matrix.mtx.gz')):
 prefix=f.name[:-len('matrix.mtx.gz')];sam=f.name.split('_')[0]
 gs=pd.read_csv(raw/(prefix+'features.tsv.gz'),sep='\t',header=None)[1].astype(str).str.upper().to_numpy();bc=pd.read_csv(raw/(prefix+'barcodes.tsv.gz'),sep='\t',header=None)[0].astype(str).to_numpy()
 with gzip.open(f,'rb') as h:x=io.mmread(h).T.tocsr()
 genes,ii=np.unique(gs,return_inverse=True);proj=sparse.csr_matrix((np.ones(len(gs)),(np.arange(len(gs)),ii)),shape=(len(gs),len(genes)));x=(x@proj).tocsr().astype('float32')
 a=ad.AnnData(x,obs=pd.DataFrame(index=[sam+':'+b for b in bc]),var=pd.DataFrame(index=genes))
 a.obs['sample']=sam;a.obs['group']=mp[sam][0];a.obs['day']=mp[sam][1];a.obs['library']=np.asarray(x.sum(1)).ravel();a.obs['n_genes']=np.asarray((x>0).sum(1)).ravel();a.obs['pct_mt']=100*np.asarray(x[:,np.char.startswith(genes.astype(str),'MT-')].sum(1)).ravel()/np.maximum(a.obs.library,1)
 keep=(a.obs.n_genes>=200)&(a.obs.library>=300)&(a.obs.pct_mt<20);nraw=len(a);a=a[keep].copy()
 print('QC',sam,nraw,len(a),flush=True)
 sc.pp.scrublet(a,expected_doublet_rate=.06,random_state=20261018,n_prin_comps=30,verbose=False)
 a.obs.to_csv(O/'tables'/(sam+'_whole_QC.tsv.gz'),sep='\t');nd=int(a.obs.predicted_doublet.sum());a=a[~a.obs.predicted_doublet].copy()
 qcs.append(dict(sample=sam,raw=nraw,qc=int(keep.sum()),doublets=nd,singlets=len(a),threshold=a.uns['scrublet'].get('threshold')))
 objs.append(a)
a=ad.concat(objs,join='inner');a.layers['counts']=a.X.copy();sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a);sc.pp.highly_variable_genes(a,n_top_genes=2500,batch_key='sample');sc.pp.pca(a,n_comps=30,random_state=20261018)
import harmonypy as hm
h=hm.run_harmony(a.obsm['X_pca'],a.obs,'sample',random_state=20261018,max_iter_harmony=30,verbose=False);z=h.Z_corr;a.obsm['X_harmony']=z.T if z.shape[0]!=a.n_obs else z
sc.pp.neighbors(a,n_neighbors=15,use_rep='X_harmony');sc.tl.leiden(a,resolution=.6,key_added='cluster',random_state=20261018);sc.tl.umap(a,random_state=20261018)
panels={'Mac':['C1QA','C1QB','C1QC','ADGRE1','MERTK','MARCO','MSR1'],'Myeloid':['LYZ2','TYROBP','FCER1G','CSF1R','CD68'],'Mono':['LY6C2','CCR2','S100A8','S100A9','VCAN'],'DC':['FLT3','ZBTB46','XCR1','CLEC9A','CD209A'],'T':['CD3D','CD3E','TRAC'],'B':['MS4A1','CD79A','CD79B'],'NK':['NKG7','NCR1','KLRD1'],'Neut':['CSF3R','S100A8','S100A9','RETNLG'],'Epi':['EPCAM','SFTPC','KRT19'],'Endo':['PECAM1','VWF','KDR'],'Fib':['COL1A1','COL1A2','DCN'],'AM':['PPARG','FABP4','CAR4','SIGLECF'],'Mast':['KIT','CPA3','TPSB2']}
markers=sorted(set(sum(panels.values(),[])));out=[]
for (cl,sam),mm in a.obs.groupby(['cluster','sample'],observed=True):
 ix=a.obs_names.get_indexer(mm.index);row=dict(cluster=str(cl),sample=sam,cells=len(mm))
 for g in markers:
  if g in a.var_names:
   gi=a.var_names.get_loc(g);v=a.X[ix,gi].toarray().ravel();row[g]=v.mean();row[g+'_pct']=float((v>0).mean())
 for name,gs in panels.items():row['score_'+name]=float(a[ix,[g for g in gs if g in a.var_names]].X.mean())
 out.append(row)
pd.DataFrame(out).to_csv(O/'tables/GSE183682_cluster_markers.tsv',sep='\t',index=False);pd.DataFrame(qcs).to_csv(O/'tables/GSE183682_QC_summary.tsv',sep='\t',index=False)
sc.tl.rank_genes_groups(a,'cluster',method='wilcoxon',n_genes=30,use_raw=False);sc.get.rank_genes_groups_df(a,group=None).to_csv(O/'tables/GSE183682_cluster_topmarkers_descriptive.tsv',sep='\t',index=False)
a.write_h5ad(O/'objects/GSE183682_whole_clustered.h5ad')
(O/'logs/whole_clustering_complete.json').write_text(json.dumps({'status':'success','cells':len(a),'clusters':int(a.obs.cluster.nunique())}),encoding='utf-8')
print('CLUSTERED',len(a),flush=True)
