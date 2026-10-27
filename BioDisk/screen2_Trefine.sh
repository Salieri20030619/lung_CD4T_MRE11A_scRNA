#!/usr/bin/env bash
set -Eeuo pipefail
export PROJECT_ROOT=/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA
mkdir -p "$PROJECT_ROOT/tmp" "$PROJECT_ROOT/02_analysis" "$PROJECT_ROOT/logs"
exec 9>"$PROJECT_ROOT/tmp/screen2_Trefine.lock"
if ! flock -n 9; then echo SCREEN2_TREFINE_ALREADY_RUNNING; exit 0; fi
trap 'rc=$?; echo "SCREEN2_TREFINE_FAIL EXIT_CODE=$rc"; exit "$rc"' ERR
TASK=$(mktemp -d "$PROJECT_ROOT/02_analysis/screen2_Trefine_XXXXXX")
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMBA_NUM_THREADS=4
cat > "$TASK/refine.py" <<'REFINE_PY'
from pathlib import Path
import os,sys,json,tarfile,traceback
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import sparse
ROOT=Path(os.environ['PROJECT_ROOT']);OUT=Path(sys.argv[1]);SOURCE=ROOT/'02_analysis/screen2_TNK_20260911_103239_nfW0yu'
CONFIG={'GSE174725':['0','1','4','5','9','10','12','13','15','16','25'],'GSE183682':['2','4','5','9','14','18']}
GENES='CD3D CD3E CD3G TRAC CD4 CD8A CD8B CD8B1 IL7R TCF7 LEF1 CCR7 SELL FOXP3 IL2RA CTLA4 ICOS TNFRSF4 GATA3 TBX21 RORC TRDC TRGC1 TRGC2 TRDV4 ZBTB16 NKG7 GNLY NCR1 KLRD1 GZMB GZMK CCL5 MKI67 TOP2A STMN1 LYZ LYZ2 LST1 TYROBP FCER1G S100A8 S100A9 PPBP PF4 MS4A1 CD79A EPCAM SFTPC SCGB1A1 HBB HBA2 HBB-BS MRE11 MRE11A'.split()
status=[]
for g,clusters in CONFIG.items():
 try:
  d=OUT/g;d.mkdir(parents=True);a=sc.read_h5ad(SOURCE/g/'TNK_reclustered.h5ad')
  a.obs['parent']=a.obs.R08.astype(str);assert set(clusters)<=set(a.obs.parent)
  allsamples=sorted(a.obs['sample'].unique());a=a[a.obs.parent.isin(clusters)].copy()
  a.obs['prior_core']=a.obs.parent.isin(['1','4','5'] if g=='GSE174725' else ['5'])
  a.obs.to_csv(d/'input_cells.tsv.gz',sep='\t')
  # Rebuild embeddings from retained counts; MRE is excluded from HVGs.
  a.X=sparse.csr_matrix(a.layers['counts']).copy();a.obsm.clear();a.obsp.clear();a.uns.clear()
  a.var=pd.DataFrame(index=a.var_names.copy());sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a)
  sc.pp.highly_variable_genes(a,flavor='seurat',n_top_genes=min(2500,a.n_vars-1))
  a.var.loc[a.var_names.str.startswith(('MT-','MT.','RPL','RPS'))|a.var_names.isin(['MRE11','MRE11A']),'highly_variable']=False
  b=a[:,a.var.highly_variable].copy();assert b.n_vars>=50
  sc.pp.scale(b,max_value=10);sc.tl.pca(b,n_comps=min(30,b.n_vars-1,b.n_obs-1),random_state=42)
  a.obsm['X_pca']=b.obsm['X_pca'];sc.pp.neighbors(a,use_rep='X_pca',n_neighbors=20,random_state=42);sc.tl.umap(a,random_state=42)
  keys=['F04','F08','F12']
  for key,res in zip(keys,[.4,.8,1.2]):sc.tl.leiden(a,resolution=res,key_added=key,random_state=42)
  panel=[k for k in GENES if k in a.var_names];C=a.layers['counts'];X=a.X
  targets=[k for k in ['MRE11','MRE11A'] if k in a.var_names];assert len(targets)==1;target=targets[0]
  v=np.asarray(C[:,a.var_names==target].sum(1)).ravel();a.obs['target_count']=v
  values={k:np.asarray(C[:,a.var_names==k].sum(1)).ravel() for k in panel};rows=[];top=[];sample_rows=[]
  for key in keys:
   pd.crosstab(a.obs[key],a.obs['sample']).reindex(columns=allsamples,fill_value=0).to_csv(d/(key+'_sample_counts.tsv'),sep='\t')
   pd.crosstab(a.obs[key],a.obs.parent).to_csv(d/(key+'_parent_counts.tsv'),sep='\t')
   for cluster in a.obs[key].cat.categories:
    z=np.asarray(a.obs[key]==cluster)
    for gene in panel:rows.append(dict(resolution=key,cluster=cluster,gene=gene,n=int(z.sum()),pct=100*np.mean(values[gene][z]>0),mean_count=float(np.mean(values[gene][z]))))
    for sample in allsamples:
     q=z&np.asarray(a.obs['sample']==sample);n=int(q.sum());sample_rows.append(dict(resolution=key,cluster=cluster,sample=sample,n=n,target_positive=int((v[q]>0).sum()),target_pct=100*np.mean(v[q]>0) if n else None,target_sum=float(v[q].sum())))
    if (~z).sum():
     avg=np.asarray(X[z].mean(0)).ravel();rest=np.asarray(X[~z].mean(0)).ravel();pct=np.asarray((C[z]>0).mean(0)).ravel();dif=avg-rest
     for j in np.argsort(-np.where(pct>=.1,dif,-np.inf))[:40]:
      if pct[j]>=.1:top.append(dict(resolution=key,cluster=cluster,gene=a.var_names[j],log_mean_diff=float(dif[j]),pct=100*pct[j]))
  pd.DataFrame(rows).to_csv(d/'marker_detection.tsv',sep='\t',index=False)
  pd.DataFrame(top).to_csv(d/'top40.tsv',sep='\t',index=False)
  pd.DataFrame(sample_rows).to_csv(d/'target_by_cluster_sample.tsv',sep='\t',index=False)
  a.obs.to_csv(d/'cell_metadata.tsv.gz',sep='\t');xy=a.obsm['X_umap']
  pd.DataFrame(xy,index=a.obs_names,columns=['UMAP1','UMAP2']).to_csv(d/'umap.tsv.gz',sep='\t')
  for key in keys+['sample','parent']:
   f,ax=plt.subplots(figsize=(7,5));cat=pd.Categorical(a.obs[key]);colors=plt.get_cmap('turbo',len(cat.categories))
   for i,c in enumerate(cat.categories):
    z=np.asarray(cat==c);label=str(c).replace('GSM53248','H').replace('GSM55675','M');ax.scatter(xy[z,0],xy[z,1],s=3,color=colors(i),label=label,rasterized=True)
   ax.set(xlabel='UMAP1',ylabel='UMAP2',title=key);ax.legend(fontsize=8,markerscale=3,bbox_to_anchor=(1.01,1));f.tight_layout();f.savefig(d/(key+'_UMAP.png'),dpi=300);plt.close(f)
  # Multi-panel markers: common layout, independent color scales.
  for start in range(0,len(panel),12):
   genes=panel[start:start+12];f,axes=plt.subplots(3,4,figsize=(14,10))
   for ax in axes.ravel():ax.axis('off')
   for ax,gene in zip(axes.ravel(),genes):
    ax.axis('on');w=np.asarray(X[:,a.var_names==gene].sum(1)).ravel();order=np.argsort(w);p=ax.scatter(xy[order,0],xy[order,1],c=w[order],s=2,cmap='viridis');ax.set(title=gene,xlabel='UMAP1',ylabel='UMAP2');ax.set_xticks([]);ax.set_yticks([]);f.colorbar(p,ax=ax,shrink=.6)
   f.tight_layout();f.savefig(d/f'markers_{start//12+1}.png',dpi=250);plt.close(f)
  (d/'captions.txt').write_text('UMAP1/2 are UMAP coordinates. F04/F08/F12 are Leiden resolutions 0.4/0.8/1.2. Parent is the original R08 cluster. Gene color: log(1+normalized counts per 10000), independent scale per panel. Sample labels: '+json.dumps({s:s.replace('GSM53248','H').replace('GSM55675','M') for s in allsamples})+'. Selected regions contain clear T clusters and separately tracked ambiguous rescue clusters; final CD4 labels are not assigned. MRE11/MRE11A excluded from HVGs. No disease-group inference or p-values generated. No batch correction or final doublet removal performed.\n')
  a.write_h5ad(d/'T_refined.h5ad',compression='gzip');status.append({'dataset':g,'status':'PASS','cells':a.n_obs,'parent_clusters':clusters});print('REFINE_PASS',g,a.n_obs,flush=True)
 except Exception as e:traceback.print_exc();status.append({'dataset':g,'status':'FAIL','error':str(e)})
(OUT/'status.json').write_text(json.dumps(status,indent=2))
with tarfile.open(ROOT/'tmp/screen2_Trefine_return.tar.gz','w:gz') as t:
 for p in OUT.rglob('*'):
  if p.is_file() and p.suffix!='.h5ad':t.add(p,arcname=str(p.relative_to(OUT)))
ok=all(s['status']=='PASS' for s in status);print('SCREEN2_TREFINE_PASS' if ok else 'SCREEN2_TREFINE_PARTIAL',flush=True);sys.exit(0 if ok else 1)
REFINE_PY
conda run -n batch_bench python -m py_compile "$TASK/refine.py"
echo SYNTAX_CHECK=PASS
echo "OUTPUT=$TASK"
conda run --no-capture-output -n batch_bench python -u "$TASK/refine.py" "$TASK"
