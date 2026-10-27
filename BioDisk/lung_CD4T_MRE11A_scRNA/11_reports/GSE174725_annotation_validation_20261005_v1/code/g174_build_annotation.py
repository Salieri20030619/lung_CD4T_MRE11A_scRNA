from pathlib import Path
import os,json,gzip,csv,collections,time
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1';os.environ['NUMBA_NUM_THREADS']='2'
import numpy as np,pandas as pd
from scipy import sparse
import anndata as ad,scanpy as sc
from threadpoolctl import threadpool_limits
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/GSE174725_annotation_validation_20261005_v1';R=O/'results';A=P/'11_reports/GSE174725_reaudit_20261005_v1/results';E=P/'01_raw_data/GSE174725/extracted'
def wr(d,n):pd.DataFrame(d).to_csv(R/n,sep='\t',index=False,na_rep='NA')
L='IL15 TNF CLSTN1 ADAM15 PLAU DSC2 CXCL2 ICAM2 VTN VEGFB GSTP1 RTN4'.split();TF=['IKZF2','STAT1']
PANELS={'T':'CD3D CD3E CD3G TRAC TRBC1 TRBC2 CD247','CD4state':'CD4 IL7R LTB TCF7 LEF1 CCR7 SELL AQP3 MAL','Treg':'FOXP3 IL2RA CTLA4 TIGIT IKZF4 LRRC32','CD8':'CD8A CD8B CCL5 GZMK GZMH PRF1','NK':'NKG7 GNLY KLRD1 NCR1 KLRF1 FCGR3A','Myeloid':'LYZ LST1 TYROBP FCER1G AIF1 CTSS SPI1','Mac':'C1QA C1QB C1QC CD68 CD163 CSF1R MRC1 MSR1 APOE','AM':'FABP4 PPARG MARCO INHBA NAPSA LPL ABCG1','Mono':'FCN1 S100A8 S100A9 S100A12 VCAN CTSD IL1B','DC':'CD1C FCER1A CLEC10A CLEC9A XCR1 LILRA4 GZMB','B':'MS4A1 CD79A CD79B CD74 CD37 CD19','Plasma':'MZB1 JCHAIN SDC1 IGKC','Epithelial':'EPCAM KRT8 KRT18 KRT19 SFTPC SFTPA1 SCGB1A1','Stromal':'COL1A1 COL1A2 DCN COL3A1 PECAM1 VWF','Platelet':'PPBP PF4 GP9','Proliferation':'MKI67 TOP2A STMN1'}
PANELS={k:v.split() for k,v in PANELS.items()};(O/'provenance/lineage_panels.json').write_text(json.dumps(PANELS,indent=2))
mp=pd.read_csv(A/'unambiguous_gene_mapping.tsv.gz',sep='\t');common=sorted(set.intersection(*[set(q.ensembl) for _,q in mp.groupby('sample')]));info=mp.drop_duplicates('ensembl').set_index('ensembl').loc[common];assert info.pooled_symbol.is_unique
symbols=info.pooled_symbol.astype(str).values;objects=[];featureqc=[]
legacy=pd.read_csv(P/'02_analysis/screen2_CD4_target_wz1lfx/GSE174725/scope_membership.tsv.gz',sep='\t').set_index(['sample','barcode'])
for path in sorted(E.glob('*.csv.gz')):
 gsm=path.name.split('_')[0];print('BUILD_CSV',gsm,flush=True)
 with gzip.open(path,'rt') as f:header=next(csv.reader(f))
 d=pd.read_csv(path,index_col=0,dtype={c:'float32' for c in header[1:]});blank=d.index.isna()&d.isna().all(axis=1);assert int(blank.sum())==(1 if gsm=='GSM5324868' else 0);d=d.loc[~blank]
 if gsm=='GSM5324866':assert np.isnan(d.loc['IFNGR2','GGCTTTCAGCCAGTAG-1']);d=d.drop(columns=['GGCTTTCAGCCAGTAG-1'])
 q=mp[mp['sample']==gsm].set_index('ensembl').loc[common];v=d.iloc[q.csv_row.to_numpy()].to_numpy().T;assert np.isfinite(v).all();X=sparse.csr_matrix(v)
 md=pd.read_csv(A/(gsm+'_cell_feasibility_SERVER_ONLY.tsv.gz'),sep='\t').set_index('barcode').loc[d.columns];md.index=[gsm+':'+b for b in d.columns];md['barcode']=d.columns;md['group']='Exposed' if gsm in ['GSM5324866','GSM5324867'] else 'Silicosis';md['old_CD4']='NotInCore';md['old_F08']='NotInTrefine'
 for b in d.columns:
  if (gsm,b) in legacy.index:
   row=legacy.loc[(gsm,b)];md.loc[gsm+':'+b,'old_F08']=str(row.F08)
   if row.scope_core:md.loc[gsm+':'+b,'old_CD4']=row.core_label
 a=ad.AnnData(X,obs=md,var=pd.DataFrame({'ensembl':common},index=symbols));objects.append(a[a.obs.QC].copy());featureqc.append(dict(sample=gsm,CSVcells=d.shape[1],common_features=len(common),QCcells=int(md.QC.sum()),old_Core_QC=int(((md.old_CD4!='NotInCore')&md.QC).sum())));del d,v,X,a
wr(featureqc,'input_common_feature_QC.tsv');a=ad.concat(objects,join='inner',merge='same');del objects;a.layers['counts']=a.X.copy();a.obs['common_counts']=np.asarray(a.X.sum(1)).ravel();a.obs['common_n_features']=np.asarray((a.X>0).sum(1)).ravel();assert a.obs_names.is_unique
a.write_h5ad(R/'common_counts_QC_SERVER_ONLY.h5ad',compression='gzip');print('COMMON_COUNTS_SAVED',a.shape,flush=True)
# Scanpy includes an in-tree implementation in installed1.11; use it only if its dependencies work.
doublets=[];a.obs['doublet_score']=np.nan;a.obs['predicted_doublet']=False;a.obs['doublet_status']='UNCOMPUTED'
for gsm in a.obs['sample'].unique():
 idx=a.obs['sample']==gsm;b=a[idx].copy();print('SCRUBLET_START',gsm,flush=True)
 try:
  with threadpool_limits(limits=2):sc.pp.scrublet(b,expected_doublet_rate=.06,random_state=42,verbose=False)
  a.obs.loc[idx,'doublet_score']=b.obs.doublet_score.values;a.obs.loc[idx,'predicted_doublet']=b.obs.predicted_doublet.values;a.obs.loc[idx,'doublet_status']='COMPUTED'
  doublets.append(dict(sample=gsm,n=b.n_obs,flagged=int(b.obs.predicted_doublet.sum()),threshold=float(b.uns['scrublet']['threshold']),status='COMPUTED'))
 except Exception as e:doublets.append(dict(sample=gsm,n=b.n_obs,flagged=None,status='UNAVAILABLE',error=str(e)))
 del b
wr(doublets,'scrublet_per_donor_audit.tsv');print('SCRUBLET_DONE',json.dumps(doublets),flush=True)
a.obs.to_csv(R/'pregraph_doublet_context_SERVER_ONLY.tsv.gz',sep='\t')
sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a);sc.pp.highly_variable_genes(a,n_top_genes=2500,flavor='seurat');a.var.loc[a.var_names.str.startswith(('MT-','RPL','RPS'))|a.var_names.isin(L+TF),'highly_variable']=False
b=a[:,a.var.highly_variable].copy()
with threadpool_limits(limits=2):
 sc.pp.scale(b,max_value=10);sc.tl.pca(b,n_comps=30,random_state=42);a.obsm['X_pca']=b.obsm['X_pca'];del b;sc.pp.neighbors(a,n_neighbors=20,use_rep='X_pca',random_state=42);sc.tl.leiden(a,resolution=.6,key_added='broad06',random_state=42,flavor='igraph',n_iterations=2);sc.tl.umap(a,random_state=42)
print('GRAPH_DONE',a.obs.broad06.value_counts().to_dict(),flush=True)
a.write_h5ad(R/'graph_checkpoint_SERVER_ONLY.h5ad',compression='gzip')
C=a.layers['counts'];X=a.X;allmarkers=sorted(set(sum(PANELS.values(),[])));panelcoverage=[]
for name,genes in PANELS.items():
 measured=[g for g in genes if g in a.var_names];pos=a.var_names.isin(measured);a.obs['n_'+name]=np.asarray((C[:,pos]>0).sum(1)).ravel();a.obs['UMI_'+name]=np.asarray(C[:,pos].sum(1)).ravel();a.obs['score_'+name]=np.asarray(X[:,pos].mean(1)).ravel()
 for g in genes:panelcoverage.append(dict(panel=name,gene=g,in_common=g in a.var_names))
wr(panelcoverage,'lineage_panel_coverage.tsv')
a.obs['strong_T_myeloid_mixed']=(a.obs.n_T>=3)&(a.obs.n_Myeloid>=3)&(a.obs.n_Mac>=2)&(a.obs.UMI_Mac>=10)
a.obs['strong_T_B_mixed']=(a.obs.n_T>=3)&(a.obs.n_B>=3)&(a.obs.UMI_B>=10)
a.obs['strong_myeloid_epithelial_mixed']=(a.obs.n_Myeloid>=3)&(a.obs.n_Epithelial>=3)&(a.obs.UMI_Epithelial>=10)
a.obs['lineage_mixed_flag']=a.obs[['strong_T_myeloid_mixed','strong_T_B_mixed','strong_myeloid_epithelial_mixed']].any(axis=1)
rows=[];tops=[];summ=[];sample=[]
for cl,obs in a.obs.groupby('broad06',observed=True):
 z=np.asarray(a.obs.broad06==cl);avg=np.asarray(X[z].mean(0)).ravel();rest=np.asarray(X[~z].mean(0)).ravel();pct=np.asarray((C[z]>0).mean(0)).ravel();diff=avg-rest
 for j in np.argsort(-np.where(pct>=.1,diff,-np.inf))[:40]:tops.append(dict(cluster=cl,gene=a.var_names[j],mean_log_diff=diff[j],pct=pct[j]))
 for g in allmarkers:
  if g in a.var_names:
   j=a.var_names.get_loc(g);rows.append(dict(cluster=cl,gene=g,cells=int(z.sum()),pct=pct[j],mean_log=avg[j],mean_count=float(C[z,j].mean())))
 donors=obs['sample'].value_counts();v=dict(cluster=cl,cells=len(obs),donors=donors.size,dominant_donor=donors.index[0],dominant_fraction=donors.iloc[0]/len(obs),median_counts=obs.total_counts.median(),median_features=obs.n_features.median(),doublet_fraction=obs.predicted_doublet.mean(),mixed_fraction=obs.lineage_mixed_flag.mean())
 for name in PANELS:v['mean_'+name]=obs['score_'+name].mean();v['n2_'+name]=(obs['n_'+name]>=2).mean()
 summ.append(v)
 for gsm,n in donors.items():sample.append(dict(cluster=cl,sample=gsm,cells=n))
wr(rows,'broad_cluster_marker_profiles.tsv');wr(tops,'broad_cluster_top40.tsv');wr(summ,'broad_cluster_summary.tsv');wr(sample,'broad_cluster_donor_counts.tsv');a.obs.to_csv(R/'broad_cell_context_SERVER_ONLY.tsv.gz',sep='\t');a.write_h5ad(R/'broad_annotation_review_SERVER_ONLY.h5ad',compression='gzip')
# Preserve reused T-cluster evidence at donor resolution in canonical full common space.
tr=[]
for (label,gsm),obs in a.obs[a.obs.old_CD4!='NotInCore'].groupby(['old_CD4','sample']):
 z=a.obs_names.isin(obs.index)
 for g in allmarkers:
  if g in a.var_names:
   j=a.var_names.get_loc(g);tr.append(dict(label=label,sample=gsm,cells=len(obs),gene=g,pct=float((C[z,j]>0).mean()),mean_count=float(C[z,j].mean())))
wr(tr,'existing_Core_donor_marker_review.tsv');(O/'provenance/build_complete.json').write_text(json.dumps({'status':'COMPLETE','cells':a.n_obs,'genes':a.n_vars,'clusters':int(a.obs.broad06.nunique()),'new_graphs':1,'reused_T_labels':True,'doublet_audit':doublets},indent=2));print('ANNOTATION_REVIEW_READY',flush=True)
