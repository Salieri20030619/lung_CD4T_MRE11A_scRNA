from pathlib import Path
import os,json
os.environ['OPENBLAS_NUM_THREADS']='1'
import pandas as pd,numpy as np,anndata as ad
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/GSE174725_annotation_validation_20261005_v1';R=O/'results';a=ad.read_h5ad(R/'broad_annotation_review_SERVER_ONLY.h5ad');C=a.layers['counts']
extra='CD14 FCGR3B CSF3R CXCR2 FPR1 FPR2 CEACAM8 RETN LTF CLC PRG2 PRG3 RNASE2 RNASE3 SIGLEC8 TPSAB1 TPSB2 KIT CPA3 CLEC4C CD1C FCER1A CLEC10A CLEC9A XCR1 LILRA4 TCF4 FABP4 PPARG MARCO INHBA NAPSA LPL ABCG1 C1QA C1QB C1QC CD68 CD163 CSF1R MRC1 MSR1 APOE CD3D CD3E TRAC CD4 CD8A CD8B EPCAM KRT8 KRT18 SCGB1A1 PPBP PF4 GP9'.split();rows=[]
for cl in a.obs.broad06.cat.categories:
 z=np.asarray(a.obs.broad06==cl)
 for g in extra:
  if g in a.var_names:
   j=a.var_names.get_loc(g);x=C[z,j].toarray().ravel();rows.append(dict(cluster=cl,gene=g,cells=int(z.sum()),pct=float((x>0).mean()),mean_count=float(x.mean())))
d=pd.DataFrame(rows);d.to_csv(R/'additional_lineage_marker_review.tsv',sep='\t',index=False)
clusters=['0','1','6','9','17','18','19'];q=d[d.cluster.isin(clusters)];print('MYELOID_PCT',q.pivot(index='gene',columns='cluster',values='pct').round(3).to_string());print('MYELOID_MEAN',q.pivot(index='gene',columns='cluster',values='mean_count').round(2).to_string())
ct=pd.crosstab(a.obs.broad06,a.obs['sample']);ct.to_csv(R/'broad_donor_table.tsv',sep='\t');print('COUNTS',ct.to_string());print('CORE_CLUSTER',pd.crosstab(a.obs.broad06,a.obs.old_CD4).to_string())
t=pd.read_csv(R/'broad_cluster_top40.tsv',sep='\t');print('TOP40',t[t.cluster.astype(str).isin(clusters)].groupby('cluster').gene.apply(lambda x:','.join(x)).to_string())
# Correct categorical value_counts includes zero levels in initial display; no labels or cells change.
s=pd.read_csv(R/'broad_cluster_summary.tsv',sep='\t');s['donors']=s.cluster.astype(str).map((ct>0).sum(axis=1).rename(index=str));s.to_csv(R/'broad_cluster_summary.tsv',sep='\t',index=False)
a.obs['T_platelet_mixed_flag']=(a.obs.n_T>=3)&(a.obs.n_Platelet>=2)&(a.obs.UMI_Platelet>=10);a.obs.to_csv(R/'annotation_review_cell_context_SERVER_ONLY.tsv.gz',sep='\t');print('CORE_FLAGS',a.obs[a.obs.old_CD4!='NotInCore'].groupby(['old_CD4','sample'],observed=True).agg(cells=('sample','size'),doublet=('predicted_doublet','sum'),mixed=('lineage_mixed_flag','sum'),platelet=('T_platelet_mixed_flag','sum')).to_string());print('EXTRA_REVIEW_COMPLETE')
