from pathlib import Path
import pandas as pd,numpy as np,scipy.sparse as sp,anndata as ad,re
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/three_cohort_four_axis_20261017_v16');P=O.parent.parent
f=pd.read_csv(O/'tables/function_membership_extended.tsv',sep='\t');d=pd.read_csv(O/'tables/DE_all_models.tsv',sep='\t',low_memory=False);genes=set(d.loc[d.model.isin(['314_disease','joint_stratified','274_time'])&(d.PValue<.05)&d.gene.isin(f.gene),'gene'])|{'IKZF2','PENK','IL1R2','TNFRSF9','GSN','DUSP5','STAT1','GSDMD'}
rows=[]
for ds,sc in [('GSE314943','Core'),('GSE183682','Core'),('GSE274372','TLike')]:
 stem=ds+'__'+sc;m=pd.read_csv(O/'inputs'/f'{stem}.cells.tsv.gz',sep='\t');g=pd.read_csv(O/'inputs'/f'{stem}.genes.txt',header=None)[0].tolist();ix={v:i for i,v in enumerate(g)};x=sp.load_npz(O/'inputs'/f'{stem}.cellcounts.npz');tidx=[ix[v] for v in ['CD3D','CD3E','CD3G','TRAC','CD247'] if v in ix];tpos=np.asarray(x[:,tidx].sum(1)).ravel()>0
 for scope in ([sc,'Naive','ActCD4','Treg'] if sc=='Core' else [sc]):
  take=np.ones(len(m),bool) if scope==sc else m.subtype.eq(scope).to_numpy()
  for sample,mm in m.loc[take].groupby('sample'):
   ids=mm.index.to_numpy();base=dict(dataset=ds,scope=scope,sample=sample,group=mm.group.iloc[0],cells=len(ids))
   for gene in sorted(genes):
    if gene not in ix:rows.append({**base,'gene':gene,'feature_present':False});continue
    vals=x[ids,ix[gene]].toarray().ravel();pos=vals>0;multi=vals>=2;bg=mm.background_fraction.to_numpy();library=mm.library.to_numpy();clean=bg<=.01
    rows.append({**base,'gene':gene,'feature_present':True,'UMI':vals.sum(),'CP10k':vals.sum()/library.sum()*1e4,'positive_cells':int(pos.sum()),'fraction_positive':pos.mean(),'multiUMI_cells':int(multi.sum()),'positive_oneUMI_fraction':float((vals[pos]==1).mean()) if pos.any() else np.nan,'positive_Tmarker_fraction':float(tpos[ids][pos].mean()) if pos.any() else np.nan,'positive_bg_median':float(np.median(bg[pos])) if pos.any() else np.nan,'BG1_cells':int(clean.sum()),'BG1_UMI':vals[clean].sum(),'BG1_CP10k':vals[clean].sum()/library[clean].sum()*1e4 if clean.any() else np.nan,'BG1_positive_fraction':pos[clean].mean() if clean.any() else np.nan})
pd.DataFrame(rows).to_csv(O/'tables/candidate_sample_expression_QC.tsv',sep='\t',index=False)
lr=pd.read_csv(O/'tables/expressed_LR_with_decoy_literature.tsv',sep='\t');ligands=set(';'.join(lr.ligand_parts).split(';'));a=ad.read_h5ad(P/'02_analysis/GSE314943/stage3_20260905_214228_VPtd0r/objects/whole_annotated.h5ad');raw=sp.csr_matrix(a.layers['counts']);sy=a.var.symbol.astype(str).str.upper().to_numpy();lib=np.asarray(raw.sum(1)).ravel();idx={v:np.flatnonzero(sy==v) for v in sorted(ligands|genes)};records=[]
for (sample,broad),mm in a.obs.groupby(['sample','Broad'],observed=True):
 ids=a.obs.index.get_indexer(mm.index);z=raw[ids];library=lib[ids].sum()
 for gene,inds in idx.items():
  vals=np.asarray(z[:,inds].sum(1)).ravel() if len(inds) else np.zeros(len(ids));records.append(dict(sample=sample,group=str(sample).split('_')[0],Broad=broad,cells=len(ids),gene=gene,feature_present=bool(len(inds)),UMI=vals.sum(),CP10k=vals.sum()/library*1e4,fraction_positive=(vals>0).mean()))
pd.DataFrame(records).to_csv(O/'tables/whole_lung_source_expression.tsv.gz',sep='\t',index=False)
print('QC_ROWS',len(rows),'SOURCE_ROWS',len(records));(O/'logs/QC_export_complete.txt').write_text('success',encoding='utf-8')