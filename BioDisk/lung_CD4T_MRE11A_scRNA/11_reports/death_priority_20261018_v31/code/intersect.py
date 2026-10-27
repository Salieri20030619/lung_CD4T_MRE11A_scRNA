from pathlib import Path
import pandas as pd,numpy as np,json
from scipy import sparse,stats
import anndata as ad
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/death_priority_20261018_v31';B=P/'11_reports/three_cohort_four_axis_20261017_v16'
def w(a,n):a.to_csv(O/'tables'/n,sep='\t',index=False)
a=pd.read_csv(O/'tables/gene_trajectory_all.tsv.gz',sep='\t');primary=a[a.variant.eq('K15_R1')].copy();go=pd.read_csv(O/'tables/GO_priority_members.tsv',sep='\t');print('GO_MISSING_SYMBOL',int(go.gene.isna().sum()));go=go[go.gene.notna()].copy();go.gene=go.gene.str.upper();mods=pd.read_csv(O/'tables/common_module_members.tsv',sep='\t');fun={term:set(z.gene) for term,z in go.groupby('GOALL')};core=set(mods.loc[mods.module.isin(['ApoExec','PyroExec']),'gene']);t=pd.read_csv(O/'tables/time_all12_D28_vs_D7.tsv',sep='\t');tm=pd.read_csv(O/'tables/time_all12_time_omnibus.tsv',sep='\t');ts=pd.read_csv(O/'tables/time_min10_D28_vs_D7.tsv',sep='\t')
t=t[['gene','logFC','PValue','FDR']].rename(columns={'logFC':'time_logFC','PValue':'time_P','FDR':'time_FDR'});tm=tm[['gene','PValue','FDR']].rename(columns={'PValue':'time_omnibus_P','FDR':'time_omnibus_FDR'});ts=ts[['gene','logFC','PValue','FDR']].rename(columns={'logFC':'time_min10_logFC','PValue':'time_min10_P','FDR':'time_min10_FDR'})
r=primary.merge(t,on='gene',how='left').merge(tm,on='gene',how='left').merge(ts,on='gene',how='left');r['apoptosis_GO']=r.gene.isin(fun['GO:0006915']);r['pyroptosis_GO']=r.gene.isin(fun['GO:0070269']);r['intrinsic_apoptosis_GO']=r.gene.isin(fun['GO:0097193']);r['autophagy_GO']=r.gene.isin(fun['GO:0006914']);r['mitophagy_GO']=r.gene.isin(fun['GO:0000422']);r['core_execution_panel']=r.gene.isin(core);r['death_function']=r.apoptosis_GO|r.pyroptosis_GO;r['trajectory_nominal']=r.P<.05;r['time_nominal']=r.time_P<.05;r['primary_intersection']=r.trajectory_nominal&r.time_nominal;r['primary_death_candidate']=r.primary_intersection&r.death_function;r['secondary_death_candidate']=r.trajectory_nominal&(r.time_omnibus_P<.05)&r.death_function;r['direction_concordant']=np.sign(r.mean_slope)==np.sign(r.time_logFC)
variants=a.groupby(['scope','gene']).apply(lambda z:pd.Series({'variants_evaluable':len(z),'variants_P_lt05':int((z.P<.05).sum()),'variant_slope_sign_agreement':float((np.sign(z.mean_slope)==np.sign(z.loc[z.variant.eq('K15_R1'),'mean_slope'].iloc[0])).mean()),'variant_P_max':z.P.max()}),include_groups=False).reset_index();r=r.merge(variants,on=['scope','gene']);w(r,'all_gene_intersection_membership.tsv.gz')
main=r[r.primary_death_candidate].copy();secondary=r[r.secondary_death_candidate].copy();allinter=r[r.primary_intersection].copy();w(allinter,'primary_intersection_all_genes.tsv');w(secondary,'secondary_time_omnibus_candidates.tsv')
main['function_tier']=np.where(main.core_execution_panel,'核心执行成员',np.where(main.intrinsic_apoptosis_GO,'内源性凋亡注释','凋亡调控/相关注释'));main['priority_state']=main.scope.eq('Treg');main=main.sort_values(['core_execution_panel','intrinsic_apoptosis_GO','priority_state','variants_P_lt05','P'],ascending=[False,False,False,False,True]);w(main,'primary_death_candidates.tsv')
# Candidate source QC and pseudotime curves. No source concern silently removes a candidate.
genes=sorted(set(main.gene)|set(secondary.gene));qc=[];curve=[]
for ds,scop in [('GSE314943','Core'),('GSE183682','Core'),('GSE274372','TLike')]:
 stem=ds+'__'+scop;x=sparse.load_npz(B/'inputs'/(stem+'.cellcounts.npz'));m=pd.read_csv(B/'inputs'/(stem+'.cells.tsv.gz'),sep='\t');g=pd.read_csv(B/'inputs'/(stem+'.genes.txt'),header=None)[0].to_list();ix={k:i for i,k in enumerate(g)};tidx=[ix[k] for k in ['CD3D','CD3E','TRAC'] if k in ix];tpos=np.asarray((x[:,tidx]>0).sum(1)).ravel()>0
 for gene in genes:
  if gene not in ix:continue
  vv=x[:,ix[gene]].toarray().ravel()
  for (sample,state),z in m.groupby(['sample','subtype']):
   ids=z.index.to_numpy();pos=ids[vv[ids]>0];clean=ids[z.background_fraction.to_numpy()<=.01];qc.append({'dataset':ds,'sample':sample,'scope':state,'group':z.group.iloc[0],'day':z.day.iloc[0],'gene':gene,'cells':len(ids),'positive':len(pos),'positive_Tmarker':int(tpos[pos].sum()),'positive_BG1':int((m.background_fraction.iloc[pos]<=.01).sum()),'CP10k':vv[ids].sum()/z.library.sum()*1e4,'BG1_CP10k':vv[clean].sum()/m.library.iloc[clean].sum()*1e4 if len(clean) else np.nan,'BG1_cells':len(clean)})
 for track in ['Naive_ActCD4','Treg'] if ds!='GSE274372' else []:
  coord=pd.read_csv(O/'tables'/(track+'_cell_trajectory.tsv.gz'),sep='\t');coord=coord[coord.dataset.eq(ds)].copy();loc=pd.Index(m.cell.astype(str)).get_indexer(coord.cell.astype(str));assert (loc>=0).all()
  for gene in genes:
   if gene not in ix:continue
   yy=np.log1p(x[loc,ix[gene]].toarray().ravel()/m.library.iloc[loc].to_numpy()*1e4)
   for state,z0 in coord.groupby('subtype'):
    valid=main[(main.gene==gene)&(main.scope==state)]
    if not len(valid):continue
    # Scale within full joint state exactly as primary slope model.
    full=pd.read_csv(O/'tables'/(track+'_cell_trajectory.tsv.gz'),sep='\t');vals=full.loc[full.subtype.eq(state),'K15_R1'];low,high=vals.min(),vals.max();pt=(coord.K15_R1.to_numpy()-low)/(high-low)
    for sample,z in coord[coord.subtype.eq(state)].groupby('sample'):
     ii=coord.index.get_indexer(z.index);bn=np.minimum((pt[ii]*4).astype(int),3)
     for k in range(4):
      jj=ii[bn==k]
      if len(jj):curve.append({'gene':gene,'scope':state,'dataset':ds,'sample':sample,'group':z.group.iloc[0],'bin':k,'DPT':pt[jj].mean(),'expression':yy[jj].mean(),'cells':len(jj)})
w(pd.DataFrame(qc),'candidate_sample_expression_QC.tsv');w(pd.DataFrame(curve),'candidate_state_curves.tsv')
# Preserve primary rooting in geometry objects, without recomputing graph or statistics.
for track in ['Naive_ActCD4','Treg']:
 f=O/'objects'/(track+'_geometry.h5ad');z=ad.read_h5ad(f);c=pd.read_csv(O/'tables'/(track+'_cell_trajectory.tsv.gz'),sep='\t').set_index('uid');z.obs['dpt_pseudotime']=c.loc[z.obs_names,'K15_R1'].to_numpy();root=pd.read_csv(O/'tables/graph_QA.tsv',sep='\t').query('track == @track and variant == "K15_R1"').iloc[0].root_cell;z.uns['iroot']=int(np.flatnonzero(z.obs_names==root)[0]);z.write_h5ad(f)
# Display-only state backbone: weighted kNN edges between DPT quantile bins, no inferred lineage claim.
from scipy.sparse.csgraph import minimum_spanning_tree
back=[]
for track in ['Naive_ActCD4','Treg']:
 z=ad.read_h5ad(O/'objects'/(track+'_geometry.h5ad'));c=pd.read_csv(O/'tables'/(track+'_cell_trajectory.tsv.gz'),sep='\t').set_index('uid').loc[z.obs_names];bins=pd.qcut(c.K15_R1,8,labels=False,duplicates='drop').to_numpy();nb=int(bins.max()+1);mat=z.obsp['connectivities'].tocoo();wgt=np.zeros((nb,nb))
 for i,j,v in zip(bins[mat.row],bins[mat.col],mat.data):
  if i!=j:wgt[i,j]+=v
 ct=np.bincount(bins);wgt=wgt/np.sqrt(ct[:,None]*ct[None,:]);dist=np.where(wgt>0,1/(wgt+1e-12),0);tree=minimum_spanning_tree(dist).tocoo();cent=[]
 for k in range(nb):
  ids=np.flatnonzero(bins==k);xy=c.iloc[ids][['UMAP1','UMAP2']].to_numpy();mid=np.median(xy,axis=0);cent.append(xy[np.argmin(((xy-mid)**2).sum(1))])
 for i,j in zip(tree.row,tree.col):
  if i>j:i,j=j,i
  back.append({'track':track,'bin_a':int(i),'bin_b':int(j),'x1':cent[i][0],'y1':cent[i][1],'x2':cent[j][0],'y2':cent[j][1]})
w(pd.DataFrame(back),'trajectory_display_backbone.tsv')
counts=[]
for sc,z in r.groupby('scope'):
 counts.append({'scope':sc,'trajectory_tested':len(z),'trajectory_P_lt05':int(z.trajectory_nominal.sum()),'time_D28D7_tested':int(z.time_P.notna().sum()),'intersection_all':int(z.primary_intersection.sum()),'primary_apoptosis':int((z.primary_death_candidate&z.apoptosis_GO).sum()),'primary_pyroptosis':int((z.primary_death_candidate&z.pyroptosis_GO).sum()),'primary_death':int(z.primary_death_candidate.sum()),'secondary_death':int(z.secondary_death_candidate.sum())})
w(pd.DataFrame(counts),'intersection_counts.tsv');summary={'status':'success','primary_death_genes':main.gene.unique().tolist(),'primary_gene_state_pairs':len(main),'secondary_death_genes':sorted(set(secondary.gene)),'primary_pyroptosis_genes':main.loc[main.pyroptosis_GO,'gene'].tolist(),'primary_core_execution_genes':main.loc[main.core_execution_panel,'gene'].tolist(),'time_primary_tested_genes':len(t),'time_primary_nominal_genes':int((t.time_P<.05).sum()),'counts':counts};(O/'logs/intersection_complete.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2));print(json.dumps(summary,ensure_ascii=False));print(main[['gene','scope','function_tier','mean_slope','P','FDR','time_logFC','time_P','time_FDR','variants_P_lt05','direction_concordant','time_min10_P']].to_string(index=False));print('MODULE_TRAJECTORY',pd.read_csv(O/'tables/module_trajectory_effects.tsv',sep='\t').query('variant == "K15_R1"').to_string(index=False))