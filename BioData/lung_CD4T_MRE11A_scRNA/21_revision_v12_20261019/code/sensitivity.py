from pathlib import Path
import os,json,itertools,hashlib
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2',NUMBA_NUM_THREADS='2')
import numpy as np,pandas as pd,anndata as ad,scanpy as sc,harmonypy as hm
from scipy import sparse,stats
from scipy.sparse.csgraph import connected_components
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');V=P/'11_reports/death_priority_20261018_v31';B=P/'11_reports/three_cohort_four_axis_20261017_v16'
O=Path('/mnt/BioData/bio/wangmingyuan/lung_CD4T_MRE11A_scRNA/21_revision_v12_20261019');(O/'tables').mkdir(parents=True,exist_ok=True);(O/'objects').mkdir(exist_ok=True)
def wr(x,n):pd.DataFrame(x).to_csv(O/'tables'/n,sep='\t',index=False)
seed=20261018;genes=['HIF1A','IL2RA','ARL6IP1'];markers=['CCR7','SELL','TCF7','LEF1','FOXP3','CTLA4','CD3D','CD3E','TRAC']
prof=[]
for ds in ['GSE314943','GSE183682']:
 stem=ds+'__Core';m=pd.read_csv(B/'inputs'/f'{stem}.cells.tsv.gz',sep='\t');g=pd.read_csv(B/'inputs'/f'{stem}.genes.txt',header=None)[0].astype(str);x=sparse.load_npz(B/'inputs'/f'{stem}.cellcounts.npz').tocsr()
 m['uid']=ds+'|'+m.cell.astype(str);m['dataset']=ds;m['sample_id']=ds+'|'+m['sample'].astype(str)
 for gene in genes+markers:
  idx=np.flatnonzero(g==gene);m['G_'+gene]=np.log1p(x[:,idx[0]].toarray().ravel()/m.library.to_numpy()*1e4) if len(idx) else np.nan
 prof.append(m)
prof=pd.concat(prof,ignore_index=True).set_index('uid');qa=pd.read_csv(V/'tables/graph_QA.tsv',sep='\t')
original=pd.read_csv(V/'tables/gene_trajectory_all.tsv.gz',sep='\t')
results=[];donors=[];cells=[];rootrows=[];graphrows=[];changes=[]
def evaluate(d,scope,gene,pt,scenario,variant):
 z=d[d.subtype==scope].copy();raw=pt[d.subtype.eq(scope).to_numpy()];tt=(raw-raw.min())/(raw.max()-raw.min());z['t']=tt
 slopes=[];ns=[];origslopes=[];maxdelta=[]
 for sam,t in z.groupby('sample_id'):
  dep=np.log1p(t.library.to_numpy());xx=t.t.to_numpy();yy=t['G_'+gene].to_numpy();C=np.column_stack([np.ones(len(t)),dep]);rx=xx-C@np.linalg.lstsq(C,xx,rcond=None)[0];lever=rx**2/(rx@rx);imax=int(lever.argmax());take=np.ones(len(t),bool)
  if scenario=='drop_top_leverage_per_donor':take[imax]=False
  if scenario=='drop_Naive_Silica2_influential' and scope=='Naive' and sam=='GSE314943|Silica_2':take[t.uid.eq('GSE314943|Silica_2:TCCGATCCAGAATCGG-1').to_numpy()]=False
  def fit(mask):
   D=np.column_stack([np.ones(mask.sum()),dep[mask]]);r=xx[mask]-D@np.linalg.lstsq(D,xx[mask],rcond=None)[0]
   return float(r@yy[mask]/(r@r))
  full=fit(np.ones(len(t),bool));val=fit(take);slopes.append(val);ns.append(int(take.sum()));origslopes.append(full)
  donors.append(dict(scope=scope,gene=gene,variant=variant,scenario=scenario,sample_id=sam,cells_before=len(t),cells_after=int(take.sum()),slope=val,full_slope=full,max_partial_leverage=lever[imax],removed_cells=';'.join(t.loc[~take,'uid'])))
 arr=np.asarray(slopes);sg=np.asarray(list(itertools.product([-1,1],repeat=len(arr))));p=float(np.mean(np.abs(sg@arr/len(arr))>=abs(arr.mean())-1e-12))
 results.append(dict(scope=scope,gene=gene,scenario=scenario,variant=variant,n=len(arr),cells=sum(ns),min_cells=min(ns),mean_slope=arr.mean(),P=p,same_direction=int((np.sign(arr)==np.sign(arr.mean())).sum()),FDR=np.nan,FDR_note='targeted sensitivity only; no claim of full-gene FDR'))
 return results[-1]
for track in ['Naive_ActCD4','Treg']:
 f=V/'objects'/f'{track}_geometry.h5ad';oldhash=hashlib.sha256(f.read_bytes()).hexdigest()
 a=ad.read_h5ad(f);d=pd.read_csv(V/'tables'/f'{track}_cell_trajectory.tsv.gz',sep='\t')
 for g in genes+markers:d['G_'+g]=prof.loc[d.uid,'G_'+g].to_numpy()
 for _,r in qa[qa.track==track].iterrows():
  row=d[d.uid==r.root_cell].iloc[0];rootrows.append(dict(track=track,variant=r.variant,**{k:row[k] for k in ['uid','sample_id','subtype','group']+['G_'+g for g in markers]}))
 for scope,g in ([('Naive','HIF1A'),('Naive','ARL6IP1')] if track=='Naive_ActCD4' else [('Treg','IL2RA')]):
  for scenario in ['baseline_reproduction','drop_top_leverage_per_donor']+(['drop_Naive_Silica2_influential'] if scope=='Naive' else []):
   rr=evaluate(d,scope,g,d.K15_R1.to_numpy(),scenario,'K15_R1')
   if scenario=='baseline_reproduction':
    old=original[(original.scope==scope)&(original.gene==g)&(original.variant=='K15_R1')].iloc[0];assert np.isclose(rr['mean_slope'],old.mean_slope,rtol=2e-5) and np.isclose(rr['P'],old.P)
 cells.append(d)
 if track=='Treg':
  # Isolate IL2RA influence: retain all other 1,999 scaled features and original root cells.
  b=a[:,a.var_names!='IL2RA'].copy()
  sc.tl.pca(b,n_comps=30,svd_solver='arpack',random_state=seed)
  ho=hm.run_harmony(b.obsm['X_pca'],b.obs,'dataset',random_state=seed,max_iter_harmony=30,verbose=False);Z=ho.Z_corr;Z=Z.T if Z.shape[0]!=b.n_obs else Z;b.obsm['X_harmony']=Z
  for k in [15,30]:
   sc.pp.neighbors(b,n_neighbors=k,use_rep='X_harmony',method='gauss',random_state=seed);nc,_=connected_components(b.obsp['connectivities'],directed=False);sc.tl.diffmap(b,n_comps=15,random_state=seed)
   for rno in [1,2]:
    variant=f'K{k}_R{rno}';root=qa[(qa.track==track)&(qa.variant==variant)].root_cell.iloc[0];b.uns['iroot']=int(np.flatnonzero(b.obs_names==root)[0]);sc.tl.dpt(b,n_dcs=10);pt=b.obs.dpt_pseudotime.to_numpy().copy();assert np.isfinite(pt).all()
    d['noIL2RA_'+variant]=pt
    rr=evaluate(d,'Treg','IL2RA',pt,'geometry_without_IL2RA_fixed_roots',variant)
    graphrows.append(dict(variant=variant,cells=b.n_obs,features=b.n_vars,components=nc,root=root,finite=int(np.isfinite(pt).sum()),rho_vs_original=stats.spearmanr(pt,d[variant]).statistic))
  b.write_h5ad(O/'objects/Treg_noIL2RA_sensitivity_last_setting.h5ad')
  wr(d,'Treg_cell_profiles.tsv.gz')
 else:wr(d,'Naive_ActCD4_cell_profiles.tsv.gz')
 assert hashlib.sha256(f.read_bytes()).hexdigest()==oldhash
wr(results,'targeted_trajectory_sensitivity.tsv');wr(donors,'sensitivity_donor_slopes.tsv');wr(rootrows,'root_marker_expression.tsv');wr(graphrows,'noIL2RA_geometry_QA.tsv')
protocol=dict(purpose='targeted correction sensitivity',primary_analysis_unchanged=True,cell_influence='fixed original DPT scale; remove top partial-DPT-leverage cell per donor; separately remove identified Silica2 cell only',geometry='Treg only; remove IL2RA from existing scaled feature matrix, retain1999; rerun PCA/Harmony/neighbors/diffmap/DPT with original fixed root identifiers, K15/K30 and R1/R2',statistics='10 independent mouse slopes;1024 sign flips; targeted diagnostic raw P only; primary full-family FDR unchanged',seed=seed)
(O/'protocol.json').write_text(json.dumps(protocol,indent=2),encoding='utf-8')
(O/'complete.json').write_text(json.dumps(dict(status='success',results=results,source_objects_unchanged=True),indent=2),encoding='utf-8')
print(pd.DataFrame(results).to_string(index=False),flush=True);print('COMPLETE',flush=True)
