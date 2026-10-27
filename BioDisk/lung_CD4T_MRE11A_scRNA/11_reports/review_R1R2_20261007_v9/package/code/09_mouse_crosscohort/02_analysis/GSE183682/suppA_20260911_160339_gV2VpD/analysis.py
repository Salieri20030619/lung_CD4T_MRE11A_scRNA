from pathlib import Path
import os,sys,json,tarfile
import numpy as np
import pandas as pd
from scipy import sparse
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(os.environ['PROJECT_ROOT']);O=Path(sys.argv[1]);O.mkdir(parents=True,exist_ok=True)
source=R/'02_analysis/screen2_CD4_target_wz1lfx/GSE183682/CD4_core.h5ad'
a=sc.read_h5ad(source);assert a.n_obs==1413 and a.obs_names.is_unique
mapping={'0':'NaiveCD4','3':'ActCD4','5':'Treg'}
assert set(a.obs.F08.astype(str))==set(mapping)
a.obs['subtype']=a.obs.F08.astype(str).map(mapping)
meta={'GSM5567501':('NS',7),'GSM5567502':('Silica',7),'GSM5567503':('NS',56),'GSM5567504':('Silica',56)}
assert set(a.obs['sample'].astype(str))==set(meta)
C=sparse.csr_matrix(a.layers['counts']);assert np.isfinite(C.data).all() and np.all(C.data>=0) and np.allclose(C.data,np.rint(C.data),rtol=0,atol=1e-6)
C=C.astype(np.int64);ix=np.flatnonzero(a.var_names=='MRE11A');assert len(ix)==1
v=np.asarray(C[:,ix].sum(1)).ravel();lib=np.asarray(C.sum(1)).ravel();assert (lib>0).all();norm=v/lib*10000
obs=a.obs.copy();obs['target_count']=v;obs['library_counts']=lib;obs['CP10k']=norm;obs.to_csv(O/'cell_target.tsv.gz',sep='\t')
rows=[];overall=[];pb={};types=list(mapping.values())
for sample,(group,day) in meta.items():
 s=np.asarray(obs['sample']==sample);n=int(s.sum())
 overall.append(dict(sample=sample,group=group,day=day,n_cells=n,target_positive=int((v[s]>0).sum()),target_counts=int(v[s].sum()),target_CPM=1e6*v[s].sum()/lib[s].sum(),mean_CP10k=float(norm[s].mean()),pct_detected=100*np.mean(v[s]>0)))
 for label in types:
  z=s&np.asarray(obs.subtype==label);nz=int(z.sum());assert nz>0,f'Missing subtype: {sample}:{label}'
  vec=np.asarray(C[z].sum(0)).ravel();pb[sample+'__'+label]=vec
  rows.append(dict(sample=sample,group=group,day=day,subtype=label,n_cells=nz,cell_fraction=nz/n,target_positive=int((v[z]>0).sum()),target_counts=int(v[z].sum()),library_counts=int(lib[z].sum()),target_CPM=1e6*v[z].sum()/lib[z].sum(),mean_CP10k=float(norm[z].mean()),pct_detected=100*np.mean(v[z]>0),median_library=float(np.median(lib[z])),small_subtype=nz<20))
t=pd.DataFrame(rows);alltab=pd.DataFrame(overall);t.to_csv(O/'subtype_target.tsv',sep='\t',index=False);alltab.to_csv(O/'overall_target.tsv',sep='\t',index=False)
pd.DataFrame(pb,index=a.var_names).to_csv(O/'subtype_pseudobulk.tsv.gz',sep='\t')
contrasts=[];parts=[];totalparts=[]
for day in [7,56]:
 x=t[(t.day==day)&(t.group=='NS')].set_index('subtype').loc[types];y=t[(t.day==day)&(t.group=='Silica')].set_index('subtype').loc[types]
 for label in types:
  xx=x.loc[label];yy=y.loc[label]
  contrasts.append(dict(day=day,subtype=label,NS_n=int(xx.n_cells),Silica_n=int(yy.n_cells),NS_CPM=xx.target_CPM,Silica_CPM=yy.target_CPM,CPM_difference=yy.target_CPM-xx.target_CPM,CPM_ratio=yy.target_CPM/xx.target_CPM if xx.target_CPM>0 else np.nan,mean_CP10k_difference=yy.mean_CP10k-xx.mean_CP10k,detection_percentage_point_difference=yy.pct_detected-xx.pct_detected))
 # Symmetric decomposition of the change in cell-weighted mean CP10k.
 # Delta = sum((p1-p0)*(m1+m0)/2) + sum((m1-m0)*(p1+p0)/2).
 p0=x.cell_fraction.to_numpy();p1=y.cell_fraction.to_numpy();m0=x.mean_CP10k.to_numpy();m1=y.mean_CP10k.to_numpy();w=(p0+p1)/2
 comp=(p1-p0)*(m1+m0)/2;within=(m1-m0)*w;delta=float(p1@m1-p0@m0)
 assert np.isclose(comp.sum()+within.sum(),delta,rtol=1e-9,atol=1e-12)
 for i,label in enumerate(types):parts.append(dict(day=day,subtype=label,standard_weight=w[i],composition_contribution=comp[i],within_contribution=within[i]))
 totalparts.append(dict(day=day,observed_NS_mean_CP10k=p0@m0,observed_Silica_mean_CP10k=p1@m1,observed_difference=delta,composition_component=comp.sum(),within_component=within.sum(),standardized_NS_mean_CP10k=w@m0,standardized_Silica_mean_CP10k=w@m1,standardized_difference=w@(m1-m0)))
pd.DataFrame(contrasts).to_csv(O/'subtype_contrasts.tsv',sep='\t',index=False);pd.DataFrame(parts).to_csv(O/'decomposition_by_subtype.tsv',sep='\t',index=False);pd.DataFrame(totalparts).to_csv(O/'decomposition_summary.tsv',sep='\t',index=False)
colors={'NS':'#2469A0','Silica':'#C64B36'}
for metric,ylabel in [('target_CPM','Mre11a CPM'),('pct_detected','Detected (%)'),('cell_fraction','Cell fraction')]:
 f,axs=plt.subplots(1,3,figsize=(12,3.7),sharey=False)
 for ax,label in zip(axs,types):
  q=t[t.subtype==label]
  for group in ['NS','Silica']:
   z=q[q.group==group].sort_values('day');ax.plot(z.day,z[metric],marker='o',color=colors[group],label=group)
   for _,r in z.iterrows():ax.annotate('n='+str(int(r.n_cells)),(r.day,r[metric]),xytext=(3,5 if group=='NS' else -12),textcoords='offset points',fontsize=7)
  ax.set(title=label,xlabel='Day',ylabel=ylabel,xticks=[7,56]);ax.legend(frameon=False,fontsize=8);ax.spines[['top','right']].set_visible(False)
 f.tight_layout();f.savefig(O/(metric+'_subtype.png'),dpi=350);f.savefig(O/(metric+'_subtype.pdf'));plt.close(f)
f,axs=plt.subplots(1,2,figsize=(8,3.7),sharey=True)
for ax,r in zip(axs,totalparts):
 vals=[r['observed_difference'],r['composition_component'],r['within_component']];ax.axhline(0,color='grey',lw=.7);ax.scatter(range(3),vals,s=55,color=['#444444','#C78A28','#276A9A']);ax.vlines(range(3),0,vals,color=['#444444','#C78A28','#276A9A']);ax.set(xticks=range(3),xticklabels=['Total','Comp','Within'],title='Day '+str(r['day']),ylabel='Delta CP10k');ax.spines[['top','right']].set_visible(False)
f.tight_layout();f.savefig(O/'decomposition.png',dpi=350);f.savefig(O/'decomposition.pdf');plt.close(f)
(O/'captions.txt').write_text('Subtype panels: x=actual sampling day (7 or 56), y=panel-specific measurement. Lines connect independent libraries at different timepoints; they are not longitudinal animal trajectories. Blue NS=saline control; red Silica=silica-exposed. n=selected subtype cell count, not biological sample size. NaiveCD4=naive-like CD4; ActCD4=activated-like CD4; Treg=regulatory T. CPM=pseudobulk target count divided by all measured selected-subtype counts times 1e6. Detected (%)=percentage of cells with positive raw Mre11a counts. Cell fraction=the subtype fraction among the three selected core CD4 subtypes. Decomposition: y=change in mean per-cell CP10k (Silica minus NS); Total=observed difference; Comp=symmetric composition component; Within=symmetric within-subtype expression component. CP10k=per-cell target count divided by measured library counts times 10000. Equal midpoint subtype weights within each timepoint are used for standardized estimates. These are descriptive components, not causal effects. CPM and mean CP10k use different weighting and are reported separately. No disease p-values are computed because each time/condition has one expression library.\n')
(O/'run_summary.json').write_text(json.dumps({'input':str(source),'cells':a.n_obs,'target_positive':int((v>0).sum()),'target_counts':int(v.sum()),'subtype_mapping':mapping,'analysis_logic':'Fixed scope; within-time comparisons; cell fractions sum to one; exact decomposition identity checked; no inferential p-values','status':'PASS'},indent=2))
with tarfile.open(R/'tmp/GSE183682_suppA_return.tar.gz','w:gz') as tar:
 for p in O.iterdir():
  if p.is_file() and 'pseudobulk' not in p.name:tar.add(p,arcname=p.name)
print('GSE183682_SUPPA_PASS',flush=True)
