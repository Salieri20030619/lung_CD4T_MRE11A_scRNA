import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2',MKL_NUM_THREADS='2',NUMBA_NUM_THREADS='2')
from pathlib import Path
import json,itertools,hashlib,time
import numpy as np,pandas as pd
from scipy import sparse,stats
from scipy.sparse.csgraph import connected_components
import scanpy as sc,anndata as ad,harmonypy as hm
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');B=P/'11_reports/three_cohort_four_axis_20261017_v16';O=P/'11_reports/death_priority_20261018_v31';seed=20261018
alias={'MRE11A':'MRE11','ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','PARK2':'PRKN'}
def canon(g):g=str(g).upper();return alias.get(g,g)
def write(a,n):a.to_csv(O/'tables'/n,sep='\t',index=False)
def bh(a,pcol,fcol,groups):
 a[fcol]=np.nan
 for key,z in a.groupby(groups,dropna=False):
  ix=z.index[z[pcol].notna()]
  if len(ix):a.loc[ix,fcol]=stats.false_discovery_control(a.loc[ix,pcol].to_numpy())
 return a
inputs={};audit=[]
for ds,scope in [('GSE314943','Core'),('GSE183682','Core'),('GSE274372','TLike')]:
 stem=ds+'__'+scope;files=[B/'inputs'/(stem+s) for s in ['.cellcounts.npz','.cells.tsv.gz','.genes.txt']];x=sparse.load_npz(files[0]).astype(np.float64);m=pd.read_csv(files[1],sep='\t');g=pd.read_csv(files[2],header=None)[0].astype(str).to_numpy();m['dataset']=ds;m['uid']=ds+'|'+m.cell.astype(str);m['sample_id']=ds+'|'+m['sample'].astype(str);inputs[ds]=(x,m,g)
 for f in files:audit.append({'file':str(f),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
common3=sorted(set.intersection(*[set(v[2]) for v in inputs.values()]));common2=sorted(set(inputs['GSE314943'][2])&set(inputs['GSE183682'][2]))
(O/'inputs/manifest.json').write_text(json.dumps(audit,indent=2));go=pd.read_csv(O/'tables/GO_priority_members.tsv',sep='\t');go['gene']=go.gene.map(canon);oldmods=pd.read_csv(P/'11_reports/append_analysis_20261017_v18/tables/module_member_audit.tsv',sep='\t');oldmods['gene']=oldmods.gene.map(canon)
sets={mod:sorted(set(oldmods.loc[oldmods.module.eq(mod),'gene'])) for mod in ['ApoExec','PyroExec','MitoISR','Mitophagy','TAct','Naive']};sets['Autophagy']=sorted(set(go.loc[go.GOALL.eq('GO:0006914'),'gene'].dropna()))
coverage=[];module_rows=[]
for mod,genes in sets.items():
 mapped=sorted(set(genes)&set(common3));coverage.append({'module':mod,'defined':len(genes),'common_mapped':len(mapped),'coverage':len(mapped)/len(genes)});assert len(mapped)>=3 and len(mapped)>=.5*len(genes),(mod,len(mapped),len(genes));sets[mod]=mapped
 module_rows.extend({'module':mod,'gene':g} for g in mapped)
write(pd.DataFrame(coverage),'module_coverage.tsv');write(pd.DataFrame(module_rows),'common_module_members.tsv')
# Sample scores use common member sets and complete original library denominators.
scores=[]
for ds,(x,m,g) in inputs.items():
 ix={v:i for i,v in enumerate(g)}
 for variant in ['All','BG1']:
  for (sample,state),mm in m.groupby(['sample','subtype']):
   if variant=='BG1':mm=mm[mm.background_fraction<=.01]
   if not len(mm):continue
   ids=mm.index.to_numpy();counts=np.asarray(x[ids].sum(0)).ravel();lib=mm.library.sum();v=np.log1p(counts/lib*1e4);r={'dataset':ds,'sample':sample,'scope':state,'group':mm.group.iloc[0],'day':int(mm.day.iloc[0]),'cells':len(mm),'variant':variant,'library':lib}
   for mod,genes in sets.items():r[mod]=v[[ix[t] for t in genes]].mean()
   scores.append(r)
scores=pd.DataFrame(scores);write(scores,'sample_module_scores.tsv');primarymods=['ApoExec','PyroExec','MitoISR','Mitophagy','Autophagy']
def design(m):return np.column_stack([np.ones(len(m)),pd.get_dummies(m.dataset.astype(str)+'_'+m.day.astype(str),drop_first=True,dtype=float).to_numpy()])
def perm_effect(m,Y):
 C=design(m);M=np.eye(len(m))-C@np.linalg.pinv(C);g=m.group.eq('Silica').to_numpy(float);x=M@g;ef=x@Y/(x@x);labels=m.dataset.astype(str)+'_'+m.day.astype(str);groups=[np.flatnonzero(labels.eq(k)) for k in sorted(labels.unique())];choices=[list(itertools.combinations(a,int(g[a].sum()))) for a in groups];null=[]
 for parts in itertools.product(*choices):
  v=np.zeros(len(m))
  for ids in parts:v[list(ids)]=1
  xx=M@v;null.append(xx@Y/(xx@xx))
 null=np.asarray(null);pv=np.mean(np.abs(null)>=np.abs(ef)-1e-12,axis=0);return ef,pv,len(null)
effects=[];ranks=[];paired=[]
for state in ['Naive','ActCD4','Treg']:
 z=scores[scores.variant.eq('All')&scores.scope.eq(state)&scores.dataset.ne('GSE274372')].sort_values(['dataset','sample']).reset_index(drop=True);assert len(z)==10;ef,pv,n=perm_effect(z,z[primarymods].to_numpy())
 for j,mod in enumerate(primarymods):effects.append({'scope':state,'module':mod,'effect':ef[j],'P':pv[j],'allocations':n,'n':10,'family':mod+'|disease_3states'});ranks.append({'scope':state,'module':mod,'Saline_mean':z.loc[z.group.eq('Saline'),mod].mean(),'Silica_mean':z.loc[z.group.eq('Silica'),mod].mean(),'Silica_median':z.loc[z.group.eq('Silica'),mod].median()})
effects=bh(pd.DataFrame(effects),'P','FDR',['family']);write(effects,'module_disease_effects.tsv');write(pd.DataFrame(ranks),'module_state_ranking.tsv')
for mod in primarymods:
 z=scores[scores.variant.eq('All')&scores.dataset.ne('GSE274372')].pivot(index=['dataset','sample','group'],columns='scope',values=mod)
 for aa,bb in itertools.combinations(['Naive','ActCD4','Treg'],2):
  dif=(z[aa]-z[bb]).to_numpy();sg=np.asarray(list(itertools.product([-1,1],repeat=len(dif))));pv=np.mean(np.abs(sg@dif/len(dif))>=abs(dif.mean())-1e-12);paired.append({'module':mod,'state_a':aa,'state_b':bb,'paired_mean_difference':dif.mean(),'P':pv,'n':len(dif),'permutations':len(sg),'family':mod+'|paired_state_3pairs'})
write(bh(pd.DataFrame(paired),'P','FDR',['family']),'module_paired_state_contrasts.tsv')
# Time module trends use day factor as a separate model, no paired animals across days.
time_results=[]
z=scores[scores.variant.eq('All')&scores.dataset.eq('GSE274372')].reset_index(drop=True);day=pd.get_dummies(z.day.astype(str),dtype=float).reindex(columns=['0','7','28','56']);C=day.to_numpy()
for mod in primarymods:
 y=z[mod].to_numpy();beta=np.linalg.lstsq(C,y,rcond=None)[0];res=y-C@beta;s2=res@res/(len(y)-4);rss0=np.sum((y-y.mean())**2);f=((rss0-res@res)/3)/s2;ct=np.array([0,-1,1,0]);se=np.sqrt(s2*(ct@np.linalg.inv(C.T@C)@ct));tt=(ct@beta)/se
 time_results.extend([{'module':mod,'contrast':'time_omnibus','effect':np.nan,'statistic':f,'P':stats.f.sf(f,3,len(y)-4),'family':'time_modules_omnibus'},{'module':mod,'contrast':'D28_vs_D7','effect':ct@beta,'statistic':tt,'P':2*stats.t.sf(abs(tt),len(y)-4),'family':'time_modules_D28_vs_D7'}])
write(bh(pd.DataFrame(time_results),'P','FDR',['family']),'time_module_effects.tsv')
print('SCORES_DONE',pd.DataFrame(coverage).to_dict('records'),flush=True)
# Combined expression geometry; no old candidate restrictions.
mat=[];obs=[]
for ds in ['GSE314943','GSE183682']:
 x,m,g=inputs[ds];idx=pd.Index(g).get_indexer(common2);mat.append(x[:,idx]);obs.append(m)
x=sparse.vstack(mat,format='csr');m=pd.concat(obs,ignore_index=True);norm=sparse.diags(1e4/m.library.to_numpy())@x;log=norm.copy().tocsr();log.data=np.log1p(log.data);geneindex={g:i for i,g in enumerate(common2)}
blockgenes=set.union(*[set(sets[k]) for k in primarymods]);scope_results=[];slopes_records=[];graphqa=[];curves=[];gene_tables=[]
for track,states in [('Naive_ActCD4',['Naive','ActCD4']),('Treg',['Treg'])]:
 take=m.subtype.isin(states).to_numpy();ids=np.flatnonzero(take);mm=m.loc[take].reset_index(drop=True);X=log[ids].tocsr();aa=ad.AnnData(X=X.copy(),obs=mm.set_index('uid'),var=pd.DataFrame(index=common2));sc.pp.highly_variable_genes(aa,n_top_genes=2500,flavor='seurat',batch_key='dataset')
 eligible=np.array([g not in blockgenes and not g.startswith(('MT-','RPL','RPS','IGH','IGK','IGL')) for g in common2]);hvg=aa.var.highly_variable.to_numpy()&eligible
 # Maintain predeclared max2000 after excluding module genes, using normalized dispersion only.
 order=np.argsort(-np.nan_to_num(aa.var.dispersions_norm.to_numpy(),nan=-np.inf));selected=[i for i in order if hvg[i]][:2000];aa=aa[:,selected].copy();sc.pp.scale(aa,max_value=10);sc.tl.pca(aa,n_comps=min(30,len(selected)-1),svd_solver='arpack',random_state=seed)
 ho=hm.run_harmony(aa.obsm['X_pca'],aa.obs,'dataset',random_state=seed,max_iter_harmony=30,verbose=False);Z=ho.Z_corr;Z=Z.T if Z.shape[0]!=len(mm) else Z;aa.obsm['X_harmony']=Z
 tact=np.asarray(X[:,[geneindex[g] for g in sets['TAct'] if g in geneindex]].mean(axis=1)).ravel();eligible_roots=mm.group.eq('Saline').to_numpy()&(mm.subtype.eq('Naive').to_numpy() if track=='Naive_ActCD4' else tact<=np.quantile(tact[mm.group.eq('Saline')],.25))
 roots=[]
 for ds in ['GSE314943','GSE183682']:
  choices=np.flatnonzero(eligible_roots&mm.dataset.eq(ds).to_numpy());assert len(choices)>0;center=Z[choices].mean(0);roots.append(int(choices[np.argmin(((Z[choices]-center)**2).sum(1))]))
 print('GRAPH',track,len(mm),'HVG',len(selected),'ROOTS',roots,flush=True)
 variants={}
 for k in [15,30]:
  sc.pp.neighbors(aa,n_neighbors=k,use_rep='X_harmony',method='gauss',random_state=seed);cc,cl=connected_components(aa.obsp['connectivities'],directed=False);sc.tl.diffmap(aa,n_comps=15,random_state=seed)
  if k==15:sc.tl.umap(aa,random_state=seed)
  for rno,root in enumerate(roots):
   aa.uns['iroot']=root;sc.tl.dpt(aa,n_dcs=10);v=aa.obs.dpt_pseudotime.to_numpy().copy();name=f'K{k}_R{rno+1}';variants[name]=v;graphqa.append({'track':track,'variant':name,'cells':len(v),'finite_cells':int(np.isfinite(v).sum()),'components':cc,'root_cell':mm.uid.iloc[root],'root_dataset':mm.dataset.iloc[root]})
  if k==15:
   # Save primary graph distances/coordinates before changing neighborhood sensitivity.
   aa.write_h5ad(O/'objects'/(track+'_geometry.h5ad'))
 coords=mm.copy();coords['track']=track;coords['UMAP1']=aa.obsm['X_umap'][:,0];coords['UMAP2']=aa.obsm['X_umap'][:,1]
 for v,pseudo in variants.items():coords[v]=pseudo
 write(coords,track+'_cell_trajectory.tsv.gz')
 for state in states:
  allstate=mm.subtype.eq(state).to_numpy();samples=sorted(mm.loc[allstate,'sample_id'].unique());assert len(samples)==10
  # Pre-test gene expression eligibility independent of trajectory P-values.
  cps=[];smeta=[]
  for sam in samples:
   pos=np.flatnonzero(allstate&mm.sample_id.eq(sam).to_numpy());m0=mm.iloc[pos[0]];cps.append(np.asarray(x[ids[pos]].sum(0)).ravel()/mm.library.iloc[pos].sum()*1e6);smeta.append({'sample_id':sam,'sample':m0['sample'],'dataset':m0.dataset,'day':m0.day,'group':m0.group,'scope':state})
  cp=np.asarray(cps);sm=pd.DataFrame(smeta);good=(cp>=1).sum(0)>=3
  for ds in ['GSE314943','GSE183682']:good&=(cp[sm.dataset.eq(ds)]>=1).sum(0)>=1
  gids=np.flatnonzero(good);gn=np.asarray(common2)[gids];sg=np.asarray(list(itertools.product([-1,1],repeat=10)),dtype=float)
  for variant,pt0 in variants.items():
   finite=np.isfinite(pt0)&allstate;pt=pt0.copy();v=pt[finite];pt[finite]=(v-v.min())/(v.max()-v.min());slope=[];met=[]
   for sam in samples:
    pos=np.flatnonzero(finite&mm.sample_id.eq(sam).to_numpy());assert len(pos)>=10,(track,state,variant,sam,len(pos));tt=pt[pos];dep=np.log1p(mm.library.iloc[pos].to_numpy());C=np.column_stack([np.ones(len(pos)),dep]);rx=tt-C@np.linalg.lstsq(C,tt,rcond=None)[0];den=rx@rx;assert den>1e-10
    svec=np.asarray(X[pos][:,gids].T@rx).ravel()/den;slope.append(svec)
    for mod in primarymods:
     inds=[geneindex[g] for g in sets[mod] if g in geneindex];yy=np.asarray(X[pos][:,inds].mean(1)).ravel();ms=float(rx@yy/den);slopes_records.append({'track':track,'scope':state,'sample_id':sam,'variant':variant,'module':mod,'slope':ms,'cells':len(pos)})
     if variant=='K15_R1':
      bins=np.minimum((tt*4).astype(int),3)
      for bn in range(4):
       jj=bins==bn
       if jj.any():curves.append({'scope':state,'sample_id':sam,'group':mm.group.iloc[pos[0]],'module':mod,'bin':bn,'DPT':tt[jj].mean(),'score':yy[jj].mean(),'cells':int(jj.sum())})
   slope=np.asarray(slope);mean=slope.mean(0);ps=np.empty(len(gids))
   for start in range(0,len(gids),1000):
    stop=min(start+1000,len(gids));null=sg@slope[:,start:stop]/10;ps[start:stop]=(np.abs(null)>=np.abs(mean[start:stop])-1e-12).mean(0)
   ef,pdisease,alloc=perm_effect(sm,slope)
   rr=pd.DataFrame({'gene':gn,'scope':state,'track':track,'variant':variant,'mean_slope':mean,'P':ps,'FDR':stats.false_discovery_control(ps),'same_direction_samples':(np.sign(slope)==np.sign(mean)[None,:]).sum(0),'n':10,'disease_slope_difference':ef,'disease_P':pdisease,'disease_FDR':stats.false_discovery_control(pdisease),'family':state+'|'+variant+'|gene_trajectory','disease_family':state+'|'+variant+'|disease_slope'})
   gene_tables.append(rr)
   if variant=='K15_R1':
    out=pd.DataFrame(slope.T,columns=samples);out.insert(0,'gene',gn);write(out,state+'_gene_sample_slopes.tsv.gz')
   print('TRENDS',state,variant,len(rr),int((rr.P<.05).sum()),flush=True)
write(pd.concat(gene_tables,ignore_index=True),'gene_trajectory_all.tsv.gz');write(pd.DataFrame(slopes_records),'module_sample_slopes.tsv');write(pd.DataFrame(curves),'module_trajectory_curves.tsv');write(pd.DataFrame(graphqa),'graph_QA.tsv')
mods=pd.DataFrame(slopes_records);rows=[]
for (scope,variant,mod),z in mods.groupby(['scope','variant','module']):
 v=z.slope.to_numpy();sg=np.asarray(list(itertools.product([-1,1],repeat=len(v))));pv=(np.abs(sg@v/len(v))>=abs(v.mean())-1e-12).mean();rows.append({'scope':scope,'variant':variant,'module':mod,'mean_slope':v.mean(),'P':pv,'n':len(v),'same_sign':int((np.sign(v)==np.sign(v.mean())).sum()),'family':scope+'|'+variant+'|module_trajectory'})
write(bh(pd.DataFrame(rows),'P','FDR',['family']),'module_trajectory_effects.tsv')
(O/'logs/python_analysis_complete.json').write_text(json.dumps({'status':'success','cells_joint':len(m),'common_genes':len(common2),'trajectory_rows':sum(len(t) for t in gene_tables),'graphs':len(graphqa)},indent=2));print('PYTHON_ANALYSIS_COMPLETE',flush=True)