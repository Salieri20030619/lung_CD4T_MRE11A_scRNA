import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2',MKL_NUM_THREADS='2',NUMBA_NUM_THREADS='2')
from pathlib import Path
import sys,json,hashlib,itertools,importlib.metadata
import numpy as np,pandas as pd,scanpy as sc,anndata as ad
from scipy import sparse
from scipy.io import mmread
from scipy.sparse.csgraph import connected_components
from scipy.stats import false_discovery_control
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=Path(sys.argv[1]);M=P/'02_analysis/mechanism_step1_20260927_152218_s4V5eF';V=P/'11_reports/customer_sequence_20261015_v2'
S=['Saline_1','Saline_2','Saline_3','Silica_1','Silica_2','Silica_3'];SCOPES=['Core','Naive','ActCD4','Treg'];MODS=['PyroExec','ApoExec','Inflam','MitoISR','MitoROS','Mitophagy','Autophagy_GO','OXPHOS']
for d in ['scores','trajectory','time','intersection','provenance','code']:(O/d).mkdir(parents=True,exist_ok=True)
spec={'task':'independent transcriptome-wide screening','no_existing_candidate_list':True,'no_disease_DE_gate':True,'scores':'recomputed matched-expression scores; no named candidate excluded from control selection','primary_graph':'K15_R1','sensitivity_graphs':['K15_R2','K30_R1','K30_R2'],'trajectory_test':'sample-wise depth-adjusted linear+quadratic curve coefficients; joint exact 64 sign flips of six sample vectors','trajectory_unit':'independent mouse','expression_filter':'CPM >=1 in at least 3 of 6 pseudobulk libraries within scope','common_support':'intersection of six sample observed DPT ranges; no extrapolated sample tails','intersection':'trajectory joint raw P<0.05 AND independent real-time omnibus raw P<0.05','direction':'annotated after intersection, not a mandatory gate','BH':'full evaluable gene family per scope/graph or time test','time_comparisons':['D7_D0','D28_D0','D56_D0','D28_D7','D56_D28'],'no_novelty_or_existing_target_priority':True}
(O/'provenance/specification.json').write_text(json.dumps(spec,ensure_ascii=False,indent=2),encoding='utf-8')
versions={x:importlib.metadata.version(x) for x in ['scanpy','anndata','numpy','pandas','scipy']};(O/'provenance/python_versions.json').write_text(json.dumps(versions,indent=2),encoding='utf-8');print('ENVIRONMENT',sys.executable,versions,flush=True)
manifest=[]
def source(f):
 f=Path(f);h=hashlib.sha256()
 with f.open('rb') as q:
  for block in iter(lambda:q.read(2**20),b''):h.update(block)
 manifest.append({'file':str(f),'sha256':h.hexdigest(),'bytes':f.stat().st_size});return f
def rd(f,**kw):return pd.read_csv(source(f),sep='\t',**kw)
def wr(d,sub,n):d.to_csv(O/sub/n,sep='\t',index=False,encoding='utf-8',na_rep='NA')
def bh(p):
 p=np.asarray(p,float);q=np.full(p.shape,np.nan);ok=np.isfinite(p)
 if ok.any():q[ok]=false_discovery_control(p[ok],method='bh')
 return q
md=rd(M/'comm_input/CD4_cells.tsv');gn=rd(M/'comm_input/CD4_genes.tsv');x=mmread(source(M/'comm_input/CD4_counts.mtx.gz')).T.tocsr();genes=pd.Index(gn.canonical)
assert x.shape==(len(md),len(genes)) and md.cell.is_unique and genes.is_unique and set(md['sample'])==set(S)
lib=md.library.to_numpy();assert np.allclose(x.data,np.rint(x.data)) and np.all(np.asarray(x.sum(1)).ravel()<=lib+1e-6)
y=x.multiply((1e4/lib)[:,None]).tocsr();y.data=np.log1p(y.data)
defs=rd(M/'tables/module_definition.tsv');go=rd(V/'inputs/GO0006914_autophagy_annotations.tsv')
alias={'MRE11A':'MRE11','ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','PARK2':'PRKN'}
canon=lambda g:alias.get(str(g).upper(),str(g).upper())
sets={m:defs.loc[defs.module==m,'gene'].tolist() for m in MODS if m!='Autophagy_GO'};sets['Autophagy_GO']=sorted(set(go.SYMBOL.dropna().map(canon)))
reference=np.mean([np.asarray(y[md['sample'].eq(s).to_numpy()].mean(0)).ravel() for s in S],axis=0)
excluded=set(defs.gene)|set(sets['Autophagy_GO']);bg=set('SFTPA1 SFTPA2 SFTPB SFTPC SFTPD SCGB1A1 SCGB3A1 SCGB3A2 ABCA3 SLC34A2 NAPSA LYZ LYZ1 LYZ2 S100A8 S100A9 CSF1R LST1 TYROBP FCER1G'.split())
pool=(reference>0)&np.array([g not in excluded and g not in bg and not g.startswith(('MT-','RPL','RPS','MRPL','MRPS','IGH','IGK','IGL')) for g in genes]);edges=np.unique(np.quantile(reference[reference>0],np.linspace(0,1,26)));bins=np.digitize(reference,edges[1:-1]);cs=md[['cell','sample','subtype','group','library']].copy();coverage=[];members=[];controls=[]
for mod in MODS:
 ix=genes.get_indexer(sets[mod]);ix=ix[ix>=0];observed=int((np.asarray(x[:,ix].sum(0)).ravel()>0).sum());cc=set()
 for j in ix:
  candidates=np.flatnonzero(pool&(bins==bins[j]));candidates=candidates if len(candidates)>=20 else np.flatnonzero(pool);order=np.argsort(abs(reference[candidates]-reference[j]),kind='stable');cc.update(candidates[order[:20]])
 cc=np.array(sorted(cc));ok=len(ix)>=3 and observed>=3 and len(ix)/len(sets[mod])>=.6 and len(cc)>=20
 cs[mod]=np.asarray(y[:,ix].mean(1)).ravel()-np.asarray(y[:,cc].mean(1)).ravel() if ok else np.nan
 coverage.append(dict(module=mod,expected=len(sets[mod]),measured=len(ix),expressed=observed,controls=len(cc),eligible=ok));members.extend(dict(module=mod,gene=g,measured=g in genes) for g in sets[mod]);controls.extend(dict(module=mod,gene=genes[j],reference_mean=reference[j]) for j in cc)
wr(pd.DataFrame(coverage),'scores','module_coverage.tsv');wr(pd.DataFrame(members),'scores','module_members.tsv');wr(pd.DataFrame(controls),'scores','control_genes.tsv');wr(cs,'scores','cell_scores.tsv.gz')
assert np.isfinite(cs[MODS].to_numpy()).all(), 'Non-evaluable score panel'
sample=[]
for scope in SCOPES:
 b=cs if scope=='Core' else cs[cs.subtype==scope]
 for s,z in b.groupby('sample'):
  sample.append(dict(scope=scope,sample=s,group=z.group.iloc[0],cells=len(z),**{mod:z[mod].mean() for mod in MODS}))
sample=pd.DataFrame(sample);wr(sample,'scores','sample_scores.tsv')
alloc=np.array(list(itertools.combinations(range(6),3)));tests=[];orderrows=[]
for scope in SCOPES:
 b=sample[sample.scope==scope].set_index('sample').loc[S]
 for mod in MODS:
  vals=b[mod].to_numpy();delta=vals[3:].mean()-vals[:3].mean();null=np.array([vals[t].mean()-vals[np.setdiff1d(range(6),t)].mean() for t in alloc]);pv=float((abs(null)>=abs(delta)-1e-12).mean()) if np.ptp(vals)>1e-12 else np.nan
  loo=[np.delete(vals[3:],j).mean()-vals[:3].mean() for j in range(3)]+[vals[3:].mean()-np.delete(vals[:3],j).mean() for j in range(3)]
  tests.append(dict(scope=scope,module=mod,n_control=3,n_case=3,control_mean=vals[:3].mean(),case_mean=vals[3:].mean(),difference=delta,P=pv,loo_same_direction=float(np.mean(np.sign(loo)==np.sign(delta))),family=('Core' if scope=='Core' else 'Subtypes')+'_'+('Death' if mod in MODS[:3] else 'MitoAuto')))
tests=pd.DataFrame(tests);tests['BH']=tests.groupby('family').P.transform(bh);wr(tests,'scores','disease_score_tests.tsv')
for mod in MODS:
 b=tests[(tests.module==mod)&tests.scope.ne('Core')].copy();b['case_score_rank']=b.case_mean.rank(ascending=False,method='min').astype(int);b['increase_rank']=b.difference.rank(ascending=False,method='min').astype(int);b['absolute_change_rank']=b.difference.abs().rank(ascending=False,method='min').astype(int);orderrows.append(b)
wr(pd.concat(orderrows),'scores','subtype_rankings.tsv');print('SCORES_COMPLETE',flush=True)

# Rebuild graphs on all Core CD4, including Treg; no original gene shortlist used.
h=ad.read_h5ad(source(P/'02_analysis/GSE314943/stage5_20260905_222457_a2dApT/objects/CD4_core.h5ad'));assert set(h.obs_names)==set(md.cell)
h=h[md.cell].copy();pcs=np.asarray(h.obsm['X_pca_harmony'])[:,:30];gidx=genes.get_indexer(['TCF7','LEF1','CCR7']);assert (gidx>=0).all();rootscore=np.asarray(y[:,gidx].mean(1)).ravel();rootinfo=[];ptrows=[];graphinfo=[];trtests=[];slopes=[];curves=[];support=[];signs=np.array(list(itertools.product([-1.,1.],repeat=6)));grid=np.linspace(-1,1,41);basis=np.column_stack([grid,grid**2]);basis-=basis.mean(0);metric=basis.T@basis/len(grid)
for k in [15,30]:
 b=ad.AnnData(sparse.csr_matrix((len(md),1)));b.obsm['X_pca_harmony']=pcs
 sc.pp.neighbors(b,n_neighbors=k,n_pcs=30,use_rep='X_pca_harmony',random_state=42)
 nc,comp=connected_components(b.obsp['connectivities'],directed=False);eligible=[]
 for c in range(nc):
  take=comp==c;z=md.loc[take]
  if set(z['sample'])==set(S) and (z.subtype=='Naive').any():eligible.append((int(take.sum()),c))
 if not eligible:graphinfo.append(dict(k=k,status='NO_COMPONENT_WITH_ALL_SAMPLES'));continue
 _,component=max(eligible);global_ids=np.flatnonzero(comp==component);bb=b[global_ids].copy();mm=md.iloc[global_ids].reset_index(drop=True);rs=rootscore[global_ids];cand=np.flatnonzero(mm.subtype.eq('Naive').to_numpy());cand=cand[rs[cand]>=np.quantile(rs[cand],.75)];pc=pcs[global_ids];center=pc[cand].mean(0);rank=cand[np.argsort(((pc[cand]-center)**2).sum(1),kind='stable')];first=int(rank[0]);second=next((int(v) for v in rank if mm.iloc[v]['sample']!=mm.iloc[first]['sample']),None);assert second is not None
 sc.tl.diffmap(bb,n_comps=15,random_state=42)
 graphinfo.append(dict(k=k,status='OK',components=nc,cells=len(mm),excluded_cells=len(md)-len(mm),state_counts=mm.subtype.value_counts().to_json()))
 for ri,root in enumerate([first,second],1):
  run=f'K{k}_R{ri}';bb.uns['iroot']=root;sc.tl.dpt(bb,n_dcs=10,n_branchings=0);pt=bb.obs.dpt_pseudotime.to_numpy();assert np.isfinite(pt).all();rootinfo.append(dict(run=run,cell=mm.iloc[root].cell,sample=mm.iloc[root]['sample'],RootScore=rs[root]));ptrows.extend(dict(run=run,cell=row.cell,sample=row['sample'],scope=row.subtype,pseudotime=float(pt[j])) for j,(_,row) in enumerate(mm.iterrows()))
  for scope in SCOPES:
   base=np.ones(len(mm),bool) if scope=='Core' else mm.subtype.eq(scope).to_numpy();sids=[np.flatnonzero(base&mm['sample'].eq(s).to_numpy()) for s in S]
   if any(len(v)==0 for v in sids):support.append(dict(run=run,scope=scope,status='MISSING_SAMPLE'));continue
   lo=max(pt[v].min() for v in sids);hi=min(pt[v].max() for v in sids)
   if hi-lo<1e-6:support.append(dict(run=run,scope=scope,status='NO_COMMON_SUPPORT'));continue
   ids=np.flatnonzero(base);gi=global_ids[ids];pb=np.vstack([np.asarray(x[global_ids[v]].sum(0)).ravel() for v in sids]);libs=np.array([lib[global_ids[v]].sum() for v in sids]);cpm=pb/libs[:,None]*1e6;keep=(cpm>=1).sum(0)>=3;geneix=np.flatnonzero(keep);names=genes[geneix].tolist()+MODS;values=np.column_stack([y[global_ids][:,geneix].toarray(),cs.iloc[global_ids][MODS].to_numpy()]);coefs=[];valid=True;counts=[]
   for s,v in zip(S,sids):
    vv=v[(pt[v]>=lo)&(pt[v]<=hi)];counts.append(len(vv));z=(pt[vv]-lo)/(hi-lo)*2-1;depth=np.log1p(mm.iloc[vv].library.to_numpy());depth-=depth.mean();A=np.column_stack([np.ones(len(vv)),depth,z,z*z]);rankA=np.linalg.matrix_rank(A)
    if len(vv)<=rankA or rankA<4:valid=False;break
    fit=np.linalg.lstsq(A,values[vv],rcond=None)[0];coefs.append(fit[2:4].T)
   support.append(dict(run=run,scope=scope,status='OK' if valid else 'INSUFFICIENT_WITHIN_SAMPLE_SUPPORT',common_pt_min=lo,common_pt_max=hi,cells_by_sample=json.dumps(dict(zip(S,counts))),genes_expression_eligible=len(geneix)))
   if not valid:continue
   coefs=np.asarray(coefs);mean=coefs.mean(0);stat=np.einsum('gi,ij,gj->g',mean,metric,mean);nullcoef=np.einsum('ps,sgd->pgd',signs,coefs)/6;null=np.einsum('pgi,ij,pgj->pg',nullcoef,metric,nullcoef);pv=(null>=stat[None]-1e-12).mean(0);pv[np.max(abs(coefs),axis=(0,2))<1e-10]=np.nan
   q=np.r_[bh(pv[:len(geneix)]),bh(pv[len(geneix):])];shape=[]
   for j,name in enumerate(names):
    a1,a2=mean[j];turn=-a1/(2*a2) if abs(a2)>1e-12 else np.nan;nonmono=np.isfinite(turn) and -1<turn<1;pred=basis@mean[j];delta=2*a1
    trtests.append(dict(run=run,scope=scope,feature=name,type='gene' if j<len(geneix) else 'program',n_samples=6,beta_linear=a1,beta_quadratic=a2,endpoint_difference=delta,curve_range=float(np.ptp(pred)),shape=('valley' if a2>0 else 'peak') if nonmono else ('increasing' if delta>0 else 'decreasing'),P=pv[j],BH=q[j],family=f'{run}_{scope}_'+('all_expression_eligible_genes' if j<len(geneix) else 'all8_programs')))
    if run=='K15_R1':
     for s,coef in zip(S,coefs[:,j]):slopes.append(dict(scope=scope,sample=s,feature=name,type='gene' if j<len(geneix) else 'program',beta_linear=coef[0],beta_quadratic=coef[1]))
   if run=='K15_R1':
    for s,v in zip(S,sids):
     vv=v[(pt[v]>=lo)&(pt[v]<=hi)];bins_q=pd.qcut(pd.Series(pt[vv]).rank(method='first'),4,labels=False).to_numpy()
     for qv in range(4):
      chosen=vv[bins_q==qv]
      for mod in MODS:curves.append(dict(scope=scope,sample=s,group=mm.iloc[chosen[0]].group,bin=qv,cells=len(chosen),mean_pt=float(pt[chosen].mean()),module=mod,mean_score=float(cs.iloc[global_ids[chosen]][mod].mean())))
   print('TRAJECTORY_DONE',run,scope,len(geneix),'nominal',int((pv[:len(geneix)]<.05).sum()),flush=True)
wr(pd.DataFrame(trtests),'trajectory','all_gene_and_program_tests.tsv.gz');wr(pd.DataFrame(slopes),'trajectory','primary_sample_coefficients.tsv.gz');wr(pd.DataFrame(curves),'trajectory','primary_program_sample_curves.tsv');wr(pd.DataFrame(ptrows),'trajectory','pseudotime_membership.tsv.gz');wr(pd.DataFrame(rootinfo),'trajectory','root_definitions.tsv');wr(pd.DataFrame(graphinfo),'trajectory','graph_coverage.tsv');wr(pd.DataFrame(support),'trajectory','common_support.tsv');wr(pd.DataFrame(manifest).drop_duplicates('file'),'provenance','input_manifest.tsv')
(O/'provenance/score_trajectory_complete.json').write_text(json.dumps({'status':'COMPLETE','samples':6,'cells':len(md),'trajectory_rows':len(trtests),'new_graphs':2,'root_settings':len(rootinfo)},indent=2),encoding='utf-8');print('SCORE_TRAJECTORY_COMPLETE',flush=True)
