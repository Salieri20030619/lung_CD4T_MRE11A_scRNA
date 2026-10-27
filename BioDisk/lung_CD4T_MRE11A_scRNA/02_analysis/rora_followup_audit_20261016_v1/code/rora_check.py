import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2',MKL_NUM_THREADS='2')
from pathlib import Path
import sys,itertools,json,hashlib
import pandas as pd,numpy as np
from scipy.io import mmread
from scipy.stats import rankdata,false_discovery_control
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');I=P/'02_analysis/independent_screen_20261016_v1';M=P/'02_analysis/mechanism_step1_20260927_152218_s4V5eF';O=Path(sys.argv[1]);O.mkdir(parents=True,exist_ok=True)
S=['Saline_1','Saline_2','Saline_3','Silica_1','Silica_2','Silica_3'];SC=['Core','Naive','ActCD4','Treg'];mods=['PyroExec','ApoExec','Inflam','MitoISR','MitoROS','Mitophagy','Autophagy_GO','OXPHOS']
panels={'activation':['CD44','ICOS','CD38','TNFRSF4','TNFRSF9'],'tissue':['ITGAE','ITGA1','CXCR6','RGS1'],'naive':['TCF7','LEF1','CCR7','SELL']}
spec={'purpose':'readonly follow-up audit; original results preserved','activation':panels['activation'],'tissue_proxy':panels['tissue'],'panel_type':'unweighted mean log1pCP10k proxy, not validated state classifier','primary_association':'within-mouse partial Pearson correlation controlling logdepth, activation and tissue proxies and DPT linear/quadratic; sample-level Fisher z exact64 signflip','association_primary_family':'Autophagy_GO and Mitophagy across three subtypes (6 tests); Core separate','secondary_family':'other six programs across three subtypes (18 tests); Core separate','no_cell_level_P':True,'adjustment_sensitivity':'original scores; target-excluded scores; target and covariate marker excluded scores','disease_test':'20 two-sided permutations of 3v3 sample means; adjusted means from within-sample centered nuisance regression, equal sample weights','important':'No non-silicosis disease comparator; specificity cannot be established; adjustment is descriptive sensitivity, not causal identification.'}
(O/'specification.json').write_text(json.dumps(spec,indent=2),encoding='utf-8')
manifest=[]
def src(p):
 p=Path(p);h=hashlib.sha256(p.read_bytes()).hexdigest();manifest.append({'path':str(p),'sha256':h});return p
def rd(p):return pd.read_csv(src(p),sep='\t')
def wr(d,n):d.to_csv(O/n,sep='\t',index=False,encoding='utf-8',na_rep='NA')
def bh(p):
 p=np.array(p,float);q=np.full(len(p),np.nan);ok=np.isfinite(p)
 if ok.any():q[ok]=false_discovery_control(p[ok])
 return q
md=rd(M/'comm_input/CD4_cells.tsv');gn=rd(M/'comm_input/CD4_genes.tsv');x=mmread(src(M/'comm_input/CD4_counts.mtx.gz')).T.tocsr();genes=pd.Index(gn.canonical);lib=md.library.to_numpy();y=x.multiply((1e4/lib)[:,None]).tocsr();y.data=np.log1p(y.data)
assert genes.is_unique and set(md['sample'])==set(S)
def expr(names):
 ix=genes.get_indexer(names);assert (ix>=0).all(),names
 return y[:,ix].toarray()
d=md.copy();d['RORA']=expr(['RORA']).ravel();d['RORA_count']=x[:,genes.get_loc('RORA')].toarray().ravel();d['logdepth']=np.log1p(lib)
for n,gs in panels.items():d[n]=expr(gs).mean(1)
pt=rd(I/'trajectory/pseudotime_membership.tsv.gz').query("run=='K15_R1'")[['cell','pseudotime']];d=d.merge(pt,on='cell',validate='one_to_one',how='left');assert d.pseudotime.notna().all();d['pt2']=d.pseudotime**2
scores=rd(I/'scores/cell_scores.tsv.gz');assert (scores.cell==d.cell).all();members=rd(I/'scores/module_members.tsv');controls=rd(I/'scores/control_genes.tsv');coverage=[]
for mode,excluded in [('original',set()),('target_excluded',{'RORA'}),('marker_excluded',{'RORA'}|set(sum(panels.values(),[])))]:
 for mod in mods:
  mg=members.loc[(members.module==mod)&members.measured,'gene'].tolist();cg=controls.loc[controls.module==mod,'gene'].tolist();mg2=[v for v in mg if v not in excluded];cg2=[v for v in cg if v not in excluded]
  vals=expr(mg2).mean(1)-expr(cg2).mean(1);d[mode+'__'+mod]=vals
  if mode=='original':assert np.allclose(vals,scores[mod],atol=1e-10)
  coverage.append(dict(mode=mode,module=mod,member_count=len(mg2),control_count=len(cg2),removed_members=';'.join(sorted(set(mg)-set(mg2))),removed_controls=';'.join(sorted(set(cg)-set(cg2)))))
wr(pd.DataFrame(coverage),'score_overlap_audit.tsv')
# Mechanism genes are assessed descriptively for detection, never assigned cell-level P.
mechanism=['CAV3','TP53','DRAM1','PINK1','PRKN','BNIP3','BNIP3L','FUNDC1','SQSTM1','OPTN','TBK1','ULK1','BECN1','MAP1LC3B','RORA','RORC']
detect=[]
for scope in SC:
 take=np.ones(len(d),bool) if scope=='Core' else d.subtype.eq(scope).to_numpy()
 for s in S:
  ids=np.flatnonzero(take&d['sample'].eq(s).to_numpy())
  for g in mechanism:
   if g not in genes:detect.append(dict(scope=scope,sample=s,gene=g,measured=False));continue
   a=x[ids,genes.get_loc(g)].toarray().ravel();detect.append(dict(scope=scope,sample=s,gene=g,measured=True,cells=len(ids),positive_cells=int((a>0).sum()),fraction_positive=float((a>0).mean()),UMI=float(a.sum()),CPM=float(a.sum()/lib[ids].sum()*1e6)))
wr(pd.DataFrame(detect),'mechanism_gene_detection.tsv')
sample=[];disease=[];alloc=list(itertools.combinations(range(6),3));coeff=[]
def group_test(v):
 delta=np.mean(v[3:])-np.mean(v[:3]);null=[np.mean(v[list(a)])-np.mean(v[np.setdiff1d(np.arange(6),a)]) for a in alloc];return delta,float(np.mean(np.abs(null)>=abs(delta)-1e-12))
for scope in SC:
 b=d.copy() if scope=='Core' else d[d.subtype==scope].copy()
 for s,z in b.groupby('sample'):
  sample.append(dict(scope=scope,sample=s,group=z.group.iloc[0],cells=len(z),RORA=z.RORA.mean(),RORA_detection=z.RORA_count.gt(0).mean(),RORA_CPM=z.RORA_count.sum()/z.library.sum()*1e6,activation=z.activation.mean(),tissue=z.tissue.mean(),naive=z.naive.mean()))
 for model,covs in [('unadjusted',[]),('depth_activation_tissue',['logdepth','activation','tissue']),('plus_pseudotime',['logdepth','activation','tissue','pseudotime','pt2'])]:
  if scope=='Core':covs=covs+(['isAct','isTreg'] if model!='unadjusted' else [])
  b['isAct']=b.subtype.eq('ActCD4').astype(float);b['isTreg']=b.subtype.eq('Treg').astype(float)
  if covs:
   centered=b[covs]-b.groupby('sample')[covs].transform('mean');yc=b.RORA-b.groupby('sample').RORA.transform('mean');w=1/np.sqrt(b.groupby('sample').cell.transform('count').to_numpy());coef=np.linalg.lstsq(centered.to_numpy()*w[:,None],yc.to_numpy()*w,rcond=None)[0];ref=b.groupby('sample')[covs].mean().mean().to_numpy();adj=b.RORA.to_numpy()-(b[covs].to_numpy()-ref)@coef
   coeff.extend(dict(scope=scope,model=model,covariate=n,beta=v) for n,v in zip(covs,coef))
  else:adj=b.RORA.to_numpy()
  b['adjusted_RORA']=adj;means=b.groupby('sample').adjusted_RORA.mean().loc[S].to_numpy();delta,p=group_test(means);loo=[np.mean(np.delete(means[3:],j))-np.mean(means[:3]) for j in range(3)]+[np.mean(means[3:])-np.mean(np.delete(means[:3],j)) for j in range(3)]
  disease.append(dict(scope=scope,model=model,control_mean=np.mean(means[:3]),case_mean=np.mean(means[3:]),difference=delta,P=p,loo_sign_agreement=np.mean(np.sign(loo)==np.sign(delta)),sample_values=json.dumps(dict(zip(S,means)))))
wr(pd.DataFrame(sample),'sample_expression.tsv');disease=pd.DataFrame(disease);disease['BH']=disease.groupby('model').P.transform(bh);wr(disease,'disease_comparisons.tsv');wr(pd.DataFrame(coeff),'nuisance_coefficients.tsv')
cor=[];support=rd(I/'trajectory/common_support.tsv').query("run=='K15_R1'");signs=np.array(list(itertools.product([-1,1],repeat=6)))
for scope in SC:
 b=d.copy() if scope=='Core' else d[d.subtype==scope].copy();sup=support[support.scope==scope].iloc[0];b=b[b.pseudotime.between(sup.common_pt_min,sup.common_pt_max)]
 for s,z in b.groupby('sample'):
  for mode in ['original','target_excluded','marker_excluded']:
   for model,covs in [('depth',['logdepth']),('full',['logdepth','activation','tissue','pseudotime','pt2'])]:
    if scope=='Core':
     z=z.copy();z['isAct']=z.subtype.eq('ActCD4').astype(float);z['isTreg']=z.subtype.eq('Treg').astype(float);covs=covs+['isAct','isTreg']
    A=np.column_stack([np.ones(len(z)),z[covs].to_numpy()]);rank=np.linalg.matrix_rank(A);valid=len(z)>rank+2
    for mod in mods+['activation','tissue']:
     if mod in covs:continue
     a=z.RORA.to_numpy();v=z[mode+'__'+mod if mod in mods else mod].to_numpy();res_a=a-A@np.linalg.lstsq(A,a,rcond=None)[0];res_v=v-A@np.linalg.lstsq(A,v,rcond=None)[0]
     r=float(np.corrcoef(res_a,res_v)[0,1]) if valid and np.std(res_a)>1e-10 and np.std(res_v)>1e-10 else np.nan
     cor.append(dict(scope=scope,sample=s,group=z.group.iloc[0],score_mode=mode,model=model,module=mod,cells=len(z),design_rank=rank,partial_r=r))
cor=pd.DataFrame(cor);wr(cor,'within_sample_correlations.tsv');summ=[]
for keys,z in cor.groupby(['scope','score_mode','model','module']):
 scope,mode,model,mod=keys;a=z.set_index('sample').reindex(S).partial_r.to_numpy();ok=np.isfinite(a).all()
 if ok:
  f=np.arctanh(np.clip(a,-.999999,.999999));null=signs@f/6;pv=np.mean(abs(null)>=abs(f.mean())-1e-12);effect=np.tanh(f.mean())
 else:pv=effect=np.nan
 fam=('Core' if scope=='Core' else 'Subtypes')+'_'+('Autophagy2' if mod in ['Autophagy_GO','Mitophagy'] else 'State2' if mod in panels else 'Other6')
 summ.append(dict(scope=scope,score_mode=mode,model=model,module=mod,n_samples=int(np.isfinite(a).sum()),mean_Fisher_r=effect,P=pv,positive_samples=int((a>0).sum()),negative_samples=int((a<0).sum()),min_r=np.nanmin(a),max_r=np.nanmax(a),saline_mean_r=np.tanh(np.arctanh(np.clip(a[:3],-.999999,.999999)).mean()),silica_mean_r=np.tanh(np.arctanh(np.clip(a[3:],-.999999,.999999)).mean()),family=fam))
summ=pd.DataFrame(summ);summ['BH']=summ.groupby(['score_mode','model','family']).P.transform(bh);wr(summ,'correlation_summary.tsv');wr(d,'audit_cell_data.tsv.gz');wr(pd.DataFrame(manifest).drop_duplicates('path'),'input_manifest.tsv')
(O/'COMPLETE.json').write_text(json.dumps({'status':'COMPLETE','cells':len(d),'samples':6,'correlation_tests':len(summ)},indent=2),encoding='utf-8')
print(disease[['scope','model','difference','P','BH']].to_string(index=False));print(summ[(summ.score_mode=='target_excluded')&(summ.model=='full')&(summ.module.isin(['Autophagy_GO','Mitophagy']))].to_string(index=False));print('RORA_AUDIT_COMPLETE')
