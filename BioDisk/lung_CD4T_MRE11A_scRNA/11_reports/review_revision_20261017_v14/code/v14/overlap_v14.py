import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
from pathlib import Path
import numpy as np,pandas as pd,itertools,json
from scipy.io import mmread
from scipy.stats import rankdata,false_discovery_control
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/review_revision_20261017_v14';Q=O/'supplement_v14';V=P/'11_reports/integrated_macrophage_CD4_20261015_v4';M=P/'02_analysis/mechanism_step1_20260927_152218_s4V5eF';BASE=P/'11_reports/methods_complete_20261016_v13/supplement'
rd=lambda p:pd.read_csv(p,sep='\t')
S=['Saline_1','Saline_2','Saline_3','Silica_1','Silica_2','Silica_3'];mods=['PyroExec','Inflam','ApoExec','OXPHOS','MitoISR','MitoROS','Mitophagy'];scopes=['Core','Naive','ActCD4','Treg'];PERM=np.array(list(itertools.permutations(range(6))));WITHIN=np.array([a+b for a in itertools.permutations(range(3)) for b in itertools.permutations(range(3,6))])
plan={'reason':'observed overlap between matched controls and frozen regulon targets','posthoc':True,'TFs':['STAT1','IKZF2'],'seed':17,'modules':mods,'scopes':scopes,'rule':'remove focal regulon targets and TF itself from existing matched controls; keep original module members, normalisation and cells; no rematching','family':'all seven modules per TF/scope/model; both pooled and condition-adjusted; fixed composition separate','no_original_replacement':True}
(O/'provenance/control_overlap_sensitivity_plan.json').write_text(json.dumps(plan,indent=2),encoding='utf-8')
def test(x,y,adjust):
 x=rankdata(x);y=rankdata(y)
 for ix in ([np.arange(3),np.arange(3,6)] if adjust else [np.arange(6)]):x[ix]-=x[ix].mean();y[ix]-=y[ix].mean()
 den=np.linalg.norm(x)*np.linalg.norm(y)
 if den<1e-12:return np.nan,np.nan
 rho=x@y/den;null=x[WITHIN if adjust else PERM]@y/den;return rho,float((abs(null)>=abs(rho)-1e-12).mean())
def bh(p):
 p=np.asarray(p,float);q=np.full(len(p),np.nan);ok=np.isfinite(p)
 if ok.any():q[ok]=false_discovery_control(p[ok])
 return q
md=rd(M/'comm_input/CD4_cells.tsv');gn=rd(M/'comm_input/CD4_genes.tsv');x=mmread(M/'comm_input/CD4_counts.mtx.gz').T.tocsr();genes=pd.Index(gn.canonical);logx=x.multiply((1e4/md.library.to_numpy())[:,None]).tocsr();logx.data=np.log1p(logx.data)
mod=rd(BASE/'S20_all_fixed_module_genes.tsv');control=rd(BASE/'S20_all_fixed_module_controls.tsv.gz');control=control[control.dataset=='GSE314943'];edge=rd(Q/'S25_reconstructed_focal_regulon_members.tsv');edge=edge[edge.seed==17];act=rd(Q/'S22_TF_all_regulon_sample_activity.tsv.gz');act=act[act.seed==17];old=rd(V/'01_phenotype/sample_program_scores.tsv');old=old[old['filter']=='All'];rows=[];controls=[];parity=[]
for tf in ['STAT1','IKZF2']:
 targets=set(edge.loc[edge.tf_canonical==tf,'target_canonical'])|{tf}
 for module in mods:
  members=set(mod.loc[mod.module==module,'gene']);c=set(control.loc[control.module==module,'gene']);clean=c-targets
  im=genes.get_indexer(sorted(members));im=im[im>=0];ic=genes.get_indexer(sorted(c));ic=ic[ic>=0];iz=genes.get_indexer(sorted(clean));iz=iz[iz>=0];assert len(iz)>=20
  score0=np.asarray(logx[:,im].mean(axis=1)-logx[:,ic].mean(axis=1)).ravel();score1=np.asarray(logx[:,im].mean(axis=1)-logx[:,iz].mean(axis=1)).ravel()
  for g in sorted(c):controls.append(dict(tf=tf,module=module,gene=g,retained=g in clean,removed_reason='regulon_or_TF_overlap' if g not in clean else ''))
  for sc in scopes:
   for s in S:
    mask=md['sample'].eq(s).to_numpy()&(np.ones(len(md),bool) if sc=='Core' else md.subtype.eq(sc).to_numpy());v0=float(score0[mask].mean());v1=float(score1[mask].mean());expected=float(old[(old.scope==sc)&(old['sample']==s)][module].iloc[0]);parity.append(abs(v0-expected));rows.append(dict(tf=tf,module=module,scope=sc,sample=s,original_score=v0,clean_control_score=v1,controls_original=len(ic),controls_retained=len(iz),cells=int(mask.sum())))
assert max(parity)<1e-6,max(parity)
z=pd.DataFrame(rows);z.to_csv(Q/'S27_control_overlap_sample_scores.tsv',sep='\t',index=False);pd.DataFrame(controls).to_csv(Q/'S27_control_overlap_members.tsv',sep='\t',index=False)
# Fixed composition uses the original six-mouse average proportions.
props=pd.crosstab(md['sample'],md.subtype,normalize='index').loc[S];weights=props.mean();fixed=[]
for (tf,module),b in z.groupby(['tf','module']):
 for s in S:
  c=b[(b['sample']==s)&b.scope.isin(weights.index)].set_index('scope');fixed.append(dict(tf=tf,module=module,scope='CoreFixed',sample=s,original_score=float(c.original_score@weights.loc[c.index]),clean_control_score=float(c.clean_control_score@weights.loc[c.index])))
z=pd.concat([z,pd.DataFrame(fixed)],ignore_index=True);tests=[]
for (tf,sc),b in z.groupby(['tf','scope']):
 if sc=='CoreFixed':
  aa=act[(act.tf==tf)&act.scope.isin(weights.index)].pivot(index='sample',columns='scope',values='mean_auc').loc[S];av=(aa@weights.loc[aa.columns]).to_numpy()
 else:av=act[(act.tf==tf)&(act.scope==sc)].set_index('sample').loc[S].mean_auc.to_numpy()
 for module,c in b.groupby('module'):
  c=c.set_index('sample').loc[S]
  for adj in [False,True]:
   for variant in ['original_score','clean_control_score']:
    rho,p=test(av,c[variant].to_numpy(),adj);tests.append(dict(tf=tf,scope=sc,module=module,variant=variant,model='ConditionAdjusted' if adj else 'Pooled',rho=rho,P=p,n_samples=6,family=f'{tf}_{sc}_{adj}_{variant}_7modules'))
t=pd.DataFrame(tests);t['BH']=t.groupby('family').P.transform(bh);t.to_csv(Q/'S27_control_overlap_associations.tsv',sep='\t',index=False);z.to_csv(Q/'S27_control_overlap_sample_scores.tsv',sep='\t',index=False)
(O/'provenance/control_overlap_validation.json').write_text(json.dumps({'original_score_max_abs_error':max(parity),'original_scores_reproduced':len(parity),'sample_count':6,'new_rows':len(t),'status':'PASS'},indent=2),encoding='utf-8')
print('OVERLAP_SENSITIVITY',t[((t.tf=='STAT1')&(t.module=='MitoROS')|(t.tf=='IKZF2')&(t.module=='ApoExec'))&t.scope.isin(['Core','CoreFixed','ActCD4','Treg'])].to_string(index=False))