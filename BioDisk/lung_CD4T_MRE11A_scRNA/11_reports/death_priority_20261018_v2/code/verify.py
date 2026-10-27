import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
from pathlib import Path
import pandas as pd,numpy as np,json,itertools,math,hashlib
from scipy import stats
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/death_priority_20261018_v2')
checks=[]
def check(name,ok,detail=''):
 checks.append(dict(check=name,passed=bool(ok),detail=str(detail)))
 if not ok:raise AssertionError(name+': '+str(detail))
def bh(p):
 p=np.asarray(p,float);n=len(p);order=np.argsort(p);v=np.minimum.accumulate((p[order]*n/np.arange(1,n+1))[::-1])[::-1];out=np.empty(n);out[order]=np.minimum(v,1);return out
for name,pc,qc in [('joint_upstream_associations.tsv','P','FDR'),('all_associations.tsv','P','FDR'),('program_effects.tsv','P','FDR'),('TF_associations.tsv','P','FDR'),('TF_disease_effects.tsv','P','FDR'),('joint_composition.tsv','P','FDR'),('joint_LR_disease_effects.tsv','P','FDR'),('joint_ORA_all.tsv','P','FDR'),('joint_GSEA_all.tsv','pval','padj'),('Core_composition_adjusted_GSEA.tsv','pval','padj'),('joint_sender_DE_all.tsv','PValue','FDR')]:
 d=pd.read_csv(O/'tables'/name,sep='\t');err=0
 for fam,z in d.groupby('family'):
  z=z[z[pc].notna()]
  if len(z):err=max(err,float(np.max(abs(bh(z[pc])-z[qc]))))
 check(name+' full family BH',err<1e-10,err)
 if 'cohort' in d and 'n' in d:check(name+' joint n10',d.loc[(d.cohort=='Joint')&d[pc].notna(),'n'].eq(10).all())
score=pd.read_csv(O/'tables/sample_scores.tsv',sep='\t');a=pd.read_csv(O/'tables/all_associations.tsv',sep='\t')
z=a[(a.cohort=='Joint')&a.kind.str.endswith('primary')&a.scope.ne('Core')&a.variant.isin(['Raw','BG1'])&a.P.notna()]
cache={};records=[]
for r in z.itertuples():
 d=score[(score.variant==r.variant)&(score.scope==r.scope)&score.dataset.ne('GSE274372')].sort_values(['dataset','sample']).reset_index(drop=True)
 x=stats.rankdata(d[r.x]);y=stats.rankdata(d[r.y]);strata=d.dataset+'_'+d.day.astype(str);C=np.column_stack([np.ones(len(d)),pd.get_dummies(strata,drop_first=True,dtype=float)])
 labels=strata.copy()
 if 'Condition' in r.model:C=np.column_stack([C,(d.group=='Silica').astype(float)]);labels=labels+'_'+d.group
 Q=np.eye(len(d))-C@np.linalg.pinv(C);rx=Q@x;ry=Q@y;rho=(rx@ry)/(np.linalg.norm(rx)*np.linalg.norm(ry))
 key=tuple(labels)
 if key not in cache:
  blocks=[np.flatnonzero(labels.to_numpy()==v) for v in sorted(labels.unique())];perm=[]
  for parts in itertools.product(*[list(itertools.permutations(b)) for b in blocks]):
   ix=np.arange(len(d))
   for bb,pp in zip(blocks,parts):ix[bb]=pp
   perm.append(ix)
  cache[key]=np.array(perm)
 vals=y[cache[key]]@Q;den=np.linalg.norm(vals,axis=1)*np.linalg.norm(rx);null=vals[den>1e-12]@rx/den[den>1e-12];p=float(np.mean(abs(null)>=abs(rho)-1e-12))
 records.append(dict(gene=r.gene,scope=r.scope,variant=r.variant,model=r.model,program=r.program,rho_error=abs(rho-r.rho),P_error=abs(p-r.P),permutations=len(null)))
num=pd.DataFrame(records);check('independent primary correlations and exact P',num.rho_error.max()<1e-10 and num.P_error.max()<1e-12,len(num))
num.to_csv(O/'tables/independent_numerical_QA.tsv',sep='\t',index=False)
cand=set(pd.read_csv(O/'tables/frozen_candidate_pairs.tsv',sep='\t').gene);members=pd.read_csv(O/'tables/module_members.tsv',sep='\t')
check('candidate exclusion from scores',not bool(cand&set(members.loc[members.retained,'gene'])))
tf=pd.read_csv(O/'tables/TF_regulon_members.tsv',sep='\t');keep=tf[tf.retained]
check('TF scores exclude candidates and all program members',not bool(set(keep.gene)&(cand|set(members.gene))))
check('TF self excluded',not (keep.tf==keep.gene).any())
pairs=pd.read_csv(O/'tables/pairwise_disjoint_module_members.tsv',sep='\t')
for pair,z in pairs.groupby('pair'):
 sets=[set(v.gene) for _,v in z.groupby('module')];check('disjoint '+pair,not bool(sets[0]&sets[1]))
cover=pd.read_csv(O/'tables/sender_coverage.tsv',sep='\t')
check('all sender samples present',all(len(z)==10 and z.cells.min()>=10 for _,z in cover.groupby('source')))
for r in json.loads((O/'inputs/manifest.json').read_text()):
 check('original input hash '+Path(r['path']).name,hashlib.sha256(Path(r['path']).read_bytes()).hexdigest()==r['sha256'])
pd.DataFrame(checks).to_csv(O/'tables/QA_checks.tsv',sep='\t',index=False)
(O/'logs/numerical_QA_complete.json').write_text(json.dumps({'status':'success','checks':len(checks),'independent_correlations':len(num)},indent=2),encoding='utf-8')
print('QA_PASS',len(checks),len(num))
