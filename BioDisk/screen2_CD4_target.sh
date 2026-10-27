#!/usr/bin/env bash
set -Eeuo pipefail
export PROJECT_ROOT=/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA
mkdir -p "$PROJECT_ROOT/tmp" "$PROJECT_ROOT/02_analysis" "$PROJECT_ROOT/logs"
exec 9>"$PROJECT_ROOT/tmp/screen2_CD4_target.lock"
if ! flock -n 9; then echo TARGET_ALREADY_RUNNING; exit 0; fi
trap 'rc=$?; echo "TARGET_FAIL EXIT_CODE=$rc"; exit "$rc"' ERR
TASK=$(mktemp -d "$PROJECT_ROOT/02_analysis/screen2_CD4_target_XXXXXX")
cat > "$TASK/target.py" <<'TARGET_PY'
from pathlib import Path
import os,sys,json,itertools,tarfile
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
R=Path(os.environ['PROJECT_ROOT']);O=Path(sys.argv[1]);S=R/'02_analysis/screen2_Trefine_nuzVQu'
config={'GSE174725':{'core':{'0':'MemCD4','1':'NaiveCD4','12':'Treg'},'extra':['2','3'],'samples':{'GSM5324866':('Exposed',0),'GSM5324867':('Exposed',0),'GSM5324868':('Silicosis',0),'GSM5324869':('Silicosis',0),'GSM5324870':('Silicosis',0)},'n':10286},'GSE183682':{'core':{'0':'NaiveCD4','3':'ActCD4','5':'Treg'},'extra':[],'samples':{'GSM5567501':('NS',7),'GSM5567502':('Silica',7),'GSM5567503':('NS',56),'GSM5567504':('Silica',56)},'n':3266}}
summary=[];tests=[];contrasts=[]
for g,c in config.items():
 d=O/g;d.mkdir(parents=True);a=sc.read_h5ad(S/g/'T_refined.h5ad');assert a.n_obs==c['n'];labels=a.obs.F08.astype(str)
 assert set(c['core'])<=set(labels);target='MRE11' if g=='GSE174725' else 'MRE11A';assert target in a.var_names
 a.obs['core_label']=labels.map(c['core']).fillna('Unresolved');a.obs['scope_core']=labels.isin(c['core']);a.obs['scope_expanded']=labels.isin(list(c['core'])+c['extra'])
 a.obs.to_csv(d/'scope_membership.tsv.gz',sep='\t')
 scopes=['Core','Expanded'] if c['extra'] else ['Core']
 for scope in scopes:
  b=a[np.asarray(a.obs['scope_'+scope.lower()])].copy();X=sparse.csr_matrix(b.layers['counts']);assert np.all(X.data>=0)
  v=np.asarray(X[:,b.var_names==target].sum(1)).ravel();lib=np.asarray(X.sum(1)).ravel();pb={}
  for sample,(group,day) in c['samples'].items():
   z=np.asarray(b.obs['sample']==sample);n=int(z.sum());vec=np.asarray(X[z].sum(0)).ravel();pb[sample]=vec
   total=vec.sum();count=v[z].sum();cpm=count/total*1e6 if total else np.nan
   summary.append(dict(dataset=g,scope=scope,sample=sample,group=group,day=day,n_cells=n,target_positive=int((v[z]>0).sum()),target_pct=100*np.mean(v[z]>0) if n else np.nan,target_counts=count,total_counts=total,target_CPM=cpm,log2_CPM1=np.log2(1+cpm),median_library=np.median(lib[z]) if n else np.nan,eligible20=n>=20))
  pd.DataFrame(pb,index=b.var_names).to_csv(d/(scope+'_pseudobulk.tsv.gz'),sep='\t')
  pd.crosstab(b.obs['sample'],b.obs.core_label).reindex(index=list(c['samples']),fill_value=0).to_csv(d/(scope+'_composition.tsv'),sep='\t')
  if scope=='Core':b.write_h5ad(d/'CD4_core.h5ad',compression='gzip')
 t=pd.DataFrame([x for x in summary if x['dataset']==g])
 if g=='GSE174725':
  for scope in scopes:
   # Samples have equal weight; no cells treated as independent subjects.
   for threshold in [20,1]:
    q=t[(t.scope==scope)&(t.n_cells>=threshold)].copy();case=q.group.eq('Silicosis').values;ncase=int(case.sum());ncontrol=int((~case).sum())
    row=dict(dataset=g,scope=scope,min_cells=threshold,n_case=ncase,n_control=ncontrol,metric='log2_CPM1',method='exact_sample_label_permutation')
    if min(ncase,ncontrol)>=2:
     vals=q.log2_CPM1.to_numpy();effect=float(vals[case].mean()-vals[~case].mean());null=[]
     for ix in itertools.combinations(range(len(vals)),ncase):
      z=np.zeros(len(vals),dtype=bool);z[list(ix)]=True;null.append(float(vals[z].mean()-vals[~z].mean()))
     row.update(effect=effect,p_nominal=float(np.mean(np.abs(null)>=abs(effect)-1e-12)),permutations=len(null),status='exploratory')
    else:row.update(status='insufficient_samples')
    tests.append(row)
 else:
  for day in [7,56]:
   q=t[(t.scope=='Core')&(t.day==day)].set_index('group');x=q.loc['NS'];y=q.loc['Silica']
   contrasts.append(dict(dataset=g,day=day,NS_CPM=x.target_CPM,Silica_CPM=y.target_CPM,log2_CPM1_difference=y.log2_CPM1-x.log2_CPM1,NS_pct=x.target_pct,Silica_pct=y.target_pct,p_status='not_tested_one_library_per_condition'))
 (d/'annotation.json').write_text(json.dumps(c,indent=2))
 print('TARGET_SCOPE_PASS',g,flush=True)
pd.DataFrame(summary).to_csv(O/'sample_target.tsv',sep='\t',index=False)
pd.DataFrame(tests).to_csv(O/'human_nominal_tests.tsv',sep='\t',index=False)
pd.DataFrame(contrasts).to_csv(O/'mouse_descriptive_contrasts.tsv',sep='\t',index=False)
(O/'README.txt').write_text('Core definitions fixed from F08 lineage markers, not target expression. Human Expanded additionally includes mixed CD4/CD8 clusters 2/3; not a purified CD4 population. Human nominal p-values: exact permutation of sample labels on mean log2(CPM+1), equal sample weight, two-sided absolute statistic; observational group exchangeability is an assumption. Core >=20 cells is the primary exploratory comparison; all-sample and Expanded analyses are sensitivity descriptions, not opportunities to select a smaller p-value. No multiple-testing adjustment included. Mouse has one library per condition/timepoint; no disease p-value calculated. CPM denominator is all measured shared-feature counts in each selected pseudobulk. MemCD4=memory-like CD4; NaiveCD4=naive-like CD4; ActCD4=activated-like CD4; Treg=regulatory T. Core represents a conservative working subset, not proof that unresolved cells lack CD4 identity.\n')
with tarfile.open(R/'tmp/screen2_CD4_target_return.tar.gz','w:gz') as tar:
 for p in O.rglob('*'):
  if p.is_file() and p.suffix!='.h5ad' and 'pseudobulk' not in p.name:tar.add(p,arcname=str(p.relative_to(O)))
print('SCREEN2_CD4_TARGET_PASS',flush=True)
TARGET_PY
conda run -n batch_bench python -m py_compile "$TASK/target.py"
echo SYNTAX_CHECK=PASS
echo "OUTPUT=$TASK"
conda run --no-capture-output -n batch_bench python -u "$TASK/target.py" "$TASK"
