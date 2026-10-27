#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA
RUN=$(mktemp -d "$ROOT/02_analysis/GSE274372/stage6_$(date +%Y%m%d_%H%M%S)_XXXXXX")
export ROOT RUN OPENBLAS_NUM_THREADS=1
trap 'rc=$?; echo "GSE274372_STAGE6_FAIL EXIT_CODE=$rc"; exit "$rc"' ERR
cat > "$RUN/analysis.py" <<'PY'
import os
from pathlib import Path
import scanpy as sc,pandas as pd,numpy as np
root=Path(os.environ['ROOT']);run=Path(os.environ['RUN'])
for p in ['tables','objects']:(run/p).mkdir()
a=sc.read_h5ad(root/'02_analysis/GSE274372/stage5_20260907_224503_cUAf52/objects/T_refined.h5ad')
labels={'0':'NaiveUnk','1':'TregLike','2':'ActCD4','3':'CytoT','4':'CD8Like'}
a.obs['identity']=a.obs.K20R08.astype(str).map(labels).astype('category')
a.obs['CD4_candidate']=a.obs.K20R08.astype(str).isin(['1','2'])
a.obs['CD4_consensus']=a.obs.CD4_candidate&a.obs.K10R08.astype(str).isin(['1','4','6','7'])
a.obs.to_csv(run/'tables/T_identity.tsv.gz',sep='\t')
b=a[a.obs.CD4_candidate].copy();assert b.n_obs==261
b.write_h5ad(run/'objects/CD4_ActReg_candidates.h5ad',compression='gzip')
a[a.obs.identity=='NaiveUnk'].copy().write_h5ad(run/'objects/Naive_unresolved.h5ad',compression='gzip')
samples=sorted(a.obs['sample'].unique());meta=[]
idx=np.flatnonzero(a.var.symbol.eq('Mre11a').to_numpy());assert len(idx)==1
for scope in ['ActReg','Consensus','TregLike','ActCD4']:
 columns=[]
 for sample in samples:
  selected=a.obs.CD4_candidate.copy()
  if scope=='Consensus':selected=a.obs.CD4_consensus.copy()
  elif scope in ['TregLike','ActCD4']:selected=a.obs.identity.eq(scope)
  mask=(selected&a.obs['sample'].eq(sample)).to_numpy();n=int(mask.sum())
  values=np.asarray(a.layers['counts'][mask].sum(axis=0)).ravel();columns.append(values)
  target=np.asarray(a.layers['counts'][mask][:,idx].sum(axis=1)).ravel()
  day=int(a.obs.loc[a.obs['sample'].eq(sample),'day'].iloc[0])
  meta.append([scope,sample,day,n,int((target>0).sum()),int(target.sum()),float(values.sum()),100*(target>0).mean() if n else np.nan,float(target.sum()/values.sum()*1e6) if values.sum()>0 else np.nan])
 pd.DataFrame(np.asarray(columns).T,index=a.var_names,columns=samples).to_csv(run/f'tables/{scope}_counts.tsv.gz',sep='\t')
pd.DataFrame(meta,columns=['scope','sample','day','cells','Mre_detected','Mre_UMI','library_UMI','Mre_pct','Mre_CPM']).to_csv(run/'tables/sample_summary.tsv',sep='\t',index=False)
a.var[['symbol']].to_csv(run/'tables/features.tsv.gz',sep='\t')
(run/'tables/notes.txt').write_text('K20R08 clusters1/2: activated/regulatory CD4 candidates; not total CD4 population. Cluster0 remains naive T unresolved, not excluded as proven CD8. Consensus uses intersection with K10R08 clusters1/4/6/7. Cluster3 cytotoxic T, cluster4 CD8-like. Pseudobulk zero columns for empty subsets are placeholders, not observed zero-expression samples; exclude from expression testing. No DEG or time trend test run. Sample-level detection and CPM are descriptive. Day0 candidate cells0/2/2; D28 3/3/94.\n')
print('GSE274372_STAGE6_ANALYSIS_PASS')
PY
bash -n "$0"
conda run -n batch_bench python -m py_compile "$RUN/analysis.py"
echo 'SYNTAX_CHECK=PASS'
echo "OUTPUT=$RUN"
conda run --no-capture-output -n batch_bench python -u "$RUN/analysis.py"
tar -czf "$ROOT/tmp/GSE274372_stage6_return.tar.gz" -C "$RUN" tables analysis.py
echo 'GSE274372_STAGE6_PASS'
