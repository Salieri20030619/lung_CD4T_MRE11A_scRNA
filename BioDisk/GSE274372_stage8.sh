#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA
RUN=$(mktemp -d "$ROOT/02_analysis/GSE274372/stage8_$(date +%Y%m%d_%H%M%S)_XXXXXX")
export ROOT RUN OPENBLAS_NUM_THREADS=1
trap 'rc=$?; echo "GSE274372_STAGE8_FAIL EXIT_CODE=$rc"; exit "$rc"' ERR
cat > "$RUN/analysis.py" <<'PY'
import os
from pathlib import Path
import numpy as np,pandas as pd,scanpy as sc
root=Path(os.environ['ROOT']);run=Path(os.environ['RUN'])
for d in ['tables','objects']:(run/d).mkdir()
out=run/'tables'
meta=pd.read_csv(root/'02_analysis/GSE274372/stage7_20260908_112650_k6qMpd/tables/author_metadata.tsv.gz',sep='\t')
meta['sample_key']=['D%02d_R%d'%(int(d),int(r)) for d,r in zip(meta.Timepoint,meta.Rep)]
meta['barcode_key']=meta['cell'].astype(str).str.strip()
assert not meta.duplicated(['sample_key','barcode_key']).any()
a=sc.read_h5ad(root/'02_analysis/GSE274372/stage3_20260907_171359_MgRRUd/objects/whole_annotated.h5ad')
obs=a.obs.copy();obs['local_cell_id']=a.obs_names
obs['sample_key']=obs['sample'].astype(str);obs['barcode_key']=obs['barcode'].astype(str).str.strip()
assert not obs.duplicated(['sample_key','barcode_key']).any()
cols=['sample_key','barcode_key','author_cell_id','fine_annotation','broad_annotation','n.umi','num_genes_expressed']
j=obs.merge(meta[cols],on=['sample_key','barcode_key'],how='left',validate='one_to_one').set_index('local_cell_id').reindex(a.obs_names)
matched=j.author_cell_id.notna()
j['umi_equal']=np.where(matched,j.total_counts.eq(j['n.umi']),False)
j.to_csv(out/'whole_author_join.tsv.gz',sep='\t')
miss=meta.merge(obs[['sample_key','barcode_key']],on=['sample_key','barcode_key'],how='left',indicator=True)
miss.loc[miss['_merge'].eq('left_only')].to_csv(out/'author_unmatched.tsv.gz',sep='\t',index=False)
summary={'author_cells':len(meta),'local_cells':a.n_obs,'matched_cells':int(matched.sum()),'umi_mismatch':int((matched&~j.umi_equal).sum())}
pd.Series(summary).to_csv(out/'join_summary.tsv',sep='\t',header=False)
if summary['umi_mismatch']:raise RuntimeError('Matched UMI mismatch: inspect whole_author_join.tsv.gz')
for col in ['Broad','C06']:
 if col in j:pd.crosstab(j[col],j.fine_annotation).to_csv(out/f'{col}_author_cross.tsv',sep='\t')
keep=j.fine_annotation.isin(['Cd4+ T cells','Regulatory T cells'])
b=a[keep.to_numpy()].copy()
b.obs['author_identity']=j.loc[keep,'fine_annotation'].astype(str).to_numpy()
b.obs['author_cell_id']=j.loc[keep,'author_cell_id'].astype(str).to_numpy()
b.write_h5ad(run/'objects/CD4_author_candidates.h5ad',compression='gzip')
sy=b.var['symbol'].astype('string').fillna('').to_numpy()
markers=['Cd4','Cd8a','Cd8b1','Zbtb7b','Runx3','Cd3d','Cd3e','Cd3g','Bcl11b','Themis','Skap1','Tcf7','Lef1','Ccr7','Sell','Il7r','Foxp3','Il2ra','Ctla4','Ikzf2','Icos','Nkg7','Ccl5','Prf1','Sftpc','Scgb1a1','Lyz2','Cd79a','Mre11a']
e=b.obs.copy()
for g in markers:
 ix=np.flatnonzero(sy==g)
 e[g]=np.asarray(b.layers['counts'][:,ix].sum(axis=1)).ravel() if len(ix) else np.nan
assert (sy=='Mre11a').sum()==1
e.to_csv(out/'candidate_marker_counts.tsv.gz',sep='\t')
pd.DataFrame({'gene':markers,'feature_matches':[int((sy==g).sum()) for g in markers]}).to_csv(out/'marker_availability.tsv',sep='\t',index=False)
rows=[]
for sample in sorted(meta.sample_key.unique()):
 for scope in ['All','Cd4+ T cells','Regulatory T cells']:
  z=e[e['sample'].astype(str).eq(sample)]
  if scope!='All':z=z[z.author_identity.eq(scope)]
  n=len(z);reads=float(z.total_counts.sum());positive=int(z.Mre11a.gt(0).sum())
  rows.append([sample,scope,n,positive,float(z.Mre11a.sum()),positive/n*100 if n else np.nan,1e6*z.Mre11a.sum()/reads if reads else np.nan,int(z.dbl_020.sum()) if 'dbl_020' in z else np.nan])
pd.DataFrame(rows,columns=['sample','scope','cells','Mre_positive','Mre_UMI','Mre_pct','Mre_CPM','doublet_flags']).to_csv(out/'sample_summary.tsv',sep='\t',index=False)
# Reference labels define candidates for review; they do not replace independent identity evidence.
(out/'notes.txt').write_text('Sample + barcode exact join with raw UMI concordance. Author CD4/Treg candidates include original QC flags without automatic deletion. Whole-author unmatched cells are listed separately. No DEG or time trend test performed.\n')
print(summary,flush=True)
print('AUTHOR_CD4_RECOVERED='+str(b.n_obs),flush=True)
print('GSE274372_STAGE8_ANALYSIS_PASS',flush=True)
PY
bash -n "$0"
conda run -n batch_bench python -m py_compile "$RUN/analysis.py"
echo SYNTAX_CHECK=PASS
echo "OUTPUT=$RUN"
conda run --no-capture-output -n batch_bench python -u "$RUN/analysis.py"
tar -czf "$ROOT/tmp/GSE274372_stage8_return.tar.gz" -C "$RUN" tables analysis.py
echo GSE274372_STAGE8_PASS
