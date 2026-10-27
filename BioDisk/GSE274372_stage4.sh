#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA
RUN=$(mktemp -d "$ROOT/02_analysis/GSE274372/stage4_$(date +%Y%m%d_%H%M%S)_XXXXXX")
export ROOT RUN OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4 MPLBACKEND=Agg
trap 'rc=$?; echo "GSE274372_STAGE4_FAIL EXIT_CODE=$rc"; exit "$rc"' ERR
cat > "$RUN/analysis.py" <<'PY'
import os
from pathlib import Path
import numpy as np,pandas as pd,scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(os.environ['ROOT']);run=Path(os.environ['RUN'])
for sub in ['tables','figures','objects']:(run/sub).mkdir()
a=sc.read_h5ad(root/'02_analysis/GSE274372/stage3_20260907_171359_MgRRUd/objects/TNK_reclustered.h5ad')
a.obs['origin_T05']=a.obs.T05.astype(str)
b=a[a.obs.origin_T05.isin(['1','2','3','5'])].copy()
assert b.n_obs>500
b.obs['candidate_status']=b.obs.origin_T05.map({'1':'MixedT','2':'ActRegT','3':'LowRNA','5':'UnconvT'}).astype('category')
# Keep flagged cells to inspect their cluster distribution before exclusion.
b.X=b.layers['counts'].copy();b.uns.clear()
for key in list(b.obsp.keys()):del b.obsp[key]
for key in list(b.obsm.keys()):del b.obsm[key]
sc.pp.normalize_total(b,target_sum=1e4);sc.pp.log1p(b)
sc.pp.highly_variable_genes(b,n_top_genes=2000,flavor='seurat',batch_key='sample')
v=b[:,b.var.highly_variable].copy();sc.pp.scale(v,max_value=10)
sc.tl.pca(v,n_comps=30,svd_solver='arpack',random_state=42)
b.obsm['X_pca']=v.obsm['X_pca'].copy()
sc.pp.neighbors(b,n_neighbors=15,n_pcs=20,random_state=42);sc.tl.umap(b,random_state=42)
for res,key in [(0.3,'T03'),(0.5,'T05'),(0.8,'T08')]:sc.tl.leiden(b,resolution=res,key_added=key,random_state=42)
sc.tl.rank_genes_groups(b,'T05',method='wilcoxon',use_raw=False,n_genes=100)
top=sc.get.rank_genes_groups_df(b,group=None);top['symbol']=top.names.map(b.var.symbol)
top.to_csv(run/'tables/TNK_top100.tsv',sep='\t',index=False)
panel=['Cd3d','Cd3e','Cd3g','Trac','Trbc1','Trbc2','Cd247','Skap1','Themis','Bcl11b','Cd4','Cd8a','Cd8b1','Tcf7','Lef1','Ccr7','Sell','Il7r','Foxp3','Il2ra','Ctla4','Ikzf2','Tnfrsf4','Icos','Pdcd1','Gata3','Rorc','Il17a','Tbx21','Ifng','Nkg7','Klrd1','Klrk1','Eomes','Prf1','Gzma','Gzmb','Ccl5','Trdc','Trgc1','Trgc2','Zbtb16','Mki67','Top2a','Cd79a','Ms4a1','Csf1r','Lyz2','Sftpc','Scgb1a1','Pecam1','Mre11a']
rows=[]
for cl in b.obs.T05.cat.categories:
 mask=b.obs.T05.eq(cl).to_numpy();c=b.layers['counts'][mask]
 for sym in panel:
  gi=np.flatnonzero(b.var.symbol.eq(sym).to_numpy())
  if not len(gi):continue
  vals=np.asarray(c[:,gi].sum(axis=1)).ravel();expr=np.asarray(b.X[mask][:,gi].sum(axis=1)).ravel()
  rows.append([cl,sym,int(mask.sum()),int((vals>0).sum()),100*(vals>0).mean(),expr.mean()])
pd.DataFrame(rows,columns=['cluster','symbol','cells','detected','pct','mean_log']).to_csv(run/'tables/TNK_panel.tsv',sep='\t',index=False)
pd.crosstab(b.obs.T05,b.obs['sample']).to_csv(run/'tables/TNK_sample.tsv',sep='\t')
pd.crosstab(b.obs.T05,b.obs.origin_T05).to_csv(run/'tables/TNK_origin.tsv',sep='\t')
b.obs.groupby('T05',observed=True)[['total_counts','n_genes','pct_mt','dbl_score','dbl_020']].mean().to_csv(run/'tables/TNK_QC.tsv',sep='\t')
b.obs.to_csv(run/'tables/TNK_cells.tsv.gz',sep='\t')
for col in ['T05','sample','Day','candidate_status']:
 sc.pl.umap(b,color=col,show=False,frameon=False,title=col)
 plt.savefig(run/f'figures/TNK_{col}.png',dpi=350,bbox_inches='tight');plt.close('all')
# Export count-level evidence; missing features are NA, never measured zeros.
cells=b.obs.copy();audit=[]
for sym in panel:
 gi=np.flatnonzero(b.var.symbol.eq(sym).to_numpy())
 audit.append([sym,len(gi)])
 cells[sym]=np.asarray(b.layers['counts'][:,gi].sum(axis=1)).ravel() if len(gi) else np.nan
cells.to_csv(run/'tables/T_cell_marker_counts.tsv.gz',sep='\t')
pd.DataFrame(audit,columns=['symbol','matched_features']).to_csv(run/'tables/marker_availability.tsv',sep='\t',index=False)
b.uns['dataset']='GSE274372';b.write_h5ad(run/'objects/T_candidates.h5ad',compression='gzip')
(run/'tables/captions.txt').write_text('T05: new T-cell candidate clusters, Leiden resolution0.5. origin_T05: Stage3 clusters. MixedT: Stage3 cluster1, mixed CD4/CD8 evidence; ActRegT: cluster2, activated/regulatory T candidate; LowRNA: cluster3, low RNA unresolved candidate; UnconvT: cluster5, unconventional T candidate. These are working labels, not final CD4 subtypes. Sample: day and submitted replicate label; Day: D00/D07/D28/D56. No cell deletion by doublet score yet. Trac/Trbc/Trdc/Trgc absent from supplied feature list: unavailable measurements, not negative lineage evidence. Raw counts retained.\n')
print(f'TNK_CELLS={b.n_obs};TNK_CLUSTERS={b.obs.T05.nunique()}')
print('GSE274372_STAGE4_ANALYSIS_PASS')
PY
bash -n "$0"
conda run -n batch_bench python -m py_compile "$RUN/analysis.py"
echo 'SYNTAX_CHECK=PASS'
echo "OUTPUT=$RUN"
conda run --no-capture-output -n batch_bench python -u "$RUN/analysis.py"
tar -czf "$ROOT/tmp/GSE274372_stage4_return.tar.gz" -C "$RUN" tables figures analysis.py
echo 'GSE274372_STAGE4_PASS'
