#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA
mkdir -p "$ROOT/02_analysis/GSE274372" "$ROOT/logs" "$ROOT/tmp"
RUN=$(mktemp -d "$ROOT/02_analysis/GSE274372/stage2_$(date +%Y%m%d_%H%M%S)_XXXXXX")
export ROOT RUN OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg
trap 'rc=$?; echo "GSE274372_STAGE2_FAIL EXIT_CODE=$rc"; exit "$rc"' ERR
cat > "$RUN/analysis.py" <<'PY'
import os
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(os.environ['ROOT']);run=Path(os.environ['RUN'])
for sub in ['tables','figures','objects']:(run/sub).mkdir()
sc.settings.n_jobs=4
src=root/'02_analysis/GSE274372/stage1_20260907_105231_SGOZ4K/objects/whole_input.h5ad'
a=sc.read_h5ad(src)
assert a.uns['dataset']=='GSE274372'
q=a.obs.copy();q['base_pass']=(q.total_counts>100)&(q.pct_mt<10)
q['genes100_pass']=q.base_pass&(q.n_genes>=100)
q['genes200_pass']=q.base_pass&(q.n_genes>=200)
q.groupby('sample',observed=True)[['base_pass','genes100_pass','genes200_pass']].sum().to_csv(run/'tables/QC_sensitivity.tsv',sep='\t')
q.to_csv(run/'tables/QC_all.tsv.gz',sep='\t')
a=a[q.base_pass.to_numpy()].copy()
assert a.n_obs==29831, 'Baseline count differs from reviewed Stage1'
a.obs['dbl_score']=np.nan;a.obs['dbl_auto']=False
for sample in sorted(a.obs['sample'].unique()):
    idx=a.obs['sample'].eq(sample);b=a[idx].copy()
    sc.pp.scrublet(b,expected_doublet_rate=0.06,random_state=42,verbose=False)
    a.obs.loc[idx,'dbl_score']=b.obs['doublet_score'].to_numpy()
    a.obs.loc[idx,'dbl_auto']=b.obs['predicted_doublet'].to_numpy()
    print(f'SCRUBLET_PASS={sample}',flush=True)
assert np.isfinite(a.obs.dbl_score).all()
a.obs['dbl_020']=a.obs.dbl_score>=0.2
# Scores are flags; no doublet-driven deletion until cluster review.
a.layers['counts']=a.X.copy()
a.obs['sample']=a.obs['sample'].astype(str).astype('category')
a.obs['Day']=a.obs.day.map(lambda v:f'D{int(v):02d}').astype('category')
sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a)
sc.pp.highly_variable_genes(a,n_top_genes=3000,flavor='seurat',batch_key='sample')
v=a[:,a.var.highly_variable].copy()
sc.pp.scale(v,max_value=10)
sc.tl.pca(v,n_comps=40,svd_solver='arpack',random_state=42)
a.obsm['X_pca']=v.obsm['X_pca'].copy()
sc.pp.neighbors(a,n_neighbors=20,n_pcs=30,random_state=42)
sc.tl.umap(a,random_state=42)
for res,key in [(0.4,'C04'),(0.6,'C06'),(0.8,'C08')]:
    sc.tl.leiden(a,resolution=res,key_added=key,random_state=42)
sc.tl.rank_genes_groups(a,'C06',method='wilcoxon',use_raw=False,n_genes=60)
deg=sc.get.rank_genes_groups_df(a,group=None)
deg['symbol']=deg.names.map(a.var.symbol)
deg.to_csv(run/'tables/cluster_top60.tsv',sep='\t',index=False)
panel=['Cd3d','Cd3e','Trac','Cd4','Cd8a','Cd8b1','Tcf7','Lef1','Il7r','Foxp3','Il2ra','Nkg7','Klrd1','Csf1r','Adgre1','Lyz2','Spp1','C1qa','Mertk','Siglecf','Pparg','Ly6c2','Flt3','Cd79a','Ms4a1','Jchain','S100a8','Csf3r','Epcam','Sftpc','Ager','Scgb1a1','Foxj1','Col1a1','Pdgfra','Pecam1','Cdh5','Prox1','Acta2','Rgs5','Mki67','Mre11a']
rows=[]
for cl in a.obs.C06.cat.categories:
    mask=a.obs.C06.eq(cl).to_numpy();counts=a.layers['counts'][mask]
    for sym in panel:
        gi=np.flatnonzero(a.var.symbol.eq(sym).to_numpy())
        if not len(gi):continue
        vals=np.asarray(counts[:,gi].sum(axis=1)).ravel()
        expr=np.asarray(a.X[mask][:,gi].sum(axis=1)).ravel()
        rows.append([cl,sym,int(mask.sum()),int((vals>0).sum()),100*(vals>0).mean(),expr.mean()])
pd.DataFrame(rows,columns=['cluster','symbol','cells','detected','pct','mean_log']).to_csv(run/'tables/cluster_panel.tsv',sep='\t',index=False)
pd.crosstab(a.obs.C06,a.obs['sample']).to_csv(run/'tables/cluster_sample.tsv',sep='\t')
a.obs.groupby('C06',observed=True)[['total_counts','n_genes','pct_mt','dbl_score','dbl_020','dbl_auto']].mean().to_csv(run/'tables/cluster_QC.tsv',sep='\t')
a.obs.to_csv(run/'tables/cells.tsv.gz',sep='\t')
for col in ['C06','Day','sample','dbl_score']:
    sc.pl.umap(a,color=col,show=False,frameon=False,title=col)
    plt.savefig(run/f'figures/UMAP_{col}.png',dpi=300,bbox_inches='tight');plt.close('all')
(run/'tables/figure_captions.txt').write_text('UMAP: uniform manifold approximation and projection. C06: provisional Leiden clusters at resolution 0.6. Day: D00/D07/D28/D56, sampled days 0/7/28/56. Sample: Dxx_Ry, day and submitted replicate label. dbl_score: per-sample Scrublet doublet score. All plots retain flagged doublet candidates pending marker review. PCA/UMAP not batch-corrected; sample and day distributions must be reviewed before integration.\n')
(run/'tables/methods.txt').write_text('Input UMI>100, mitochondrial percentage<10; no additional gene-count cutoff. Scrublet expected rate 0.06 is a working prior, not measured sci-RNA-seq collision rate. Report automatic and >=0.2 flags; no deletion yet. HVG: 3000, batch_key sample; scale capped10; PCA40, neighbors20 using30PC; Leiden0.4/0.6/0.8. Marker P values are cell-level descriptive annotation statistics. Raw counts preserved in counts layer.\n')
a.write_h5ad(run/'objects/whole_screened.h5ad',compression='gzip')
print(f'CELLS={a.n_obs};CLUSTERS={a.obs.C06.nunique()};DBL_FLAG={int(a.obs.dbl_020.sum())}')
print('GSE274372_STAGE2_ANALYSIS_PASS')
PY
bash -n "$0"
conda run -n batch_bench python -m py_compile "$RUN/analysis.py"
echo 'SYNTAX_CHECK=PASS'
echo "OUTPUT=$RUN"
conda run --no-capture-output -n batch_bench python -u "$RUN/analysis.py"
tar -czf "$ROOT/tmp/GSE274372_stage2_return.tar.gz" -C "$RUN" tables figures analysis.py
echo 'GSE274372_STAGE2_PASS'
