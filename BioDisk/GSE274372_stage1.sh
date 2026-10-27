#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA
mkdir -p "$ROOT/02_analysis/GSE274372" "$ROOT/logs" "$ROOT/tmp"
RUN=$(mktemp -d "$ROOT/02_analysis/GSE274372/stage1_$(date +%Y%m%d_%H%M%S)_XXXXXX")
export ROOT RUN OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
trap 'rc=$?; echo "GSE274372_STAGE1_FAIL EXIT_CODE=$rc"; exit "$rc"' ERR
cat > "$RUN/analysis.py" <<'PY'
import os, re, gzip
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from scipy.io import mmread
root=Path(os.environ['ROOT']); run=Path(os.environ['RUN'])
raw=root/'00_raw/GSE274372'; out=run/'tables';out.mkdir()
(run/'objects').mkdir()
assert 'GSE274372_DOWNLOAD_PASS' in (raw/'download_audit.txt').read_text()
paths=sorted((raw/'extracted').glob('*.umi_counts.mtx.gz'))
assert len(paths)==12
objects=[]; summaries=[]; markers=[]; designs=[]
reference=None
panel=['Mre11a','Mre11','Cd3d','Cd3e','Trac','Cd4','Cd8a','Cd8b1','Tcf7','Lef1','Il7r','Foxp3','Il2ra','Nkg7','Csf1r','Spp1']
for path in paths:
    stem=path.name.removesuffix('.umi_counts.mtx.gz')
    match=re.fullmatch(r'(GSM\d+)_\d+\.silicosis\.day(\d+)\.(\d+)',stem)
    assert match,stem
    gsm,day,rep=match.groups();day=int(day)
    cells=pd.read_csv(path.with_name(stem+'.cell_annotations.txt.gz'),sep='\t',header=None,dtype=str,keep_default_na=False)
    genes=pd.read_csv(path.with_name(stem+'.gene_annotations.txt.gz'),sep='\t',header=None,dtype=str,keep_default_na=False)
    assert cells.shape[1]==1 and genes.shape[1]==2, 'Unexpected annotation columns'
    assert cells[0].is_unique and genes[0].is_unique
    if reference is None: reference=genes.copy()
    assert genes.equals(reference), 'Gene order differs across samples'
    with gzip.open(path,'rb') as f: matrix=sparse.csr_matrix(mmread(f))
    if matrix.shape==(len(genes),len(cells)):
        matrix=matrix.T.tocsr(); orientation='genes_by_cells'
    elif matrix.shape==(len(cells),len(genes)):
        orientation='cells_by_genes'
    else: raise ValueError((stem,matrix.shape,len(cells),len(genes)))
    assert np.isfinite(matrix.data).all() and (matrix.data>=0).all()
    assert np.equal(matrix.data,np.floor(matrix.data)).all(), 'Noninteger counts'
    matrix.sum_duplicates();matrix.eliminate_zeros()
    umi=np.asarray(matrix.sum(axis=1)).ravel(); ng=np.diff(matrix.indptr)
    mt=genes[1].str.lower().str.startswith('mt-').to_numpy()
    assert mt.any(), 'No mitochondrial symbols; inspect feature mapping'
    mtumi=np.asarray(matrix[:,mt].sum(axis=1)).ravel()
    pct=np.divide(mtumi*100,umi,out=np.zeros_like(umi,dtype=float),where=umi>0)
    sample=f'D{day:02d}_R{rep}'
    obs=pd.DataFrame({'barcode':cells[0].values,'sample':sample,'GSM':gsm,'day':day,'replicate_label':rep,'total_counts':umi,'n_genes':ng,'pct_mt':pct},index=gsm+':'+cells[0].astype(str))
    obs.index.name='cell_id'
    var=pd.DataFrame({'symbol':genes[1].values},index=genes[0].values);var.index.name='gene_id'
    objects.append(ad.AnnData(X=matrix,obs=obs,var=var))
    row={'sample':sample,'GSM':gsm,'day':day,'replicate_label':rep,'cells':len(cells),'genes':len(genes),'orientation':orientation,'zero_UMI_cells':int((umi==0).sum()),'mt_genes':int(mt.sum())}
    for name,values in [('UMI',umi),('genes',ng),('pct_mt',pct)]:
        for q in [0,.05,.25,.5,.75,.95,1]:row[f'{name}_q{int(q*100):02d}']=float(np.quantile(values,q))
    summaries.append(row)
    for symbol in panel:
        idx=np.flatnonzero(genes[1].eq(symbol).to_numpy())
        vals=np.asarray(matrix[:,idx].sum(axis=1)).ravel() if len(idx) else np.zeros(len(cells))
        markers.append({'sample':sample,'day':day,'symbol':symbol,'matched_features':len(idx),'detected_cells':int((vals>0).sum()),'pct_detected':100*float((vals>0).mean()),'total_UMI':float(vals.sum()),'scope':'whole_lung_before_new_QC'})
    designs.append({'sample':sample,'GSM':gsm,'day':day,'replicate_label':rep,'source_file':path.name})
    print(f'READ_PASS={sample};CELLS={len(cells)};MEDIAN_UMI={np.median(umi):.1f}',flush=True)
assert sorted(x['day'] for x in designs)==[0]*3+[7]*3+[28]*3+[56]*3
merged=ad.concat(objects,axis=0,join='inner',merge='same',index_unique=None)
assert merged.obs_names.is_unique and merged.n_obs==sum(x.n_obs for x in objects)
merged.uns['dataset']='GSE274372'
merged.uns['stage']='Input counts; no additional QC filtering or cell identity assignment'
merged.write_h5ad(run/'objects/whole_input.h5ad',compression='gzip')
pd.DataFrame(summaries).sort_values(['day','sample']).to_csv(out/'sample_QC_summary.tsv',sep='\t',index=False)
pd.DataFrame(markers).to_csv(out/'marker_detection.tsv',sep='\t',index=False)
pd.DataFrame(designs).sort_values(['day','sample']).to_csv(out/'sample_design.tsv',sep='\t',index=False)
merged.obs.to_csv(out/'cell_QC.tsv.gz',sep='\t')
merged.var.to_csv(out/'features.tsv.gz',sep='\t')
(out/'notes.txt').write_text('Published count matrices; original prefiltering must be considered. No new QC thresholds applied. Marker detection is whole-lung, not CD4-specific. Day and replicate labels follow supplied filenames; matching replicate numbers across days do not establish paired animals. The author CDS object is retained separately.\n')
print(f'TOTAL_CELLS={merged.n_obs};FEATURES={merged.n_vars}',flush=True)
print('GSE274372_STAGE1_ANALYSIS_PASS',flush=True)
PY
bash -n "$0"
conda run -n batch_bench python -m py_compile "$RUN/analysis.py"
echo 'SYNTAX_CHECK=PASS'
echo "OUTPUT=$RUN"
conda run --no-capture-output -n batch_bench python -u "$RUN/analysis.py"
tar -czf "$RUN/return.tar.gz" -C "$RUN" tables analysis.py
cp "$RUN/return.tar.gz" "$ROOT/tmp/GSE274372_stage1_return.tar.gz"
echo 'GSE274372_STAGE1_PASS'
