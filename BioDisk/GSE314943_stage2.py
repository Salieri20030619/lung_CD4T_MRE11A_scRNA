import os, sys, gzip, importlib, json
from pathlib import Path
required = ['numpy', 'pandas', 'scipy', 'h5py', 'scanpy', 'anndata',
            'skimage', 'igraph', 'leidenalg', 'matplotlib']
missing = []
for name in required:
    try:
        importlib.import_module(name)
    except ImportError:
        missing.append(name)
if missing:
    raise SystemExit('MISSING_PACKAGES=' + ','.join(missing))
import numpy as np
import pandas as pd
import h5py
import scanpy as sc
import anndata as ad
from scipy import sparse, io
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
if len(sys.argv) == 2 and sys.argv[1] == '--check':
    assert hasattr(sc.pp, 'scrublet'), 'SCANPY_SCRUBLET_UNAVAILABLE'
    print('DEPENDENCY_CHECK=PASS')
    sys.exit(0)
out = Path(sys.argv[1])
raw = Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/00_raw/GSE314943/scRNA')
for folder in ['objects', 'tables', 'figures']:
    (out / folder).mkdir(exist_ok=True)
sc.settings.n_jobs = 8
sc.settings.verbosity = 2
samples = [('GSM9473144','Saline_1'), ('GSM9529135','Saline_2'),
           ('GSM9529136','Saline_3'), ('GSM9529140','Silica_1'),
           ('GSM9529141','Silica_2'), ('GSM9529142','Silica_3')]
objects, qc, mappings = [], [], []
for gsm, sample in samples:
    print('SAMPLE_START=' + sample, flush=True)
    p = str(raw / f'{gsm}_{sample}')
    if sample == 'Saline_2':
        x = io.mmread(p + '_matrix.mtx.gz').tocsc()
        f = pd.read_csv(p + '_features.tsv.gz', sep='\t', header=None)
        ids, symbols = f[0].astype(str).to_numpy(), f[1].astype(str).to_numpy()
        types = f[2].to_numpy() if f.shape[1] > 2 else np.repeat('Gene Expression', len(f))
        b = pd.read_csv(p + '_barcodes.tsv.gz', header=None)[0].astype(str).to_numpy()
    else:
        with h5py.File(p + '_filtered_feature_bc_matrix.h5') as h:
            m = h['matrix']
            x = sparse.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]), shape=tuple(m['shape'][:]))
            ids = m['features/id'].asstr()[:]
            symbols = m['features/name'].asstr()[:]
            types = m['features/feature_type'].asstr()[:]
            b = m['barcodes'].asstr()[:]
    keep = types == 'Gene Expression'
    x, ids, symbols = x[keep].T.tocsr(), ids[keep], symbols[keep]
    x.sum_duplicates(); x.eliminate_zeros()
    assert len(set(ids)) == len(ids) and len(set(b)) == len(b)
    assert np.isfinite(x.data).all() and (x.data >= 0).all()
    assert np.equal(x.data, np.floor(x.data)).all()
    a = ad.AnnData(x, obs=pd.DataFrame(index=[sample + ':' + v for v in b]),
                   var=pd.DataFrame({'symbol':symbols}, index=ids))
    a.obs['sample'], a.obs['group'], a.obs['gsm'] = sample, sample.split('_')[0], gsm
    a.var['mt'] = a.var.symbol.str.lower().str.startswith('mt-')
    assert a.var.mt.any()
    sc.pp.calculate_qc_metrics(a, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    passed = (a.obs.n_genes_by_counts >= 200) & (a.obs.total_counts >= 500) & (a.obs.pct_counts_mt <= 15)
    a.obs['basic_qc_pass'] = passed
    a.obs.to_csv(out / 'tables' / f'{sample}_basic_qc.tsv.gz', sep='\t')
    original = a.n_obs
    mappings.append(pd.DataFrame({'sample':sample,'gene_id':ids,'symbol':symbols}))
    a = a[passed].copy()
    sc.pp.scrublet(a, expected_doublet_rate=0.06, random_state=42, verbose=True)
    if 'predicted_doublet' not in a.obs or a.obs.predicted_doublet.isna().any():
        raise RuntimeError('SCRUBLET_THRESHOLD_FAILED=' + sample)
    a.obs.to_csv(out / 'tables' / f'{sample}_doublets.tsv.gz', sep='\t')
    a.write_h5ad(out / 'objects' / f'{sample}_qc_doublets.h5ad', compression='gzip')
    nd = int(a.obs.predicted_doublet.sum())
    qc.append(dict(sample=sample, input_cells=original, basic_qc=a.n_obs,
                   doublets=nd, singlets=a.n_obs-nd,
                   threshold=float(a.uns['scrublet']['threshold'])))
    objects.append(a[~a.obs.predicted_doublet].copy())
    print(f'SAMPLE_DONE={sample};SINGLETS={a.n_obs-nd}', flush=True)
pd.DataFrame(qc).to_csv(out/'tables/qc_summary.tsv',sep='\t',index=False)
mapping = pd.concat(mappings, ignore_index=True)
mapping.to_csv(out/'tables/gene_mapping.tsv.gz',sep='\t',index=False)
common = set(objects[0].var_names)
for a in objects[1:]:
    common &= set(a.var_names)
conflicts = mapping.groupby('gene_id').symbol.nunique()
conflicts = set(conflicts[conflicts > 1].index)
ids = [v for v in objects[0].var_names if v in common and v not in conflicts]
pd.DataFrame([{'sample':s,'input_features':a.n_vars,
               'shared_consistent':len(ids),'excluded':a.n_vars-len(ids)}
              for (_,s),a in zip(samples,objects)]).to_csv(out/'tables/feature_alignment.tsv',sep='\t',index=False)
assert 'ENSMUSG00000031928' in ids
v = objects[0].var.loc[ids].copy()
a = ad.concat([obj[:,ids] for obj in objects], join='inner', merge='same')
a.var = v
a.layers['counts'] = a.X.copy()
sc.pp.normalize_total(a, target_sum=1e4)
sc.pp.log1p(a)
sc.pp.highly_variable_genes(a, flavor='seurat', n_top_genes=3000, batch_key='sample')
h = a[:,a.var.highly_variable].copy()
h.layers.clear()
sc.pp.scale(h, max_value=10)
sc.tl.pca(h, n_comps=40, svd_solver='arpack', random_state=42)
a.obsm['X_pca'] = h.obsm['X_pca'].copy()
del h, objects
sc.pp.neighbors(a, n_neighbors=20, n_pcs=30, random_state=42)
for res,key in [(0.4,'C04'),(0.6,'C06'),(0.8,'C08')]:
    sc.tl.leiden(a, resolution=res, key_added=key, random_state=42, flavor='leidenalg')
sc.tl.umap(a, random_state=42)
print('CLUSTER_DONE;CELLS=' + str(a.n_obs), flush=True)
sc.tl.rank_genes_groups(a, groupby='C06', method='wilcoxon',
                       use_raw=False, n_genes=100, pts=True)
markers = sc.get.rank_genes_groups_df(a, group=None)
markers['symbol'] = markers.names.map(a.var.symbol)
markers.to_csv(out/'tables/cluster_top100.tsv',sep='\t',index=False)
pd.crosstab(a.obs.C06,a.obs['sample']).to_csv(out/'tables/cluster_sample.tsv',sep='\t')
meta = a.obs.copy()
meta[['UMAP1','UMAP2']] = a.obsm['X_umap']
meta.to_csv(out/'tables/cell_metadata.tsv.gz',sep='\t')
panel = ['Ptprc','Cd3d','Cd3e','Trbc2','Cd4','Il7r','Cd8a','Nkg7',
         'Ms4a1','Cd79a','Jchain','Lyz2','Adgre1','C1qa','Siglecf',
         'Itgax','Flt3','S100a8','S100a9','Cpa3','Kit','Epcam',
         'Sftpc','Ager','Scgb1a1','Foxj1','Pecam1','Cdh5','Prox1',
         'Col1a1','Dcn','Pdgfra','Rgs5','Acta2','Mki67','Mre11a']
rows = []
for cl in a.obs.C06.cat.categories:
    sub = a[a.obs.C06 == cl]
    for gene in panel:
        positions = np.flatnonzero(a.var.symbol.to_numpy() == gene)
        if not len(positions):
            continue
        z = sub.X[:,positions]
        c = sub.layers['counts'][:,positions]
        rows.append(dict(cluster=cl,gene=gene,cells=sub.n_obs,
            mean_log=float(z.sum()/sub.n_obs),
            pct=100*float((np.asarray(c.sum(axis=1)).ravel()>0).mean())))
pd.DataFrame(rows).to_csv(out/'tables/lineage_evidence.tsv',sep='\t',index=False)
for field in ['C06','sample','group']:
    fig, ax = plt.subplots(figsize=(9,7))
    sc.pl.umap(a,color=field,ax=ax,show=False,size=5,frameon=False,title=field)
    fig.savefig(out/'figures'/f'UMAP_{field}.png',dpi=350,bbox_inches='tight')
    fig.savefig(out/'figures'/f'UMAP_{field}.pdf',bbox_inches='tight')
    plt.close(fig)
a.write_h5ad(out/'objects/GSE314943_whole_lung.h5ad',compression='gzip')
with open(out/'tables/methods.txt','w') as f:
    f.write('QC: genes>=200; UMI>=500; MT<=15%.\nScrublet per sample; expected rate=0.06; automatic threshold.\n')
    f.write('Ensembl intersection; symbol conflicts excluded. X=log1p(CP10K); counts layer=original UMI.\n')
    f.write('Uncorrected PCA; sample-aware HVG; no group regression. C04/C06/C08=Leiden resolutions 0.4/0.6/0.8.\n')
    f.write('UMAP plots: clusters, samples, treatment groups; not final cell-type annotation.\n')
    f.write('Cluster markers are annotation evidence, not mouse-level disease DEG.\n')
with open(out/'tables/versions.json','w') as f:
    json.dump({m:str(getattr(importlib.import_module(m),'__version__','unknown')) for m in required},f,indent=2)
print('RAW_DATA_MODIFIED=NO')
print('GSE314943_STAGE2_PASS', flush=True)
import tarfile
archive = raw.parents[2] / 'tmp/GSE314943_stage2_return.tar.gz'
archive.parent.mkdir(exist_ok=True)
with tarfile.open(str(archive) + '.part', 'w:gz') as t:
    t.add(out/'tables', arcname='tables')
    t.add(out/'figures', arcname='figures')
os.replace(str(archive) + '.part', archive)
print('RETURN_PACKAGE=' + str(archive), flush=True)
print('GSE314943_STAGE2_DELIVERY_PASS', flush=True)
