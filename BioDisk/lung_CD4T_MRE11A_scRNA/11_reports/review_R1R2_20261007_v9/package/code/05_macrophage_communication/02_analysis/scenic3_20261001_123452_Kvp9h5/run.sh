#!/usr/bin/env bash
set -Eeuo pipefail
ROOT="${ROOT:-/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA}"
PREV1="$ROOT/02_analysis/mechanism_step1_20260927_152218_s4V5eF"
PREV2="$ROOT/02_analysis/mechanism_step2_20260927_154919_EsN4gk"
PREV8="$ROOT/02_analysis/mechanism_step8_20260929_101407_jWFgLl"
SCENIC1="$ROOT/02_analysis/scenic1_20260929_103513_xRoEnr"
SCENIC2="$ROOT/02_analysis/scenic2_20261001_092655_ZKNWD6"
export ROOT PREV1 PREV2 PREV8 SCENIC1 SCENIC2
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=4
bash -n "$0"
mkdir -p "$ROOT/logs" "$ROOT/tmp" "$ROOT/02_analysis"
exec 9>"$ROOT/tmp/lung_scenic3.lock"
flock -n 9 || { echo SCENIC3_ALREADY_RUNNING; exit 1; }
OUT=$(mktemp -d "$ROOT/02_analysis/scenic3_$(date +%Y%m%d_%H%M%S)_XXXXXX")
export OUT
trap 'rc=$?; if ((rc != 0)); then printf "LUNG_SCENIC3_FAIL EXIT_CODE=%s OUTPUT=%s\n" "$rc" "$OUT"; fi' EXIT
mkdir -p "$OUT/input" "$OUT/tables" "$OUT/objects" "$OUT/figures"
cp -- "$0" "$OUT/run.sh"
printf 'OUTPUT=%s\n' "$OUT"
cat > "$OUT/export.R" <<'CODE_EXPORT'
args <- commandArgs(TRUE); p8 <- args[1]; out <- args[2]
suppressPackageStartupMessages(library(Matrix))
ck <- readRDS(file.path(p8,'objects/marker_audit_checkpoint.rds'))
md <- ck$metadata; col <- 'RNA_snn_res.0.4'
stopifnot(col%in%names(md),!anyDuplicated(md$cell))
sel <- as.character(md[[col]])%in%c('3','6','16','22')
md <- md[sel,,drop=FALSE]; wanted <- md$cell
stopifnot(length(wanted)>0,all(is.finite(ck$library)),file.exists(ck$source))
oldlib <- ck$library[match(wanted,ck$metadata$cell)]
cat('HUMAN_SELECTED_CELLS=',length(wanted),'\n',sep='');flush.console()
suppressPackageStartupMessages(library(Seurat))
a <- readRDS(ck$source)
stopifnot(inherits(a,'Seurat'),'RNA'%in%Assays(a),all(wanted%in%colnames(a)))
layers <- if(inherits(a[['RNA']],'Assay5'))Layers(a[['RNA']],search='^counts') else 'counts'
pieces <- list(); seen <- integer(length(wanted)); audit <- list()
for(ly in layers){
 x <- if(inherits(a[['RNA']],'Assay5'))LayerData(a,assay='RNA',layer=ly) else GetAssayData(a,assay='RNA',slot='counts')
 keep <- which(colnames(x)%in%wanted)
 if(!length(keep))next
 x <- as(x[,keep,drop=FALSE],'CsparseMatrix'); ix <- match(colnames(x),wanted)
 stopifnot(!anyNA(ix),all(seen[ix]==0),all(is.finite(x@x)),all(x@x>=0),all(abs(x@x-round(x@x))<1e-6))
 stopifnot(isTRUE(all.equal(as.numeric(Matrix::colSums(x)),oldlib[ix],tolerance=1e-8)))
 genes <- toupper(rownames(x));stopifnot(!anyNA(genes),all(nzchar(genes)))
 rownames(x) <- genes;pieces[[ly]] <- x;seen[ix] <- seen[ix]+1L
 audit[[ly]] <- data.frame(layer=ly,cells=ncol(x),features=nrow(x))
}
stopifnot(all(seen==1));rm(a);gc(FALSE)
genes <- sort(unique(unlist(lapply(pieces,rownames),use.names=FALSE)))
# Aggregate duplicate symbols and place every cell in its original checkpoint order.
ii <- jj <- vv <- list()
for(k in seq_along(pieces)){
 x <- pieces[[k]];s <- summary(x)
 ii[[k]] <- match(rownames(x)[s$i],genes);jj[[k]] <- match(colnames(x)[s$j],wanted);vv[[k]] <- s$x
}
x <- sparseMatrix(i=unlist(ii,use.names=FALSE),j=unlist(jj,use.names=FALSE),x=unlist(vv,use.names=FALSE),dims=c(length(genes),length(wanted)),dimnames=list(genes,wanted))
stopifnot(isTRUE(all.equal(as.numeric(Matrix::colSums(x)),oldlib,tolerance=1e-8)))
md$parent <- as.character(md[[col]]);md$library <- oldlib
wr <- function(d,p)write.table(d,p,sep='\t',quote=FALSE,row.names=FALSE,na='NA')
wr(md,file.path(out,'input/human_mac_cells.tsv'))
wr(data.frame(gene=genes),file.path(out,'input/human_mac_genes.tsv'))
writeMM(x,file.path(out,'input/human_mac_counts.mtx'))
wr(do.call(rbind,audit),file.path(out,'tables/human_layer_audit.tsv'))
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
cat('ANALYSIS_LOGIC_PASS HUMAN_RAW_COUNTS_AND_MAPPING\n')

CODE_EXPORT
cat > "$OUT/preflight.py" <<'CODE_PREFLIGHT'
import os, json, importlib, importlib.util
from pathlib import Path
for variable,expected in [('PREV1','mechanism_step1'),('PREV2','mechanism_step2_human_and_driver_audit'),('PREV8','mechanism_step8_human_mac_input_audit'),('SCENIC1','lung_scenic1_network'),('SCENIC2','lung_scenic2_activity_sensitivity')]:
    p=Path(os.environ[variable]);s=json.loads((p/'status.json').read_text())
    assert s['status']=='PASS' and s['stage']==expected,(variable,s)
for variable,files in [
 ('PREV1',['comm_input/CD4_counts.mtx.gz','comm_input/CD4_genes.tsv','comm_input/CD4_cells.tsv','tables/module_definition.tsv']),
 ('PREV2',['input/human_counts.mtx','input/human_genes.tsv','input/human_cells.tsv']),
 ('PREV8',['objects/marker_audit_checkpoint.rds'])]:
    for name in files:
        p=Path(os.environ[variable])/name;assert p.is_file() and p.stat().st_size>0,p
for arm in ['All','BG1']:
    for seed in [17,29,43]:
        p=Path(os.environ['SCENIC2'])/'objects'/f'Mouse_{arm}_{seed}_DropOverlap_auc.csv.gz'
        assert p.is_file() and p.stat().st_size>0,p
    for sp in ['Mouse','Human']:assert (Path(os.environ['SCENIC1'])/sp/arm/'cells.tsv').is_file()
versions={}
for pkg in ['numpy','pandas','scipy','scanpy','anndata','matplotlib']:
    mod=importlib.import_module(pkg);versions[pkg]=getattr(mod,'__version__','NA')
assert importlib.util.find_spec('leidenalg') or importlib.util.find_spec('louvain'),'Clustering dependency unavailable'
(Path(os.environ['OUT'])/'python_versions.json').write_text(json.dumps(versions,indent=2))
print('ANALYSIS_LOGIC_PASS FROZEN_INPUTS_AND_DEPENDENCIES',flush=True)

CODE_PREFLIGHT
cat > "$OUT/analyze.py" <<'CODE_ANALYZE'
import os, json, tarfile, hashlib, importlib.util
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.io import mmread
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(os.environ['OUT']);P1=Path(os.environ['PREV1']);P2=Path(os.environ['PREV2']);S1=Path(os.environ['SCENIC1']);S2=Path(os.environ['SCENIC2'])
SEEDS=[17,29,43]
RECEPTORS=['IL2RB','IL2RG','IL15RA','IL18R1','IL18RAP','CXCR6','TNFRSF1A','TNFRSF1B']
PANEL='CD3D CD3E TRAC CD4 CD8A CD8B1 TRDC TRGC1 NKG7 KLRD1 FOXP3 IL2RA CTLA4 IKZF2 IKZF4 LRRC32 ENTPD1 TIGIT CCR7 SELL TCF7 LEF1 IL7R GZMA GZMB PRF1 IFNG IL1R2 STAT1 ETV7 LST1 TYROBP FCER1G MS4A1 EPCAM KRT19 SFTPC'.split()
PROGRAMS={
 'AM':'FABP4 PPARG MARCO INHBA NAPSA',
 'C1qMac':'C1QA C1QB C1QC APOE MERTK',
 'SPP1Mac':'SPP1 GPNMB TREM2 CTSK CD9',
 'FOLR2Mac':'FOLR2 LYVE1 MRC1 SELENOP',
 'Mono':'FCN1 S100A8 S100A9 VCAN CD14',
 'cDC2':'CD1C FCER1A CLEC10A',
 'cDC1':'CLEC9A XCR1 BATF3',
 'pDC':'GZMB JCHAIN GZMB TCF4 CLEC4C',
 'Mast':'TPSAB1 TPSB2 KIT CPA3',
 'Neut':'FCGR3B CSF3R CXCR2',
 'Tcell':'CD3D CD3E TRAC',
 'Bcell':'MS4A1 CD79A CD79B',
 'Epi':'EPCAM KRT19 KRT8 KRT18 SFTPC',
 'Stroma':'COL1A1 COL1A2 DCN COL3A1',
 'PanMyel':'TYROBP LST1 FCER1G CTSS',
 'EpiBG':'SFTPA1 SFTPA2 SFTPB SFTPC SFTPD SCGB1A1 SCGB3A1 SCGB3A2 ABCA3 SLC34A2 NAPSA EPCAM KRT19'}
TARGETS=set('IL1R2 PENK TTC39C FOXP1 TNFRSF9 LAG3 MRE11 MRE11A IKZF2 STAT1 ETV7'.split())
def rd(p):return pd.read_csv(p,sep='\t',low_memory=False)
def wr(d,name):pd.DataFrame(d).to_csv(OUT/'tables'/name,sep='\t',index=False)
def project(x,genes):
    genes=np.array([str(g).upper() for g in genes]);assert all(genes!='NAN')
    unique,idx=np.unique(genes,return_inverse=True)
    proj=sparse.csr_matrix((np.ones(len(genes)),(np.arange(len(genes)),idx)),shape=(len(genes),len(unique)))
    return (x@proj).tocsr(),unique
def load(prefix,gcol):
    x=mmread(str(prefix)+'_counts.mtx'+('.gz' if Path(str(prefix)+'_counts.mtx.gz').exists() else '')).T.tocsr()
    g=rd(str(prefix)+'_genes.tsv');m=rd(str(prefix)+'_cells.tsv')
    assert x.shape==(len(m),len(g)) and m.cell.is_unique
    assert m[['cell','library']].notna().all().all() and np.isfinite(x.data).all() and (x.data>=0).all() and np.allclose(x.data,np.rint(x.data))
    assert (m.library>0).all() and np.all(np.asarray(x.sum(axis=1)).ravel()<=m.library.to_numpy()+1e-6)
    x,genes=project(x,g[gcol].astype(str));return x,genes,m
def gene_values(x,g,m,gene,kind='norm'):
    hit=np.flatnonzero(g==gene)
    if not len(hit):return np.full(len(m),np.nan)
    v=x[:,hit[0]].toarray().ravel()
    return v if kind=='raw' else v*1e4/m.library.to_numpy()
def signature(x,g,m,genes):
    genes=sorted(set(genes)&set(g));assert genes
    v=np.log1p(np.column_stack([gene_values(x,g,m,z) for z in genes]));mu=v.mean(axis=0);sd=v.std(axis=0)
    z=np.divide(v-mu,sd,out=np.zeros_like(v),where=sd>0)
    return z.mean(axis=1),genes
def summ_expr(x,g,m,genes,keys):
    values={z:gene_values(x,g,m,z) for z in genes};rows=[]
    for key,v in m.groupby(keys,observed=True,dropna=False):
        key=key if isinstance(key,tuple) else (key,);idx=v.index.to_numpy();base=dict(zip(keys,key))
        for gene in genes:
            y=values[gene][idx];present=gene in set(g)
            rows.append(dict(**base,gene=gene,cells=len(v),eligible20=len(v)>=20,mean_CP10k=float(np.mean(y)) if present else np.nan,pct=float(np.mean(y>0)*100) if present else np.nan,present=present))
    return pd.DataFrame(rows)
def fit(y,d,covariates):
    y=np.asarray(y,float);case=d.group.isin(['Silica','CWP']).astype(float).to_numpy()
    X=np.column_stack([np.ones(len(d)),case]+[d[c].to_numpy(float) for c in covariates])
    r=dict(n=len(d),n_control=int((case==0).sum()),n_case=int((case==1).sum()),coef=np.nan,p=np.nan,df=np.nan,reason='')
    if not np.isfinite(X).all() or not np.isfinite(y).all():r['reason']='NONFINITE';return r
    if min(r['n_control'],r['n_case'])<2:r['reason']='FEWER_THAN_TWO_PER_GROUP';return r
    # Standardize covariates for numerical stability; retain intercept and binary group.
    for j in range(2,X.shape[1]):
        sd=X[:,j].std()
        if sd>0:X[:,j]=(X[:,j]-X[:,j].mean())/sd
    rank=np.linalg.matrix_rank(X);r['df']=len(y)-rank
    if rank<X.shape[1] or r['df']<=0:r['reason']='UNIDENTIFIABLE';return r
    beta=np.linalg.lstsq(X,y,rcond=None)[0];r['coef']=beta[1]
    se=np.sqrt(max(0,float(np.sum((y-X@beta)**2)/r['df'])*np.linalg.inv(X.T@X)[1,1]))
    if se>1e-15:r['p']=2*stats.t.sf(abs(beta[1]/se),r['df'])
    else:r['reason']='ZERO_RESIDUAL_VARIANCE'
    return r
def fit_loo(y,d,covariates):
    r=fit(y,d,covariates);vals=[]
    for j in range(len(d)):
        k=np.arange(len(d))!=j;v=fit(np.asarray(y)[k],d.iloc[k],covariates)['coef']
        if np.isfinite(v):vals.append(v)
    r.update(loo_valid=len(vals),loo_min=min(vals) if vals else np.nan,loo_max=max(vals) if vals else np.nan);return r
def frozen(x,g,m,species,arm):
    a=rd(S1/species/arm/'cells.tsv');idx=pd.Index(m.cell).get_indexer(a.cell);assert (idx>=0).all()
    z=m.iloc[idx].reset_index(drop=True)
    assert np.array_equal(z.cell,a.cell) and np.array_equal(z['sample'],a['sample']) and np.array_equal(z.group,a.group) and np.allclose(z.library,a.library)
    z['subtype']=a.subtype.to_numpy();z['bg']=a.bg.to_numpy();return x[idx],z
def receptor_summary(x,g,m,species,arm):
    res=[];co=[]
    for scope in ['Core']+sorted(m.subtype.unique()):
        idx=np.arange(len(m)) if scope=='Core' else np.flatnonzero(m.subtype.eq(scope))
        z=m.iloc[idx].reset_index(drop=True);xx=x[idx]
        t=summ_expr(xx,g,z,RECEPTORS,['sample','group','subtype'] if scope!='Core' else ['sample','group'])
        res.append(t.assign(species=species,filter=arm,scope=scope))
        for name,genes in [('IL15_rec',['IL2RB','IL2RG']),('IL18_rec',['IL18R1','IL18RAP'])]:
            available=all(v in g for v in genes)
            det=np.logical_and.reduce([gene_values(xx,g,z,v,'raw')>0 for v in genes]) if available else None
            for ss,v in z.groupby('sample',observed=True):
                co.append(dict(species=species,filter=arm,scope=scope,sample=ss,group=v.group.iloc[0],complex=name,cells=len(v),eligible20=len(v)>=20,present=available,pct=float(det[v.index].mean()*100) if available else np.nan))
    return pd.concat(res,ignore_index=True),pd.DataFrame(co)
def run_mouse():
    x,g,m=load(P1/'comm_input/CD4','canonical')
    frozenmods=rd(P1/'tables/module_definition.tsv')
    treg=[z for z in frozenmods.loc[frozenmods.module.eq('TReg'),'gene'] if z not in TARGETS]
    t,used=signature(x,g,m,treg);m['TregNoTF']=t
    wr([dict(program='TregNoTF',gene=z) for z in used],'mouse_identity_signature.tsv')
    det={z:gene_values(x,g,m,z,'raw') for z in ['FOXP3','IL2RA','IL1R2']}
    m['FOXP3_det']=np.where(np.isnan(det['FOXP3']),np.nan,det['FOXP3']>0)
    m['FOXP3_IL2RA_det']=np.where(np.isnan(det['FOXP3'])|np.isnan(det['IL2RA']),np.nan,(det['FOXP3']>0)&(det['IL2RA']>0))
    marker=[];rec=[];co=[];effects=[];cellouts=[];strata=[];sampleouts=[]
    for arm in ['All','BG1']:
        xx,z=frozen(x,g,m,'Mouse',arm)
        marker.append(summ_expr(xx,g,z,PANEL,['sample','group','subtype']).assign(filter=arm))
        r,c=receptor_summary(xx,g,z,'Mouse',arm);rec.append(r);co.append(c)
        for seed in SEEDS:
            a=pd.read_csv(S2/'objects'/f'Mouse_{arm}_{seed}_DropOverlap_auc.csv.gz',index_col=0)
            assert a.index.is_unique and set(a.index)==set(z.cell)
            a=a.loc[z.cell];v=z.copy()
            for tf in ['Ikzf2','Stat1']:
                reg=tf+'(+)';assert reg in a
                v[tf+'_auc']=a[reg].to_numpy()
            v['filter']=arm;v['seed']=seed;cellouts.append(v)
            act=v[v.subtype.eq('ActCD4')].copy()
            d=act.groupby('sample',observed=True).agg(group=('group','first'),cells=('cell','size'),library=('library','mean'),TregNoTF=('TregNoTF','mean'),FOXP3_det=('FOXP3_det','mean'),FOXP3_IL2RA_det=('FOXP3_IL2RA_det','mean'),Ikzf2_auc=('Ikzf2_auc','mean'),Stat1_auc=('Stat1_auc','mean')).reset_index()
            d=d[d.cells>=20].reset_index(drop=True);d['loglib']=np.log1p(d.library);sampleouts.append(d.assign(filter=arm,seed=seed))
            for tf in ['Ikzf2','Stat1']:
                for name,cov in [('Group',[]),('Depth',['loglib']),('TregAdj',['TregNoTF']),('DepthTreg',['loglib','TregNoTF'])]:
                    effects.append(dict(filter=arm,seed=seed,scope='ActCD4',tf=tf,model=name,**fit_loo(d[tf+'_auc'].to_numpy(),d,cov)))
            # Detection strata are descriptions, not new cell-type definitions or exclusion gates.
            act['FoxDet']=np.where(act.FOXP3_det==1,'Detected','NotDetect')
            for (ss,st),vv in act.groupby(['sample','FoxDet'],observed=True):
                strata.append(dict(filter=arm,seed=seed,sample=ss,group=vv.group.iloc[0],stratum=st,cells=len(vv),Ikzf2_auc=vv.Ikzf2_auc.mean(),Stat1_auc=vv.Stat1_auc.mean(),TregNoTF=vv.TregNoTF.mean()))
    wr(pd.concat(marker),'mouse_marker_sample.tsv');wr(effects,'mouse_TF_identity_adjustment.tsv');wr(strata,'mouse_FOXP3_strata.tsv');wr(pd.concat(sampleouts),'mouse_identity_sample.tsv')
    pd.concat(cellouts).to_csv(OUT/'tables/mouse_identity_cells.tsv.gz',sep='\t',index=False,compression='gzip')
    return rec,co
def run_human_receptors():
    x,g,m=load(P2/'input/human','symbol');rec=[];co=[]
    for arm in ['All','BG1']:
        xx,z=frozen(x,g,m,'Human',arm);r,c=receptor_summary(xx,g,z,'Human',arm);rec.append(r);co.append(c)
    return rec,co
def plot_umap(a,key,name):
    xy=a.obsm['X_umap'];labs=a.obs[key].fillna('Unknown').astype(str);levels=sorted(labs.unique());fig,ax=plt.subplots(figsize=(8,6))
    cmap=plt.get_cmap('tab20')
    for j,lab in enumerate(levels):
        k=labs.eq(lab).to_numpy();label=lab if len(lab)<=10 else 'L'+str(j)
        ax.scatter(xy[k,0],xy[k,1],s=2,c=[cmap(j%20)],label=label,rasterized=True)
    ax.set(xlabel='UMAP1',ylabel='UMAP2',title=name);ax.legend(markerscale=4,bbox_to_anchor=(1.01,1),loc='upper left',frameon=False,fontsize=7)
    fig.tight_layout();fig.savefig(OUT/'figures'/f'{name}.png',dpi=300);fig.savefig(OUT/'figures'/f'{name}.pdf');plt.close(fig)
    return [dict(figure=name,label=(lab if len(lab)<=10 else 'L'+str(j)),original=lab) for j,lab in enumerate(levels)]
def run_mac():
    import scanpy as sc
    import anndata as ad
    x,g,m=load(OUT/'input/human_mac','gene');assert np.allclose(np.asarray(x.sum(axis=1)).ravel(),m.library)
    for k in ['sample','group','series']:m[k]=m[k].fillna('Unknown').astype(str)
    m['parent']=m.parent.astype(str);obs=m.copy();obs.index=obs.cell
    a=ad.AnnData(x.astype('float32'),obs=obs,var=pd.DataFrame(index=g));a.layers['counts']=x
    a.X=x.multiply(1e4/m.library.to_numpy()[:,None]).tocsr().astype('float32');sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a,flavor='seurat',n_top_genes=3000,batch_key='sample')
    excluded=np.array([z in TARGETS or z.startswith(('MT-','RPL','RPS')) for z in g]);a.var['highly_variable'] &= ~excluded
    assert a.var.highly_variable.sum()>30
    b=a[:,a.var.highly_variable].copy();sc.pp.scale(b,max_value=10);sc.tl.pca(b,n_comps=30,svd_solver='arpack',random_state=17)
    a.obsm['X_pca']=b.obsm['X_pca'];del b
    sc.pp.neighbors(a,n_neighbors=15,n_pcs=30,random_state=17);sc.tl.umap(a,random_state=17)
    algorithm='leiden' if importlib.util.find_spec('leidenalg') else 'louvain'
    assert importlib.util.find_spec('leidenalg') or importlib.util.find_spec('louvain')
    for res in [0.2,0.4,0.8]:
        key='r'+str(int(res*10));getattr(sc.tl,algorithm)(a,resolution=res,random_state=17,key_added=key)
    sc.tl.rank_genes_groups(a,groupby='r4',method='wilcoxon',use_raw=False,n_genes=40,pts=True)
    sc.get.rank_genes_groups_df(a,group=None).to_csv(OUT/'tables/human_cluster_markers.tsv',sep='\t',index=False)
    m=a.obs.reset_index(drop=True).copy();coverage=[]
    for name,gs in PROGRAMS.items():
        score,used=signature(x,g,m,gs.split());m[name]=score;coverage.extend([dict(program=name,gene=z,present=z in g) for z in sorted(set(gs.split()))])
    epigenes=sorted(set(PROGRAMS['EpiBG'].split())&set(g));m['epi_fraction']=np.asarray(x[:,np.isin(g,epigenes)].sum(axis=1)).ravel()/m.library.to_numpy()
    m['IL15_IL15RA_det']=((gene_values(x,g,m,'IL15','raw')>0)&(gene_values(x,g,m,'IL15RA','raw')>0)).astype(float) if all(z in g for z in ['IL15','IL15RA']) else np.nan
    rows=[];mix=[]
    for resolution in ['r2','r4','r8']:
        for cluster,v in m.groupby(resolution,observed=True):
            base=dict(resolution=resolution,cluster=str(cluster),cells=len(v),mapped_cells=int(v['sample'].ne('Unknown').sum()),subjects=v.loc[v['sample'].ne('Unknown'),'sample'].nunique(),dominant_sample_fraction=v['sample'].value_counts(normalize=True).max(),epi_fraction_mean=v.epi_fraction.mean())
            means={p:float(v[p].mean()) for p in PROGRAMS}
            ordered=sorted([p for p in PROGRAMS if p not in ['PanMyel','EpiBG']],key=lambda k:means[k],reverse=True)
            rows.append(dict(**base,program_rank1=ordered[0],program_rank2=ordered[1],**means))
            for (ss,gr,se,par),vv in v.groupby(['sample','group','series','parent'],observed=True):mix.append(dict(resolution=resolution,cluster=str(cluster),sample=ss,group=gr,series=se,parent=par,cells=len(vv)))
    wr(rows,'human_cluster_scorecard.tsv');wr(mix,'human_cluster_sample_composition.tsv');wr(coverage,'human_marker_coverage.tsv')
    genes=sorted(set(sum([v.split() for v in PROGRAMS.values()],[])+['IL15','IL15RA','TNF','IL18','CXCL16','IL1R2']))
    wr(summ_expr(x,g,m,genes,['r4']),'human_cluster_marker_expression.tsv')
    allsend=summ_expr(x,g,m,['IL15','IL15RA','TNF','IL18','CXCL16'],['r4','sample','group','series'])
    wr(allsend,'human_sender_cluster_expression.tsv')
    sen=[]
    # Paired sensitivity within each refined cluster excludes only the preflagged parent16;
    # does not remove cells using a newly chosen expression threshold.
    for arm in ['All','NoParent16']:
        keep=np.ones(len(m),bool) if arm=='All' else m.parent.ne('16').to_numpy()
        zz=m.loc[keep].reset_index(drop=True);xx=x[keep]
        t=summ_expr(xx,g,zz,['IL15','IL15RA','TNF','IL18','CXCL16'],['r4','sample','group','series']);t['sender_filter']=arm;sen.append(t)
    wr(pd.concat(sen),'human_sender_parent_sensitivity.tsv')
    comps=m.groupby(['r4','sample','group','series'],observed=True).agg(cells=('cell','size'),library=('library','mean'),epi_fraction=('epi_fraction','mean'),IL15_IL15RA_pct=('IL15_IL15RA_det',lambda v:100*v.mean())).reset_index()
    wr(comps,'human_sender_sample_covariates.tsv')
    # Save all input genes and all candidate cells for reviewed annotation; program ranks are not final labels.
    a.obs=m.set_index('cell',drop=False);a.uns['clustering_algorithm']=algorithm;a.write_h5ad(OUT/'objects/human_myeloid_refined.h5ad',compression='gzip')
    wr(m,'human_refined_cells.tsv.gz')
    labels=[]
    for key,name in [('r4','MacClust'),('parent','MacParent'),('group','MacGroup'),('series','MacSource')]:labels.extend(plot_umap(a,key,name))
    wr(labels,'figure_label_map.tsv')
    (OUT/'figures/captions.txt').write_text('MacClust: human candidate myeloid cells, unsupervised resolution 0.4 clusters; not final macrophage annotations. MacParent: original broad clusters 3/6/16/22. MacGroup: Control/CWP/Unknown. MacSource: original data series. UMAP1/UMAP2: embedding axes; same embedding across figures. Full label mapping: tables/figure_label_map.tsv. Clustering algorithm: '+algorithm+'. No sample integration or group regression applied to clustering; sample/source composition is exported for review.\n')
    print('HUMAN_MYeloid_REFINEMENT_PASS cells='+str(len(m)),flush=True)
def main():
    for p in ['tables','objects','figures','input']:(OUT/p).mkdir(exist_ok=True)
    rec,co=run_mouse();r,c=run_human_receptors();rec+=r;co+=c
    wr(pd.concat(rec),'CD4_receptor_detection.tsv');wr(pd.concat(co),'CD4_receptor_codetection.tsv')
    print('MOUSE_IDENTITY_AND_RECEPTOR_PASS',flush=True)
    run_mac()
    (OUT/'methods.txt').write_text('Frozen SCENIC1 cells and SCENIC2 DropOverlap AUC; no GRN rerun. Mouse identity: original TReg module excluding tested TF/candidate genes, gene-wise logCP10k z scores frozen over all original CD4 cells; same scale across All/BG1. FOXP3 detection strata are descriptions, not Treg labels. Sample-level ActCD4 group models with/without depth and TregNoTF; exact residual df and LOO reported. n=6 is exploratory, no multiplicity-confirmed claim. Human raw RNA myeloid candidates are original 0.4 clusters 3/6/16/22; all cells retained, unknown subjects retained for annotation only. Selection follows preceding marker audit, not disease or ligand expression. HVGs selected by sample; remove candidate/TF/MT/RPL/RPS from clustering features only. LogCP10k; PCA30, neighbors15; resolutions .2/.4/.8; no integration. Cluster markers use cell-level Wilcoxon for annotation, not subject-level disease inference. Program ranking is exploratory and never a final cell label. Whole sender raw-count data checkpoint retained. Parent16 exclusion is a paired sensitivity because pre-audit showed epithelial co-signal, not a proven doublet filter. Epithelial burden retained continuously. Per-subject receptor expression and IL2RB/IL2RG or IL18R1/IL18RAP co-detection exported; transcript co-detection does not establish receptor function. Ligand IL15/IL15RA co-detection evaluated separately on human sender. Frozen >=20-cell flag retained, no new cell-expression cutoff.\n')
    (OUT/'status.json').write_text(json.dumps(dict(status='PASS',stage='lung_scenic3_identity_sender_refinement',output=str(OUT),scenic1=str(S1),scenic2=str(S2),pending=['review_human_myeloid_subtype_assignments','final_sender_TF_models_using_reviewed_labels']),indent=2))
    ret=Path(os.environ['ROOT'])/'tmp/lung_scenic3_return.tar.gz';part=Path(str(ret)+'.part')
    with tarfile.open(part,'w:gz') as t:
        for p in sorted(OUT.rglob('*')):
            if p.is_file() and not {'input','objects','__pycache__'}&set(p.relative_to(OUT).parts):t.add(p,arcname=str(p.relative_to(OUT)))
    with tarfile.open(part) as t:
        for m in t:
            if m.isfile():
                with t.extractfile(m) as f:
                    while f.read(2**20):pass
    part.replace(ret);print('RETURN='+str(ret),flush=True);print('RETURN_SHA256='+hashlib.sha256(ret.read_bytes()).hexdigest(),flush=True);print('LUNG_SCENIC3_PASS',flush=True)
if __name__=='__main__':main()

CODE_ANALYZE
conda run --no-capture-output -n batch_bench python -m py_compile "$OUT/preflight.py" "$OUT/analyze.py"
conda run --no-capture-output -n R_4.5.2 Rscript -e 'invisible(parse(file=commandArgs(TRUE)[1]));cat("R_PARSE_PASS\n")' "$OUT/export.R"
echo SYNTAX_CHECK=PASS
conda run --no-capture-output -n batch_bench python -u "$OUT/preflight.py"
conda run --no-capture-output -n R_4.5.2 Rscript "$OUT/export.R" "$PREV8" "$OUT"
conda run --no-capture-output -n batch_bench python -u "$OUT/analyze.py"
