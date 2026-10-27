#!/usr/bin/env bash
set -Eeuo pipefail
ROOT="${ROOT:-/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA}"
PREV="${PREV:-$ROOT/02_analysis/mechanism_step1_20260927_152218_s4V5eF}"
export ROOT PREV
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-2}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
export MPLBACKEND=Agg
bash -n "$0"
mkdir -p "$ROOT/tmp" "$ROOT/logs" "$ROOT/02_analysis"
exec 9>"$ROOT/tmp/lung_mechanism_step2.lock"
flock -n 9 || { echo 'MECHANISM_STEP2_ALREADY_RUNNING'; exit 1; }
OUT=$(mktemp -d "$ROOT/02_analysis/mechanism_step2_$(date +%Y%m%d_%H%M%S)_XXXXXX")
export OUT
trap 'rc=$?; if ((rc != 0)); then printf "LUNG_MECHANISM_STEP2_FAIL EXIT_CODE=%s OUTPUT=%s\n" "$rc" "$OUT"; fi' EXIT
mkdir -p "$OUT/input" "$OUT/tables" "$OUT/objects" "$OUT/figures"
cp -- "$0" "$OUT/run.sh"
cat > "$OUT/shared.py" <<'PY_SHARED'
import os, re, json, gzip, tarfile, hashlib, platform, itertools
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata
from scipy.io import mmwrite

TARGETS = ['IL1R2','PENK','TTC39C','FOXP1','TNFRSF9','LAG3','MRE11']
SETS = {
 'OXPHOS':'NDUFA1 NDUFA2 NDUFA9 NDUFB8 NDUFS1 NDUFS2 NDUFS3 NDUFV1 SDHA SDHB SDHC UQCRC1 UQCRC2 UQCRFS1 CYC1 COX4I1 COX5A COX5B COX6A1 ATP5F1A ATP5F1B ATP5F1C ATP5MC1 ATP5PO',
 'MitoISR':'ATF4 ATF5 DDIT3 DDIT4 CHAC1 ASNS TRIB3 PPP1R15A SESN2 DELE1 EIF2AK1 HSPA9 HSPD1 HSPE1 LONP1 CLPP OMA1',
 'MitoROS':'SOD2 PRDX3 PRDX5 GPX4 TXN2 TXNRD2 HMOX1 NQO1 GCLC GCLM SRXN1 FTH1',
 'Mitophagy':'PINK1 PRKN BNIP3 BNIP3L FUNDC1 SQSTM1 OPTN TBK1 ULK1 BECN1 MAP1LC3B',
 'Inflam':'NLRP3 AIM2 NLRC4 PYCARD CASP1',
 'PyroExec':'PYCARD CASP1 GSDMD IL1B IL18',
 'ApoExec':'BAX BAK1 APAF1 CASP9 CASP3 CASP7 CYCS',
 'TAct':'CD69 IL2RA ICOS TNFRSF4 BATF IRF4 NR4A1 NR4A2 EGR1',
 'TReg':'FOXP3 IL2RA CTLA4 IKZF2 IKZF4 LRRC32 ENTPD1 TIGIT',
 'Inhib':'PDCD1 CTLA4 TIGIT HAVCR2 TOX TOX2 ENTPD1',
 'Naive':'CCR7 SELL TCF7 LEF1 KLF2 IL7R MAL',
 'Cycle':'MKI67 TOP2A PCNA MCM2 MCM3 MCM4 MCM5 MCM6 MCM7 UBE2C'
}
SETS={k:v.split() for k,v in SETS.items()}
EPI='SFTPA1 SFTPA2 SFTPB SFTPC SFTPD SCGB1A1 SCGB3A1 SCGB3A2 ABCA3 SLC34A2 NAPSA'.split()
MYE='LYZ LYZ1 LYZ2 S100A8 S100A9 CSF1R LST1 TYROBP FCER1G'.split()
SOURCE_GENES='IL1A IL1B IL18 CXCL16 TNF IL6 IL10 TGFB1 IL12A IL12B IL23A CCL3 CCL4 CCL5 SPP1 LGALS9 CD80 CD86 CD274 PDCD1LG2 TNFSF9'.split()
RECEPTORS='IL1R1 IL1RAP IL1R2 IL18R1 IL18RAP CXCR6 TNFRSF1A TNFRSF1B IL6R IL6ST IL10RA IL10RB TGFBR1 TGFBR2 IL12RB1 IL12RB2 IL23R CCR5 CD44 ITGAV ITGB1 ITGB3 HAVCR2 CD28 CTLA4 PDCD1 TNFRSF9'.split()
SPECS=[
 ('GSE314943','CD4','GSE314943/stage5_20260905_222457_a2dApT/objects/CD4_core.h5ad',2311),
 ('GSE314943','Mac','GSE314943/stage11_20260906_100314_pALwAB/objects/macrophages_core.h5ad',1796),
 ('GSE183682','CD4','screen2_CD4_target_wz1lfx/GSE183682/CD4_core.h5ad',1413),
 ('GSE274372','TLike','GSE274372/stage9_20260908_140721_Yj56oq/objects/CD4_author_reclustered.h5ad',888),
 ('GSE174725','CD4','screen2_Trefine_nuzVQu/GSE174725/T_refined.h5ad',10286)]

def canonical(s):
    s=str(s).strip().upper()
    return {'MRE11A':'MRE11','ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','PARK2':'PRKN'}.get(s,s)

def rho(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float);ok=np.isfinite(a)&np.isfinite(b);a=a[ok];b=b[ok]
    if len(a)<3 or np.ptp(a)==0 or np.ptp(b)==0:return np.nan
    return float(np.clip(np.corrcoef(rankdata(a),rankdata(b))[0,1],-1,1))

def partial(a,b,z):
    a=np.asarray(a,float);b=np.asarray(b,float);z=np.asarray(z,float)
    if z.ndim==1:z=z[:,None]
    ok=np.isfinite(a)&np.isfinite(b)&np.isfinite(z).all(axis=1);a=a[ok];b=b[ok];z=z[ok]
    if len(a)<6 or np.ptp(a)==0 or np.ptp(b)==0:return np.nan
    d=np.column_stack([np.ones(len(a))]+[rankdata(z[:,j]) for j in range(z.shape[1])])
    if len(a)-np.linalg.matrix_rank(d)<3:return np.nan
    u=rankdata(a);v=rankdata(b);u-=d@np.linalg.lstsq(d,u,rcond=None)[0];v-=d@np.linalg.lstsq(d,v,rcond=None)[0]
    if np.std(u)<1e-10 or np.std(v)<1e-10:return np.nan
    return float(np.clip(np.corrcoef(u,v)[0,1],-1,1))

def effect(b,c):
    b=np.asarray(b,float);c=np.asarray(c,float);b=b[np.isfinite(b)];c=c[np.isfinite(c)]
    d=float(c.mean()-b.mean()) if len(b) and len(c) else np.nan
    loo=[]
    if min(len(b),len(c))>=2:
        loo=[c.mean()-np.delete(b,i).mean() for i in range(len(b))]+[np.delete(c,i).mean()-b.mean() for i in range(len(c))]
    return dict(n_control=len(b),n_case=len(c),control_mean=b.mean() if len(b) else np.nan,case_mean=c.mean() if len(c) else np.nan,
      difference=d,loo_same=np.mean(np.sign(loo)==np.sign(d)) if loo and d!=0 else np.nan,
      loo_min=min(loo) if loo else np.nan,loo_max=max(loo) if loo else np.nan)

def bgmask(symbols):
    return np.array([bool(re.match(r'^(IGH[VDJMAEG]|IGK[VJC]|IGL[VJC])',s)) or s in EPI or s in MYE for s in symbols])

def scores(y, x, symbols, samples):
    # Equal-sample reference expression; controls fixed before group comparisons.
    reference=np.mean([np.asarray(y[np.asarray(samples)==s].mean(axis=0)).ravel() for s in sorted(set(samples))],axis=0)
    excluded=set(TARGETS)|set(itertools.chain.from_iterable(SETS.values()))|set(EPI)|set(MYE)
    pool=np.array([s not in excluded and not re.match(r'^(MT-|RPL|RPS|MRPL|MRPS|IGH|IGK|IGL)',s) for s in symbols])&(reference>0)
    edges=np.unique(np.quantile(reference[reference>0],np.linspace(0,1,26)))
    bins=np.digitize(reference,edges[1:-1]);coverage=[];controls=[];values={};members={}
    for name,expected in SETS.items():
        j=np.flatnonzero(np.isin(symbols,expected));observed=j[np.asarray(x[:,j].sum(axis=0)).ravel()>0]
        present=list(symbols[j]);eligible=len(j)>=3 and len(j)/len(expected)>=.6 and len(observed)>=3
        selected=set()
        for ix in j:
            same=np.flatnonzero(pool & (bins==bins[ix]))
            options=same if len(same)>=20 else np.flatnonzero(pool)
            order=np.argsort(np.abs(reference[options]-reference[ix]),kind='stable')
            selected.update(options[order[:20]].tolist())
        cc=np.array(sorted(selected),dtype=int)
        eligible=eligible and len(cc)>=20
        coverage.append(dict(module=name,expected=len(expected),present=len(j),expressed=len(observed),coverage=len(j)/len(expected),eligible=eligible,
          genes=';'.join(present),missing=';'.join(sorted(set(expected)-set(present))),controls=len(cc)))
        for ix in cc:controls.append(dict(module=name,gene=symbols[ix],reference_mean=reference[ix],bin=int(bins[ix])))
        mean=np.asarray(y[:,j].mean(axis=1)).ravel() if len(j) else np.full(y.shape[0],np.nan)
        values[name]=mean-np.asarray(y[:,cc].mean(axis=1)).ravel() if eligible else np.full(y.shape[0],np.nan)
        members[name]=np.asarray((x[:,j]>0).sum(axis=1)).ravel() if len(j) else np.zeros(y.shape[0])
    return values,members,coverage,controls


PY_SHARED
cat > "$OUT/export.R" <<'R_EXPORT'
args <- commandArgs(TRUE)
root <- args[1]; out <- args[2]
suppressPackageStartupMessages({library(Seurat); library(Matrix)})
p <- file.path(root,'06_CD4T_analysis/CD4_final/GSE288887_CD4_final.rds')
stopifnot(file.exists(p))
a <- readRDS(p); DefaultAssay(a) <- 'RNA'
if(inherits(a[['RNA']], 'Assay5')) {
  lyr <- Layers(a[['RNA']], search='^counts')
  if(length(lyr)!=1L || lyr!='counts') a <- JoinLayers(a,assay='RNA')
}
x <- GetAssayData(a,assay='RNA',layer='counts')
x <- as(x,'CsparseMatrix'); md <- a[[]][colnames(x),,drop=FALSE]
stopifnot(identical(rownames(md),colnames(x)),!anyDuplicated(colnames(x)),!anyDuplicated(rownames(x)))
stopifnot(all(is.finite(x@x)),all(x@x>=0),all(abs(x@x-round(x@x))<1e-6),all(Matrix::colSums(x)>0))
pick <- function(candidates) {
  found <- intersect(candidates,colnames(md))
  if(!length(found)) stop('Missing metadata: ',paste(candidates,collapse=','))
  found[1]
}
ss <- as.character(md[[pick(c('subject_id','subject','Subject','donor_id','patient_id'))]])
ss <- toupper(ss); hits <- regexpr('(CTL|CWP)[0-9]+',ss)
stopifnot(all(hits>0)); ss <- regmatches(ss,hits)
group <- tolower(as.character(md[[pick(c('group','Group','condition','Condition'))]]))
group <- ifelse(grepl('cwp',group),'CWP',ifelse(grepl('control|ctl',group),'Control',NA_character_))
ty <- as.character(md[[pick(c('cd4_label','CD4_anno','cd4_anno','Subtype','subtype'))]])
g227 <- c('CWP85','CWP86','CTL092','CTL094','CTL098','CTL103','CTL105','CTL106','CTL108','CTL111','CTL117')
g135 <- c('CTL068','CTL070');g288 <- c('CWP89','CWP94','CWP109','CWP127','CWP144','CWP149','CWP152')
series <- ifelse(ss%in%g227,'G227',ifelse(ss%in%g135,'G135',ifelse(ss%in%g288,'G288',NA_character_)))
stopifnot(!anyNA(series),!anyNA(group),!anyNA(ty),all(nzchar(ty)))
stopifnot(all((grepl('^CTL',ss)&group=='Control')|(grepl('^CWP',ss)&group=='CWP')))
m <- data.frame(cell=colnames(x),sample=ss,group=group,series=series,subtype=ty,core=ty!='CD4state',library=Matrix::colSums(x))
u <- unique(m[,c('sample','group','series')]);stopifnot(!anyDuplicated(u$sample))
write.table(m,file.path(out,'input/human_cells.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(data.frame(symbol=rownames(x)),file.path(out,'input/human_genes.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
writeMM(x,file.path(out,'input/human_counts.mtx'))
write.table(u,file.path(out,'tables/human_subject_design.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(as.data.frame(table(m$subtype,m$sample)),file.path(out,'tables/human_subtype_counts.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(data.frame(package=c('Seurat','Matrix','nichenetr','multinichenetr'),available=vapply(c('Seurat','Matrix','nichenetr','multinichenetr'),requireNamespace,logical(1),quietly=TRUE)),file.path(out,'tables/package_readiness.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(c(p,paste('cells',ncol(x),'genes',nrow(x))),file.path(out,'input/human_source.txt'))
cat('HUMAN_EXPORT_PASS\n')

R_EXPORT
cat > "$OUT/analyze.py" <<'PY_ANALYZE'
import os,json,gzip,hashlib,itertools
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from shared import TARGETS,SETS,canonical,bgmask,scores,effect

root=Path(os.environ['ROOT']);out=Path(os.environ['OUT']);prev=Path(os.environ['PREV'])
def wr(df,name):df.to_csv(out/'tables'/name,sep='\t',index=False,na_rep='NA')
def matrix(path):
    with gzip.open(path,'rb') if str(path).endswith('.gz') else open(path,'rb') as h:return sparse.csr_matrix(mmread(h).T,dtype=float)
def check(x,m):
    x.sum_duplicates();x.eliminate_zeros()
    assert x.shape[0]==len(m) and m.cell.is_unique
    assert np.isfinite(x.data).all() and (x.data>=0).all() and np.allclose(x.data,np.rint(x.data),atol=1e-6,rtol=0)
    assert (m.library>0).all()
    assert (np.asarray(x.sum(axis=1)).ravel()<=m.library.to_numpy()+1e-6).all()

# Human counts retain original full-library denominators and exact subject/source labels.
hm=pd.read_csv(out/'input/human_cells.tsv',sep='\t');raw=matrix(out/'input/human_counts.mtx');check(raw,hm)
hm['core']=hm.core.astype(str).str.upper().map({'TRUE':True,'FALSE':False});assert hm.core.notna().all()
original=pd.read_csv(out/'input/human_genes.tsv',sep='\t').symbol.astype('string').fillna('')
sy=np.array([canonical(s) for s in original]);valid=(sy!='')&~np.char.startswith(sy,'ENS')
symbols=np.array(sorted(set(sy[valid])));assert len(symbols)>5000
proj=sparse.csr_matrix((np.ones(valid.sum()),(np.flatnonzero(valid),pd.Index(symbols).get_indexer(sy[valid]))),shape=(len(sy),len(symbols)))
hx=(raw@proj).tocsr();lib=hm.library.to_numpy()
wr(pd.DataFrame({'original':original,'canonical':sy,'retained':valid}),'human_gene_mapping.tsv.gz')
hm['bg']=np.asarray(hx[:,bgmask(symbols)].sum(axis=1)).ravel()/lib
hm['mt']=np.asarray(hx[:,np.char.startswith(symbols,'MT-')].sum(axis=1)).ravel()/lib
hm['nfeature']=np.asarray((hx>0).sum(axis=1)).ravel()
y=hx.multiply((1e4/lib)[:,None]).tocsr();y.data=np.log1p(y.data)
vals,det,cov,ctrl=scores(y,hx,symbols,hm['sample'].to_numpy())
wr(pd.DataFrame(cov),'human_module_coverage.tsv');wr(pd.DataFrame(ctrl),'human_score_controls.tsv.gz')
for mod in SETS:hm['score_'+mod]=vals[mod];hm['detect_'+mod]=det[mod]
for gene in TARGETS:
    j=np.flatnonzero(symbols==gene);a=np.asarray(hx[:,j].sum(axis=1)).ravel() if len(j) else np.full(len(hm),np.nan)
    hm['expr_'+gene]=a/lib*1e4;hm['positive_'+gene]=a>0 if len(j) else np.nan
rows=[]
for filt,cut in [('All',1.),('BG1',.01)]:
    for scope in ['Primary','Core']+['Subtype:'+s for s in sorted(hm.subtype.unique())]:
        mask=np.ones(len(hm),bool) if scope=='Primary' else hm.core.to_numpy() if scope=='Core' else (hm.subtype==scope.split(':',1)[1]).to_numpy()
        for ss in sorted(hm['sample'].unique()):
            allm=hm[hm['sample']==ss];z=hm[mask&(hm['sample']==ss)&(hm.bg<=cut)];n=len(z)
            d=dict(scope=scope,filter=filt,sample=ss,group=allm.group.iloc[0],series=allm.series.iloc[0],cells=n,original_cells=int(np.sum(mask&(hm['sample']==ss))),eligible20=n>=20,library=z.library.mean(),bg=z.bg.mean(),mt=z.mt.mean())
            for k in SETS:d['score_'+k]=z['score_'+k].mean();d['detect_'+k]=z['detect_'+k].mean()
            for k in TARGETS:d['expr_'+k]=z['expr_'+k].mean();d['pct_'+k]=z['positive_'+k].mean()*100
            rows.append(d)
ht=pd.DataFrame(rows);wr(ht,'human_sample_scores.tsv')
hm.to_csv(out/'objects/human_cell_scores.tsv.gz',sep='\t',index=False,na_rep='NA')
# Fixed subtype weights shared across groups/filters within each human scope.
std=[]
for sc in ['Primary','Core']:
    b=hm if sc=='Primary' else hm[hm.core]
    ct=pd.crosstab(b['sample'],b.subtype);weights=ct.div(ct.sum(axis=1),axis=0).mean()
    for filt in ['All','BG1']:
        for ss in sorted(b['sample'].unique()):
            z=ht[(ht['filter']==filt)&(ht['sample']==ss)&ht.scope.str.startswith('Subtype:')].copy();z['type']=z.scope.str.replace('Subtype:','',regex=False);z=z.set_index('type').reindex(weights.index)
            for metric in ['score_'+k for k in SETS]+['expr_'+g for g in TARGETS]:
                std.append(dict(scope=sc,filter=filt,sample=ss,group=z.group.iloc[0],series=z.series.iloc[0],metric=metric,value=float(z[metric]@weights) if z[metric].notna().all() else np.nan,all_subtypes_eligible20=bool((z.cells>=20).all()),weights=json.dumps(weights.to_dict())))
wr(pd.DataFrame(std),'human_fixed_composition.tsv')

# Decompose main mouse and human panel scores into measured gene contributions.
mm=pd.read_csv(prev/'comm_input/CD4_cells.tsv',sep='\t');mx=matrix(prev/'comm_input/CD4_counts.mtx.gz');check(mx,mm)
ms=pd.read_csv(prev/'comm_input/CD4_genes.tsv',sep='\t').canonical.astype(str).to_numpy();assert mx.shape[1]==len(ms)
mm['core']=True;mm['series']='G314';mm['bg']=np.asarray(mx[:,bgmask(ms)].sum(axis=1)).ravel()/mm.library.to_numpy()
rawrows=[];contrib=[];joint=[];lomo=[]
extra=['GSDME','CASP3','CASP1','GSDMD','PYCARD','IL1B','IL18','NLRP3']
for ds,x,meta,genes,ct in [('GSE314943',mx,mm,ms,pd.read_csv(prev/'tables/score_controls.tsv.gz',sep='\t').query("dataset == 'GSE314943'")),('HumanLung',hx,hm,symbols,pd.DataFrame(ctrl))]:
    log=x.multiply((1e4/meta.library.to_numpy())[:,None]).tocsr();log.data=np.log1p(log.data)
    targets=sorted(set(TARGETS+extra+list(itertools.chain.from_iterable(SETS.values()))));counts={}
    for g in targets:
        j=np.flatnonzero(genes==g);counts[g]=np.asarray(x[:,j].sum(axis=1)).ravel() if len(j) else np.full(len(meta),np.nan)
    available={g:np.any(genes==g) for g in targets}
    ctrls={k:np.flatnonzero(np.isin(genes,ct.loc[ct.module==k,'gene'])) for k in SETS}
    ctrlvalues={k:np.asarray(log[:,j].mean(axis=1)).ravel() if len(j) else np.full(len(meta),np.nan) for k,j in ctrls.items()}
    subtypes=sorted(meta.subtype.unique());scopes=['Primary','Core']+['Subtype:'+k for k in subtypes] if ds=='HumanLung' else ['Core']+['Subtype:'+k for k in subtypes]
    for filt,cut in [('All',1.),('BG1',.01)]:
        for sc in scopes:
            sm=np.ones(len(meta),bool) if sc=='Primary' else meta.core.to_numpy() if sc=='Core' else (meta.subtype==sc.split(':',1)[1]).to_numpy()
            for ss in sorted(meta['sample'].unique()):
                ix=np.flatnonzero(sm&(meta['sample']==ss)&(meta.bg<=cut));n=len(ix);base=meta[meta['sample']==ss].iloc[0]
                key=dict(dataset=ds,scope=sc,filter=filt,sample=ss,group=base.group,series=base.series,cells=n,eligible20=n>=20)
                for g in targets:
                    v=counts[g][ix];ok=available[g] and n>0
                    rawrows.append(dict(**key,gene=g,present=available[g],positive=int(np.sum(v>0)) if ok else np.nan,pct=100*np.mean(v>0) if ok else np.nan,CP10k=np.mean(v/meta.library.to_numpy()[ix]*1e4) if ok else np.nan,mean_log=np.mean(np.log1p(v/meta.library.to_numpy()[ix]*1e4)) if ok else np.nan))
                for label,gs in [('CaspGsd',['CASP1','GSDMD']),('Pyro3',['CASP1','GSDMD','PYCARD']),('AltPyro',['CASP3','GSDME'])]:
                    ok=all(available[g] for g in gs) and n>0
                    pos=np.logical_and.reduce([counts[g][ix]>0 for g in gs]) if ok else None
                    joint.append(dict(**key,pair=label,genes=';'.join(gs),positive=int(pos.sum()) if ok else np.nan,pct=100*pos.mean() if ok else np.nan))
                for mod in ['MitoISR','MitoROS','PyroExec','Inflam','OXPHOS']:
                    gs=[g for g in SETS[mod] if available[g]];nmod=len(gs)
                    if not n or not nmod:continue
                    bg=float(np.mean(ctrlvalues[mod][ix]));means={g:np.mean(np.log1p(counts[g][ix]/meta.library.to_numpy()[ix]*1e4)) for g in gs};full=float(np.mean(list(means.values()))-bg)
                    for g,v in means.items():
                        contrib.append(dict(**key,module=mod,gene=g,mean_log=v,control_mean=bg,contribution=(v-bg)/nmod,module_score=full,module_genes=nmod))
                        lomo.append(dict(**key,module=mod,removed_gene=g,full_score=full,reduced_score=(sum(means.values())-v)/(nmod-1)-bg if nmod>1 else np.nan))
wr(pd.DataFrame(rawrows),'panel_gene_expression.tsv');wr(pd.DataFrame(contrib),'module_gene_contributions.tsv');wr(pd.DataFrame(joint),'pyro_codetection.tsv');wr(pd.DataFrame(lomo),'module_leave_one_gene.tsv')
# Exact additive decomposition audit against previous mouse panel scores.
c=pd.DataFrame(contrib);prevtab=pd.read_csv(prev/'tables/sample_scores.tsv',sep='\t');errors=[]
for (scope,filt,ss,mod),v in c[c.dataset=='GSE314943'].groupby(['scope','filter','sample','module'],observed=True):
    old_scope='All' if scope=='Core' else scope.split(':',1)[1]
    q=prevtab[(prevtab.dataset=='GSE314943')&(prevtab.scope==old_scope)&(prevtab['filter']==filt)&(prevtab['sample']==ss)]
    assert len(q)==1
    observed=float(v.contribution.sum());expected=float(q.iloc[0]['score_'+mod])
    errors.append(dict(scope=scope,filter=filt,sample=ss,module=mod,absolute_error=abs(observed-expected)))
assert max(z['absolute_error'] for z in errors)<1e-8,'Score decomposition mismatch'
wr(pd.DataFrame(errors),'decomposition_audit.tsv')
# Correlation leverage in the six mouse samples: direction and range after deletion.
from scipy.stats import rankdata

def residual_corr(a,b,group):
    a=np.asarray(a,float);b=np.asarray(b,float);group=np.asarray(group,float)
    ok=np.isfinite(a)&np.isfinite(b)&np.isfinite(group);a=a[ok];b=b[ok];group=group[ok]
    if len(a)<5 or np.ptp(a)==0 or np.ptp(b)==0 or len(set(group))!=2:return np.nan
    d=np.column_stack([np.ones(len(a)),group]);a=rankdata(a);b=rankdata(b);a-=d@np.linalg.lstsq(d,a,rcond=None)[0];b-=d@np.linalg.lstsq(d,b,rcond=None)[0]
    if np.std(a)<1e-10 or np.std(b)<1e-10:return np.nan
    return float(np.clip(np.corrcoef(a,b)[0,1],-1,1))
leverage=[]
for (scope,filt),v in prevtab[(prevtab.dataset=='GSE314943')&prevtab.eligible20].groupby(['scope','filter'],observed=True):
    for g in TARGETS:
        for mod in SETS:
            z=v[['sample','group','expr_'+g,'score_'+mod]].dropna()
            if len(z)!=6:continue
            full=residual_corr(z['expr_'+g],z['score_'+mod],(z.group=='Silica').astype(int))
            for ss in z['sample']:
                q=z[z['sample']!=ss];val=residual_corr(q['expr_'+g],q['score_'+mod],(q.group=='Silica').astype(int))
                leverage.append(dict(scope=scope,filter=filt,gene=g,module=mod,removed_sample=ss,full_rho=full,loo_rho=val,same_sign=bool(np.sign(val)==np.sign(full)) if np.isfinite(val) and np.isfinite(full) and full!=0 else np.nan))
wr(pd.DataFrame(leverage),'mouse_correlation_leverage.tsv')
# Audit known NicheNet cache files without network or package changes.
cache=[]
for directory in [root/'00_raw',root/'02_analysis',root/'08_MRE11A',root/'resources']:
    if directory.exists():
        for dp,dirs,files in os.walk(directory):
            dirs[:]=[d for d in dirs if d not in ['comm_input','node_modules','.git']]
            for name in files:
                if name.endswith(('.rds','.RDS')) and any(k in name.lower() for k in ['ligand_target','lr_network','weighted_network','nichenet']):
                    p=Path(dp)/name;cache.append(dict(path=str(p),size_bytes=p.stat().st_size))
wr(pd.DataFrame(cache,columns=['path','size_bytes']),'nichenet_cache_inventory.tsv')
(out/'python_status.json').write_text(json.dumps({'status':'PASS','human_cells':len(hm),'human_subjects':hm['sample'].nunique(),'mouse_cells':len(mm),'max_decomposition_error':max(z['absolute_error'] for z in errors)},indent=2))
print('ANALYSIS_LOGIC_PASS=HUMAN_SCOPE_SOURCE_COUNTS_AND_GENE_DECOMPOSITION',flush=True)

PY_ANALYZE
cat > "$OUT/models.R" <<'R_MODELS'
out <- commandArgs(TRUE)[1]
t <- read.delim(file.path(out,'tables/human_sample_scores.tsv'),check.names=FALSE)
wr <- function(x,name) write.table(x,file.path(out,'tables',name),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
metrics <- grep('^(score_|expr_|pct_)',colnames(t),value=TRUE)
targets <- c('IL1R2','PENK','TTC39C','FOXP1','TNFRSF9','LAG3','MRE11')
mods <- sub('^score_','',grep('^score_',colnames(t),value=TRUE))
results <- list();assocs <- list();idx <- 0L;jdx <- 0L
for(scope in unique(t$scope)) for(filt in c('All','BG1')) {
  z <- t[t$scope==scope & t$filter==filt & t$cells>=20,,drop=FALSE]
  for(arm in c('Source','SourceDep','WithinG227')) {
    v <- if(arm=='WithinG227') z[z$series=='G227',,drop=FALSE] else z
    v$group <- factor(v$group,levels=c('Control','CWP'));v$series <- factor(v$series)
    for(metric in metrics) {
      w <- v[is.finite(v[[metric]]) & is.finite(v$library),,drop=FALSE]
      w <- droplevels(w);nc <- sum(w$group=='Control');nd <- sum(w$group=='CWP')
      row <- data.frame(scope=scope,filter=filt,model=arm,metric=metric,n_control=nc,n_case=nd,
        coefficient=NA_real_,se=NA_real_,CI_low=NA_real_,CI_high=NA_real_,p_nominal=NA_real_,residual_df=NA_real_,status='INSUFFICIENT_SAMPLES',stringsAsFactors=FALSE)
      if(nc>=2 && nd>=2) {
        terms <- c('group',if(arm!='WithinG227' && nlevels(w$series)>1) 'series',if(arm=='SourceDep') 'log1p(library)')
        form <- reformulate(terms,response=metric);X <- model.matrix(form,w)
        if(qr(X)$rank!=ncol(X) || nrow(X)-ncol(X)<2) row$status <- 'NON_IDENTIFIABLE' else {
          fit <- lm(form,data=w);co <- summary(fit)$coefficients
          row$status <- 'FIT'
          if('groupCWP'%in%rownames(co)) {
            b <- co['groupCWP',];df <- df.residual(fit);delta <- qt(.975,df)*b[2]
            row$coefficient <- b[1];row$se <- b[2];row$CI_low <- b[1]-delta;row$CI_high <- b[1]+delta;row$p_nominal <- b[4];row$residual_df <- df
          } else row$status <- 'GROUP_COEFFICIENT_MISSING'
        }
      }
      idx <- idx+1L;results[[idx]] <- row
    }
  }
  # Sample-level partial rank correlations adjust condition, source and mean library depth.
  for(g in targets) for(mod in mods) {
    ex <- paste0('expr_',g);sc <- paste0('score_',mod)
    v <- z[is.finite(z[[ex]]) & is.finite(z[[sc]]) & is.finite(z$library),,drop=FALSE]
    row <- data.frame(scope=scope,filter=filt,gene=g,module=mod,n=nrow(v),positive_samples=sum(v[[ex]]>0),rho=NA_real_,p_nominal=NA_real_,residual_df=NA_real_,status='INSUFFICIENT',stringsAsFactors=FALSE)
    if(nrow(v)>=10 && sum(v[[ex]]>0)>=4 && length(unique(v$group))==2 && sd(v[[ex]])>0 && sd(v[[sc]])>0) {
      v$group <- factor(v$group);v$series <- factor(v$series);v$depth_rank <- rank(v$library)
      terms <- c('group',if(nlevels(v$series)>1) 'series','depth_rank')
      X <- model.matrix(reformulate(terms),v);rankX <- qr(X)$rank;df <- nrow(v)-rankX-1L
      if(rankX==ncol(X) && df>=3) {
        a <- qr.resid(qr(X),rank(v[[ex]]));b <- qr.resid(qr(X),rank(v[[sc]]))
        if(sd(a)>1e-10 && sd(b)>1e-10) {
          r <- max(-1,min(1,cor(a,b)));row$rho <- r;row$residual_df <- df
          row$p_nominal <- if(abs(r)>=1) 0 else 2*pt(-abs(r)*sqrt(df/(1-r*r)),df)
          row$status <- 'EXPLORATORY_PARTIAL_RANK'
        }
      } else row$status <- 'NON_IDENTIFIABLE'
    }
    jdx <- jdx+1L;assocs[[jdx]] <- row
  }
}
r <- do.call(rbind,results);r$family <- sub('_.*$','',r$metric);r$FDR <- NA_real_
r$FDR_status <- 'NOT_COMPUTED_GROUP_N_LT10'
for(ix in split(seq_len(nrow(r)),interaction(r$scope,r$filter,r$model,r$family,drop=TRUE))) {
  keep <- ix[r$n_control[ix]>=10 & r$n_case[ix]>=10 & is.finite(r$p_nominal[ix])]
  if(length(keep)) {r$FDR[keep] <- p.adjust(r$p_nominal[keep],'BH');r$FDR_status[keep] <- 'BH_WITHIN_SCOPE_FILTER_MODEL_FAMILY'}
}
wr(r,'human_source_models.tsv');a <- do.call(rbind,assocs);a$FDR <- NA_real_
for(ix in split(seq_len(nrow(a)),interaction(a$scope,a$filter,drop=TRUE))) {
  keep <- ix[is.finite(a$p_nominal[ix])]
  if(length(keep)) a$FDR[keep] <- p.adjust(a$p_nominal[keep],'BH')
}
wr(a,'human_partial_correlations.tsv')
# Fixed-composition source-adjusted effects, explicitly retaining low-count subtype flags.
s <- read.delim(file.path(out,'tables/human_fixed_composition.tsv'),check.names=FALSE)
fres <- list();j <- 0L
for(ix in split(seq_len(nrow(s)),interaction(s$scope,s$filter,s$metric,drop=TRUE))) {
  v <- s[ix,,drop=FALSE];v <- v[is.finite(v$value),,drop=FALSE]
  if(!nrow(v))next
  row <- data.frame(scope=v$scope[1],filter=v$filter[1],metric=v$metric[1],n_control=sum(v$group=='Control'),n_case=sum(v$group=='CWP'),all_samples_subtypes_ge20=all(v$all_subtypes_eligible20),coefficient=NA_real_,status='INSUFFICIENT',stringsAsFactors=FALSE)
  if(row$n_control>=2 && row$n_case>=2) {
    v$group <- factor(v$group,levels=c('Control','CWP'));v$series <- factor(v$series)
    f <- reformulate(c('group',if(nlevels(v$series)>1)'series'),response='value');X <- model.matrix(f,v)
    if(qr(X)$rank==ncol(X) && nrow(X)>ncol(X)) {row$coefficient <- unname(coef(lm(f,v))['groupCWP']);row$status <- 'DESCRIPTIVE_FIXED_COMPOSITION'}
  }
  j <- j+1L;fres[[j]] <- row
}
if(length(fres))wr(do.call(rbind,fres),'human_fixed_composition_effects.tsv')
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
cat('HUMAN_SOURCE_MODELS_PASS\n')

R_MODELS
cat > "$OUT/finish.py" <<'PY_FINISH'
import os,json,tarfile
from pathlib import Path
import pandas as pd
import numpy as np
from shared import effect
r=Path(os.environ['ROOT']);o=Path(os.environ['OUT']);prev=Path(os.environ['PREV'])
def wr(d,n):d.to_csv(o/'tables'/n,sep='\t',index=False,na_rep='NA')
# Mouse module membership effects and leave-one-gene-out direction summaries.
c=pd.read_csv(o/'tables/module_gene_contributions.tsv',sep='\t');l=pd.read_csv(o/'tables/module_leave_one_gene.tsv',sep='\t')
rows=[]
for (scope,filt,mod,g),v in c[(c.dataset=='GSE314943')&c.eligible20].groupby(['scope','filter','module','gene'],observed=True):
    d=effect(v.loc[v.group=='Saline','contribution'],v.loc[v.group=='Silica','contribution'])
    rows.append(dict(scope=scope,filter=filt,module=mod,gene=g,**d))
wr(pd.DataFrame(rows),'mouse_gene_effects.tsv')
rows=[]
for (scope,filt,mod,g),v in l[(l.dataset=='GSE314943')&l.eligible20].groupby(['scope','filter','module','removed_gene'],observed=True):
    a=effect(v.loc[v.group=='Saline','full_score'],v.loc[v.group=='Silica','full_score'])
    b=effect(v.loc[v.group=='Saline','reduced_score'],v.loc[v.group=='Silica','reduced_score'])
    rows.append(dict(scope=scope,filter=filt,module=mod,removed_gene=g,baseline_difference=a['difference'],same_direction=np.sign(a['difference'])==np.sign(b['difference']),**b))
wr(pd.DataFrame(rows),'mouse_leave_one_gene_effects.tsv')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
h=pd.read_csv(o/'tables/human_source_models.tsv',sep='\t');q=h[(h.scope=='Core')&(h['filter']=='All')&(h.model=='Source')&h.metric.str.startswith('score_')].copy()
assert (q.status=='FIT').any(),'No estimable human Core source-adjusted module models'
q['module']=q.metric.str.replace('score_','',regex=False)
fig,ax=plt.subplots(figsize=(6,5));q=q.iloc[::-1]
for i,z in enumerate(q.itertuples()):
    if np.isfinite(z.coefficient):ax.errorbar(z.coefficient,i,xerr=[[z.coefficient-z.CI_low],[z.CI_high-z.coefficient]],fmt='o',color='#176B9B',capsize=2)
ax.axvline(0,color='grey',lw=.8);ax.set(yticks=range(len(q)),yticklabels=q.module.tolist(),xlabel='Adj effect');ax.spines[['top','right']].set_visible(False);fig.tight_layout();fig.savefig(o/'figures/human_modules.png',dpi=300);plt.close(fig)
(o/'captions.txt').write_text('human_modules.png: y=module abbreviations inherited from step1 module_definition.tsv; x=source-adjusted CWP-Control mean-score coefficient; bars=model-based95% confidence intervals, each subject weighted equally, Core scope, >=20 cells, All filter. Models use group+series; confidence intervals describe these exploratory models. No significance stars. Adj effect=adjusted effect.\nCaspGsd=CASP1 and GSDMD co-detected; Pyro3=CASP1/GSDMD/PYCARD co-detected; AltPyro=CASP3/GSDME co-detected. These are transcript co-detection prerequisites, not protein cleavage measurements.\n')
(o/'methods.txt').write_text('Human lung uses the original final Seurat object and subject/source mapping from existing project scripts. Primary=all final CD4; Core excludes CD4state. No external sample added. The source/condition design and subtype cell counts are exported.\nSame12 frozen modules and matched-expression scoring formula as step1, human controls selected within human cohort then fixed across scopes/filters. Targets excluded. Human and mouse absolute scores are not pooled. Full raw library denominators are retained after symbol aggregation.\nPer-subject mean scores/target CP10k/detection. Source model=response~group+series; SourceDep additionally controls log1p(mean library); WithinG227 compares within the reused mixed-condition source. Minimum20 cells. Non-identifiable or insufficient comparisons explicitly flagged. Original edgeR DEG is not refitted here. Linear models of target expression are exploratory mean-CP10k models, not count-model DEG.\nGroup model BH correction is calculated only when each group has>=10 subjects; all nominal values retained for audit. Partial rank associations require>=10 total subjects,>=4 target-positive subjects and both conditions; condition, source and depth controlled; BH across all target-module associations within scope/filter. Small-sample and source-confounding limitations apply.\nFixed-composition weights are original equal-subject proportions and remain fixed across filters. Missing subtype scores produce NA; low-count subtypes flagged. BG1 uses the same predefined panel fraction<=1% rule.\nGene contributions exactly sum to module mean minus fixed control mean. Leave-one-gene-out keeps original control genes; no new significance testing. Mouse correlation leave-one-sample-out is descriptive with five remaining samples; original six-sample correlations are not selected for significance.\nNicheNet package/cache availability recorded only; ligand-target inference follows verified receiver program definition. No network inference claimed in this stage.\n')
required=['human_source_models.tsv','human_partial_correlations.tsv','module_gene_contributions.tsv','pyro_codetection.tsv','mouse_leave_one_gene_effects.tsv']
assert all((o/'tables'/p).is_file() for p in required)
status=json.loads((o/'python_status.json').read_text());status.update(stage='mechanism_step2_human_and_driver_audit',previous=str(prev),output=str(o),pending=['ligand_target_inference','candidate_regulatory_network'])
(o/'status.json').write_text(json.dumps(status,indent=2))
dest=r/'tmp/lung_mechanism_step2_return.tar.gz';tmp=dest.with_name(dest.name+'.part')
with tarfile.open(tmp,'w:gz') as t:
    for folder in ['tables','figures']:
        for p in sorted((o/folder).iterdir()):t.add(p,arcname=str(p.relative_to(o)))
    for p in sorted(o.iterdir()):
        if p.is_file() and p.suffix in ['.R','.py','.txt','.json','.sh']:t.add(p,arcname=p.name)
assert tmp.stat().st_size<25*1024*1024,'Package exceeds25MiB; all results retained on server'
os.replace(tmp,dest)
(r/'tmp/lung_mechanism_step2_latest.txt').write_text(str(o)+'\n')
print('RETURN_PACKAGE='+str(dest));print('RETURN_BYTES='+str(dest.stat().st_size));print('LUNG_MECHANISM_STEP2_PASS',flush=True)

PY_FINISH
conda run --no-capture-output -n batch_bench python -m py_compile "$OUT/shared.py" "$OUT/analyze.py" "$OUT/finish.py"
conda run --no-capture-output -n R_4.5.2 Rscript -e 'for(f in commandArgs(TRUE)) invisible(parse(file=f)); cat("R_PARSE_PASS\n")' "$OUT/export.R" "$OUT/models.R"
echo 'SYNTAX_CHECK=PASS'
echo "OUTPUT=$OUT"
for f in status.json comm_input/CD4_counts.mtx.gz comm_input/CD4_genes.tsv comm_input/CD4_cells.tsv tables/score_controls.tsv.gz tables/sample_scores.tsv; do
  test -s "$PREV/$f" || { echo "MISSING_STEP1_INPUT=$PREV/$f"; exit 1; }
done
conda run --no-capture-output -n batch_bench python -c 'import json,os; from pathlib import Path; s=json.loads((Path(os.environ["PREV"])/"status.json").read_text()); assert s["status"]=="PASS" and s["stage"]=="mechanism_step1"; print("PREVIOUS_STAGE_CHECK=PASS")'
conda run --no-capture-output -n R_4.5.2 Rscript "$OUT/export.R" "$ROOT" "$OUT"
conda run --no-capture-output -n batch_bench python -u "$OUT/analyze.py"
conda run --no-capture-output -n R_4.5.2 Rscript "$OUT/models.R" "$OUT"
conda run --no-capture-output -n batch_bench python -u "$OUT/finish.py"
