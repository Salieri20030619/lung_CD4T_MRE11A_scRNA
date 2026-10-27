#!/usr/bin/env bash
# Step3: replicate-aware mouse NicheNet + human candidate correlation sensitivity.
# Run: bash lung_mechanism_step3.sh
# Optional: PROXY_URL=http://127.0.0.1:<active_forwarded_port> bash lung_mechanism_step3.sh
set -Eeuo pipefail
ROOT="${ROOT:-/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA}"
PREV="${PREV:-$ROOT/02_analysis/mechanism_step1_20260927_152218_s4V5eF}"
PREV2="${PREV2:-$ROOT/02_analysis/mechanism_step2_20260927_154919_EsN4gk}"
CACHE="${CACHE:-$ROOT/00_raw/nichenet_v2_mouse}"
export ROOT PREV PREV2 CACHE
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-2}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
export MPLBACKEND=Agg
bash -n "$0"
mkdir -p "$ROOT/tmp" "$ROOT/logs" "$ROOT/02_analysis"
exec 9>"$ROOT/tmp/lung_mechanism_step3.lock"
flock -n 9 || { echo MECHANISM_STEP3_ALREADY_RUNNING; exit 1; }
OUT=$(mktemp -d "$ROOT/02_analysis/mechanism_step3_$(date +%Y%m%d_%H%M%S)_XXXXXX")
export OUT
trap 'rc=$?; if ((rc != 0)); then printf "LUNG_MECHANISM_STEP3_FAIL EXIT_CODE=%s OUTPUT=%s\n" "$rc" "$OUT"; fi' EXIT
mkdir -p "$OUT/tables" "$OUT/figures" "$OUT/objects"
cp -- "$0" "$OUT/run.sh"
printf 'OUTPUT=%s\n' "$OUT"
cat > "$OUT/fetch.py" <<'PY_FETCH'
import os, json, hashlib, subprocess, time
from pathlib import Path
out=Path(os.environ['OUT']); cache=Path(os.environ['CACHE']); cache.mkdir(parents=True,exist_ok=True)
prev=Path(os.environ['PREV']); prev2=Path(os.environ['PREV2'])
for p,stage in [(prev,'mechanism_step1'),(prev2,'mechanism_step2_human_and_driver_audit')]:
    s=json.loads((p/'status.json').read_text()); assert s['status']=='PASS' and s['stage']==stage,(p,s)
for line in ['CD4','Mac']:
    for suf in ['counts.mtx.gz','genes.tsv','cells.tsv']:
        assert (prev/'comm_input'/f'{line}_{suf}').is_file()
base=['curl','-L','--fail','--show-error','--silent','--connect-timeout','30','--speed-time','300','--speed-limit','128']
if os.environ.get('PROXY_URL'):base+=['--proxy',os.environ['PROXY_URL'],'--noproxy','']
def get(url,path,resume=True):
    for attempt in range(1,7):
        args=base+(['--continue-at','-'] if resume else [])+['-o',str(path),url]
        rc=subprocess.run(args).returncode
        if rc==0:return
        print(f'DOWNLOAD_RETRY file={path.name} attempt={attempt} curl_exit={rc}',flush=True)
        time.sleep(10)
    raise RuntimeError('Download failed: '+path.name+'; check server proxy/tunnel, then rerun.')
meta=cache/'zenodo_7074291.json'
if not meta.exists():
    tmp=meta.with_suffix('.json.part');get('https://zenodo.org/api/records/7074291',tmp,False)
    j=json.loads(tmp.read_text());assert 'files' in j;tmp.replace(meta)
j=json.loads(meta.read_text()); records={v['key']:v for v in j['files']}
rows=[]
for name in ['lr_network_mouse_21122021.rds','ligand_target_matrix_nsga2r_final_mouse.rds']:
    spec=records[name]; expected=spec['checksum'];algorithm,digest=expected.split(':',1)
    path=cache/name
    def verify(p):
        h=hashlib.new(algorithm)
        with p.open('rb') as f:
            for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
        return h.hexdigest()==digest and p.stat().st_size==int(spec['size'])
    if path.exists():
        if not verify(path):raise RuntimeError('Cached file failed checksum: '+str(path))
    else:
        part=path.with_suffix('.rds.part')
        print('DOWNLOAD_START='+name,flush=True)
        get('https://zenodo.org/records/7074291/files/'+name+'?download=1',part)
        if not verify(part):raise RuntimeError('Downloaded file failed checksum: '+str(part))
        part.replace(path)
    rows.append(dict(file=name,bytes=path.stat().st_size,checksum=expected,status='PASS'))
    print('RESOURCE_PASS='+name,flush=True)
(out/'resource_manifest.json').write_text(json.dumps(rows,indent=2))
print('NICHE_NET_RESOURCE_CHECK=PASS',flush=True)

PY_FETCH
cat > "$OUT/human_audit.R" <<'R_HUMAN'
args <- commandArgs(TRUE);prev2 <- args[1];out <- args[2]
rd <- function(f)read.delim(f,check.names=FALSE,stringsAsFactors=FALSE)
wr <- function(x,n)write.table(x,file.path(out,'tables',n),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
t <- rd(file.path(prev2,'tables/human_sample_scores.tsv'))
# Correct the previous reporting gate: BH is mathematically valid with fewer than10 samples/group.
m <- rd(file.path(prev2,'tables/human_source_models.tsv'));m$FDR <- NA_real_
for(ix in split(seq_len(nrow(m)),interaction(m$scope,m$filter,m$model,m$family,drop=TRUE))) {
 j <- ix[is.finite(m$p_nominal[ix])];if(length(j))m$FDR[j] <- p.adjust(m$p_nominal[j],'BH')
}
m$FDR_status <- 'BH_ALL_FINITE_WITHIN_SCOPE_FILTER_MODEL_FAMILY'
wr(m,'human_source_models_BH.tsv')
calc <- function(v,g,mod,adjust) {
 base <- c('group',if(length(unique(v$series))>1)'series','depth_rank')
 if(!is.null(adjust))base <- c(base,adjust)
 v$depth_rank <- rank(v$library);v$group <- factor(v$group);v$series <- factor(v$series)
 X <- model.matrix(reformulate(base),v);df <- nrow(v)-qr(X)$rank-1
 if(qr(X)$rank<ncol(X)||df<3)return(c(rho=NA,p=NA,df=df))
 u <- qr.resid(qr(X),rank(v[[paste0('expr_',g)]]));w <- qr.resid(qr(X),rank(v[[paste0('score_',mod)]]))
 if(sd(u)<1e-10||sd(w)<1e-10)return(c(rho=NA,p=NA,df=df))
 r <- max(-1,min(1,cor(u,w)));p <- if(abs(r)>=1)0 else 2*pt(-abs(r)*sqrt(df/(1-r*r)),df)
 c(rho=r,p=p,df=df)
}
res <- list();loo <- list()
for(filt in c('All','BG1')) {
 v <- t[t$scope=='Core' & t$filter==filt & t$cells>=20,,drop=FALSE]
 # Use all observed cells in each Core subtype, including <20-cell subtypes, to compute proportions.
 sub <- t[grepl('^Subtype:',t$scope)&t$filter==filt & t$scope!='Subtype:CD4state',,drop=FALSE]
 for(label in unique(sub$scope)) {
   q <- sub[sub$scope==label,,drop=FALSE];cv <- paste0('frac_',make.names(sub('^Subtype:','',label)))
   v[[cv]] <- q$cells[match(v$sample,q$sample)];v[[cv]][is.na(v[[cv]])] <- 0;v[[cv]] <- v[[cv]]/v$cells
 }
 vars <- grep('^frac_',names(v),value=TRUE)
 stopifnot(max(abs(rowSums(v[,vars,drop=FALSE])-1))<1e-8)
 for(g in c('IL1R2','PENK','TTC39C','FOXP1','TNFRSF9','LAG3','MRE11'))for(mod in sub('^score_','',grep('^score_',names(t),value=TRUE)))for(adj in c('Base','JointComp',vars)) {
  w <- v[is.finite(v[[paste0('expr_',g)]]) & is.finite(v[[paste0('score_',mod)]]),,drop=FALSE]
  rr <- calc(w,g,mod,if(adj=='Base')NULL else if(adj=='JointComp')head(vars,-1) else adj)
  res[[length(res)+1]] <- data.frame(filter=filt,gene=g,module=mod,adjust=adj,n=nrow(w),rho=rr['rho'],p_nominal=rr['p'],df=rr['df'])
  if(adj=='Base')for(i in seq_len(nrow(w))) {
    r <- calc(w[-i,,drop=FALSE],g,mod,NULL)
    loo[[length(loo)+1]] <- data.frame(filter=filt,gene=g,module=mod,removed_sample=w$sample[i],rho=r['rho'],full_rho=rr['rho'])
  }
 }
}
r <- do.call(rbind,res);r$FDR <- NA_real_
for(ix in split(seq_len(nrow(r)),interaction(r$filter,r$adjust,drop=TRUE))) {
 j <- ix[is.finite(r$p_nominal[ix])];if(length(j))r$FDR[j] <- p.adjust(r$p_nominal[j],'BH')
}
wr(r,'human_composition_correlations.tsv');wr(do.call(rbind,loo),'human_correlation_LOO.tsv')
cat('HUMAN_CORRELATION_AUDIT_PASS\n')

R_HUMAN
cat > "$OUT/pairing.R" <<'R_PAIR'
# Sourced after analyze.R: uses sc/tl generated above, independent-library descriptive association.
partial_group <- function(x,y,grp) {
 ok <- is.finite(x)&is.finite(y);x <- x[ok];y <- y[ok];grp <- grp[ok]
 if(length(x)<5||length(unique(grp))!=2||min(table(grp))<2||sd(x)==0||sd(y)==0)return(NA_real_)
 X <- model.matrix(~factor(grp));u <- qr.resid(qr(X),rank(x));v <- qr.resid(qr(X),rank(y))
 if(sd(u)<1e-10||sd(v)<1e-10)return(NA_real_)
 cor(u,v)
}
pr <- list()
if(nrow(tl)) {
 q <- unique(tl[tl$test=='FDR05' & tl$ligand_rank<=20,c('receiver','filter','ligand','target','canonical')])
 for(i in seq_len(nrow(q))) {
   rz <- sc[[paste('CD4',q$receiver[i],q$filter[i],sep='_')]]
   if(!q$target[i]%in%rownames(rz$avg))next
   for(sk in names(sc)[vapply(sc,function(z)z$lineage=='Mac',logical(1))]) {
     sz <- sc[[sk]];lg <- q$ligand[i]
     if(!lg%in%rownames(sz$avg))next
     keep <- rz$sm$eligible20&sz$sm$eligible20
     x <- sz$avg[lg,keep];y <- rz$avg[q$target[i],keep];gr <- rz$sm$group[keep]
     r <- if(sum(keep)>=6)partial_group(x,y,gr) else NA_real_
     lr <- if(sum(keep)>=6)vapply(seq_along(x),function(j)partial_group(x[-j],y[-j],gr[-j]),numeric(1)) else NA_real_
     pr[[length(pr)+1]] <- data.frame(receiver=q$receiver[i],filter=q$filter[i],sender=sz$scope,ligand=lg,target=q$target[i],canonical=q$canonical[i],n=sum(keep),n_control=sum(gr=='Saline'),n_case=sum(gr=='Silica'),rho_group_adjusted=r,loo_min=if(any(is.finite(lr)))min(lr,na.rm=TRUE) else NA_real_,loo_max=if(any(is.finite(lr)))max(lr,na.rm=TRUE) else NA_real_,loo_same_fraction=if(is.finite(r)&&r!=0&&any(is.finite(lr)))mean(sign(lr[is.finite(lr)])==sign(r)) else NA_real_,evidence='DESCRIPTIVE_NOT_CAUSAL')
   }
 }
}
save_rows(pr,'paired_ligand_target_correlations.tsv',c('receiver','filter','sender','ligand','target','n','rho_group_adjusted'))
cat('MATCHED_SAMPLE_CHECK_PASS\n')

R_PAIR
cat > "$OUT/analyze.R" <<'R_NICHE'
args <- commandArgs(TRUE)
prev <- args[1]; prev2 <- args[2]; out <- args[3]; cache <- args[4]
suppressPackageStartupMessages({library(Matrix);library(edgeR);library(nichenetr)})
set.seed(314943)
wr <- function(d,n) write.table(d,file.path(out,'tables',n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
rd <- function(p) read.delim(p,check.names=FALSE,stringsAsFactors=FALSE)
canon <- function(x) {
 x <- toupper(trimws(x)); a <- c(MRE11A='MRE11',ATP5A1='ATP5F1A',ATP5B='ATP5F1B',ATP5C1='ATP5F1C',ATP5G1='ATP5MC1',ATP5O='ATP5PO',PARK2='PRKN')
 i <- x%in%names(a);x[i] <- unname(a[x[i]]);x
}
targets <- c('IL1R2','PENK','TTC39C','FOXP1','TNFRSF9','LAG3','MRE11')
# Fixed background panel identical to step1. CD4 only; do not remove macrophage genes from macrophages.
bg_panel <- strsplit('SFTPA1 SFTPA2 SFTPB SFTPC SFTPD SCGB1A1 SCGB3A1 SCGB3A2 ABCA3 SLC34A2 NAPSA LYZ LYZ1 LYZ2 S100A8 S100A9 CSF1R LST1 TYROBP FCER1G',' ')[[1]]
load_input <- function(lineage,expected) {
 pre <- file.path(prev,'comm_input',lineage)
 con <- gzfile(paste0(pre,'_counts.mtx.gz'),'rt');x <- as(readMM(con),'CsparseMatrix');close(con)
 g <- rd(paste0(pre,'_genes.tsv'));m <- rd(paste0(pre,'_cells.tsv'))
 stopifnot(ncol(x)==expected,ncol(x)==nrow(m),nrow(x)==nrow(g),!anyDuplicated(m$cell),!anyNA(m),all(m$library>0),all(x@x>=0),all(is.finite(x@x)),all(abs(x@x-round(x@x))<1e-6))
 stopifnot(setequal(unique(m$group),c('Saline','Silica')),length(unique(m$sample))==6,all(m$group==sub('_.*','',m$sample)))
 stopifnot(all(Matrix::colSums(x)<=m$library+1e-6))
 b <- g$canonical%in%bg_panel | grepl('^(IGH[VDJMAEG]|IGK[VJC]|IGL[VJC])',g$canonical)
 m$bg <- as.numeric(Matrix::colSums(x[b,,drop=FALSE]))/m$library
 sym <- as.character(nichenetr::convert_alias_to_symbols(g$mouse_symbol,'mouse',verbose=FALSE))
 stopifnot(length(sym)==nrow(g));bad <- is.na(sym)|!nzchar(sym);sym[bad] <- g$mouse_symbol[bad]
 u <- sort(unique(sym));pr <- sparseMatrix(i=match(sym,u),j=seq_along(sym),x=1,dims=c(length(u),length(sym)))
 x <- pr%*%x;rownames(x) <- u;colnames(x) <- m$cell
 wr(data.frame(lineage,canonical=g$canonical,original=g$mouse_symbol,mapped=sym),paste0(lineage,'_gene_mapping.tsv'))
 list(x=x,m=m)
}
a <- load_input('CD4',2311);b <- load_input('Mac',1796)
lr <- as.data.frame(readRDS(file.path(cache,'lr_network_mouse_21122021.rds')))
stopifnot(all(c('from','to')%in%names(lr)))
lr <- unique(lr[,c('from','to')]);lr$from <- as.character(lr$from);lr$to <- as.character(lr$to)
lt <- readRDS(file.path(cache,'ligand_target_matrix_nsga2r_final_mouse.rds'))
stopifnot(length(dim(lt))==2,!is.null(rownames(lt)),!is.null(colnames(lt)),!anyDuplicated(rownames(lt)),!anyDuplicated(colnames(lt)))
stopifnot(sum(rownames(a$x)%in%rownames(lt))>3000,sum(rownames(b$x)%in%unique(lr$from))>50)
wr(data.frame(resource=c('target_rows','ligand_columns','LR_pairs'),n=c(nrow(lt),ncol(lt),nrow(lr))),'prior_dimensions.tsv')
# All six samples remain in descriptive matrices; only >=20 cells enter inference/expression support.
samples <- sort(unique(a$m$sample)); stopifnot(setequal(samples,unique(b$m$sample)))
make_scope <- function(obj,scope,filt,lineage) {
 m <- obj$m;x <- obj$x
 keep <- if(scope%in%c('Core','AllMac')) rep(TRUE,nrow(m)) else m$subtype==scope
 if(filt=='BG1')keep <- keep&m$bg<=.01
 x <- x[,keep,drop=FALSE];m <- m[keep,,drop=FALSE]
 agg <- sparseMatrix(i=seq_len(nrow(m)),j=match(m$sample,samples),x=1,dims=c(nrow(m),length(samples)))
 pb <- as.matrix(x%*%agg);det <- as.matrix((x>0)%*%agg)
 n <- as.numeric(table(factor(m$sample,levels=samples)));lib <- as.numeric(tapply(m$library,factor(m$sample,levels=samples),sum));lib[is.na(lib)] <- 0
 pct <- sweep(det,2,pmax(n,1),'/');pct[,n==0] <- NA_real_
 # Per-sample mean per-cell CP10k for descriptive comparisons, not the DE input.
 norm <- x%*%Diagonal(x=1e4/m$library);avg <- sweep(as.matrix(norm%*%agg),2,pmax(n,1),'/');avg[,n==0] <- NA_real_
 colnames(pb) <- colnames(det) <- colnames(pct) <- colnames(avg) <- samples
 sm <- data.frame(sample=samples,group=sub('_.*','',samples),cells=n,library=lib,eligible20=n>=20,scope=scope,filter=filt,lineage=lineage)
 list(pb=pb,pct=pct,avg=avg,sm=sm,scope=scope,filter=filt,lineage=lineage)
}
sc <- list()
for(s in c('Core',sort(unique(a$m$subtype)))) for(f in c('All','BG1')) sc[[paste('CD4',s,f,sep='_')]] <- make_scope(a,s,f,'CD4')
for(s in c('AllMac',sort(unique(b$m$subtype)))) sc[[paste('Mac',s,'All',sep='_')]] <- make_scope(b,s,'All','Mac')
wr(do.call(rbind,lapply(sc,`[[`,'sm')),'scope_sample_counts.tsv')
# Independent-library edgeR; primary inference requires three samples per group.
fits <- list();degs <- list();de_status <- list()
for(key in names(sc)) {
 z <- sc[[key]];sel <- z$sm$eligible20;md <- z$sm[sel,,drop=FALSE]
 nc <- sum(md$group=='Saline');nd <- sum(md$group=='Silica')
 st <- data.frame(key,n_control=nc,n_case=nd,status='SKIP_NOT_3V3',tested=0)
 if(nc>=3 && nd>=3) {
   grp <- factor(md$group,levels=c('Saline','Silica'));design <- model.matrix(~grp)
   y <- DGEList(z$pb[,sel,drop=FALSE],lib.size=md$library,group=grp)
   keep <- filterByExpr(y,design=design)
   if(sum(keep)<500)stop('Insufficient genes in ',key)
   y <- y[keep,,keep.lib.sizes=TRUE];y <- calcNormFactors(y)
   y <- estimateDisp(y,design,robust=TRUE);fit <- glmQLFit(y,design,robust=TRUE)
   tt <- topTags(glmQLFTest(fit,coef=2),n=Inf,sort.by='none')$table
   d <- data.frame(gene=rownames(tt),tt,check.names=FALSE);d$canonical <- canon(d$gene)
   d$gene_role <- ifelse(d$canonical%in%targets,'Candidate','Other')
   wr(d,paste0('DE_',key,'.tsv'));degs[[key]] <- d
   fits[[key]] <- list(y=y,design=design,fit=fit,metadata=md)
   st$status <- 'PASS_3V3';st$tested <- nrow(d)
 }
 de_status[[key]] <- st
}
wr(do.call(rbind,de_status),'DE_status.tsv');saveRDS(fits,file.path(out,'objects/pseudobulk_models.rds'))
stopifnot(all(c('CD4_Core_All','CD4_ActCD4_All','CD4_Core_BG1','CD4_ActCD4_BG1')%in%names(degs)))
# Expression support: >=10% cells in >=2 qualifying disease samples (5% sensitivity is exported).
expressed <- function(z,group='Silica',cut=.1) {
 j <- which(z$sm$eligible20 & z$sm$group==group)
 if(length(j)<2)return(character())
 rownames(z$pct)[rowSums(z$pct[,j,drop=FALSE]>=cut,na.rm=TRUE)>=2]
}
expr_rows <- list()
for(key in names(sc)) {
 z <- sc[[key]];gs <- if(z$lineage=='Mac') union(lr$from,rownames(z$pb)[canon(rownames(z$pb))%in%targets]) else union(lr$to,rownames(z$pb)[canon(rownames(z$pb))%in%targets])
 gs <- intersect(gs,rownames(z$pb))
 for(j in seq_along(samples))expr_rows[[length(expr_rows)+1]] <- data.frame(key=key,lineage=z$lineage,scope=z$scope,filter=z$filter,gene=gs,canonical=canon(gs),sample=samples[j],group=z$sm$group[j],cells=z$sm$cells[j],eligible20=z$sm$eligible20[j],pct=100*z$pct[gs,j],mean_CP10k=z$avg[gs,j],raw_sum=z$pb[gs,j])
}
wr(do.call(rbind,expr_rows),'ligand_receptor_sample_expression.tsv')
# Sender support attaches subtype-specific replicate evidence; sparse control sources are never zero-imputed DE.
lr_rows <- list();activities <- list();set_status <- list();target_links <- list();network_edges <- list()
receivers <- c('Core','ActCD4');prior_audit <- list()
for(receiver in receivers)for(filt in c('All','BG1')) {
 key <- paste('CD4',receiver,filt,sep='_');z <- sc[[key]];d <- degs[[key]]
 # Background is tested genes AND sufficiently expressed in either condition, inside the prior model.
 bg <- Reduce(intersect,list(d$gene,union(expressed(z,'Saline'),expressed(z,'Silica')),rownames(lt)))
 receptor_sets <- list(Strict10=intersect(expressed(z,cut=.1),lr$to),Sens5=intersect(expressed(z,cut=.05),lr$to))
 for(cutlabel in names(receptor_sets)) {
   cut <- if(cutlabel=='Strict10').1 else .05; rec <- receptor_sets[[cutlabel]]
   for(sk in names(sc)[vapply(sc,function(v)v$lineage=='Mac',logical(1))]) {
     sender <- sc[[sk]]; lig <- intersect(expressed(sender,cut=cut),colnames(lt))
     edges <- lr[lr$from%in%lig & lr$to%in%rec,,drop=FALSE]
     if(!nrow(edges))next
     nj <- which(sender$sm$eligible20 & sender$sm$group=='Silica');rj <- which(z$sm$eligible20 & z$sm$group=='Silica')
     edges$receiver <- receiver;edges$filter <- filt;edges$threshold <- cutlabel;edges$sender <- sender$scope
     edges$ligand_case_samples <- rowSums(sender$pct[edges$from,nj,drop=FALSE]>=cut,na.rm=TRUE)
     edges$receptor_case_samples <- rowSums(z$pct[edges$to,rj,drop=FALSE]>=cut,na.rm=TRUE)
     edges$sender_control_eligible <- sum(sender$sm$eligible20&sender$sm$group=='Saline')
     edges$source_status <- ifelse(edges$sender_control_eligible>=3,'BALANCED_3V3','CASE_SOURCE_ONLY_NO_BALANCED_DE')
     lr_rows[[length(lr_rows)+1]] <- edges
   }
 }
 # Primary candidates from the union of supported macrophage subtypes, including disease-enriched Spp1Mac.
 edges <- do.call(rbind,lr_rows)
 ee <- if(is.null(edges))data.frame() else edges[edges$receiver==receiver & edges$filter==filt & edges$threshold=='Strict10',,drop=FALSE]
 potential <- if(nrow(ee)) sort(unique(ee$from)) else character()
 cs <- rownames(lt)[canon(rownames(lt))%in%targets]
 prior_audit[[key]] <- data.frame(receiver,filter=filt,canonical=targets,prior_present=targets%in%canon(rownames(lt)),expression_present=targets%in%canon(rownames(z$pb)),DE_tested=targets%in%d$canonical)
 # Gene-set gates are frozen; nominal lists remain labeled exploratory, never substituted for FDR sets.
 for(test in c('FDR05','Nominal05'))for(direction in c('Up','Down')) {
   pv <- if(test=='FDR05')d$FDR else d$PValue
   gs <- d$gene[is.finite(pv)&pv<.05 & (if(direction=='Up')d$logFC>=.5 else d$logFC<= -.5)]
   gs <- intersect(gs,bg)
   st <- data.frame(receiver,filter=filt,test,direction,background=length(bg),geneset=length(gs),ligands=length(potential),status='PASS')
   if(length(bg)<500)st$status <- 'SKIP_SMALL_BACKGROUND'
   if(length(gs)<10)st$status <- 'SKIP_GENESET_LT10'
   if(length(setdiff(bg,gs))<10)st$status <- 'SKIP_NEGATIVE_GENES_LT10'
   if(length(potential)<2)st$status <- 'SKIP_LIGANDS_LT2'
   set_status[[length(set_status)+1]] <- st
   wr(data.frame(gene=gs),paste0('geneset_',key,'_',test,'_',direction,'.tsv'))
   if(st$status!='PASS')next
   act <- as.data.frame(nichenetr::predict_ligand_activities(geneset=gs,background_expressed_genes=bg,ligand_target_matrix=lt,potential_ligands=potential))
   stopifnot(all(c('test_ligand','aupr_corrected','aupr')%in%names(act)),any(is.finite(act$aupr_corrected)))
   act <- act[order(-act$aupr_corrected,act$test_ligand),,drop=FALSE];act$rank <- rank(-act$aupr_corrected,ties.method='min')
   act$receiver <- receiver;act$filter <- filt;act$test <- test;act$direction <- direction;act$geneset_n <- length(gs);act$background_n <- length(bg)
   activities[[length(activities)+1]] <- act
   best <- head(act$test_ligand[is.finite(act$aupr_corrected)&act$aupr_corrected>0],20)
   # Full selected-candidate prior weights: no arbitrary top-target cutoff can hide absent support.
   for(lig in potential)for(g in cs) {
     weights <- lt[bg,lig];w <- as.numeric(lt[g,lig]);rr <- (1+sum(weights>w))/(1+length(weights))
     row <- data.frame(receiver,filter=filt,test,direction,ligand=lig,target=g,canonical=canon(g),weight=w,rank_fraction=rr,target_in_background=g%in%bg,target_in_geneset=g%in%gs,ligand_rank=act$rank[match(lig,act$test_ligand)],aupr_corrected=act$aupr_corrected[match(lig,act$test_ligand)])
     target_links[[length(target_links)+1]] <- row
   }
   # Supported ligand-receptor + target prior form hypotheses, not a proven receptor-to-target path.
   for(lig in best) {
     tg <- unique(c(cs,head(gs[order(-lt[gs,lig])],30)))
     lrsel <- ee[ee$from==lig,,drop=FALSE]
     for(g in tg) {
       row <- lrsel;row$target <- g;row$canonical <- canon(g);row$weight <- as.numeric(lt[g,lig]);row$test <- test;row$direction <- direction
       row$target_in_geneset <- g%in%gs;row$ligand_rank <- act$rank[match(lig,act$test_ligand)];row$aupr_corrected <- act$aupr_corrected[match(lig,act$test_ligand)]
       row$evidence <- 'LR_EXPRESSION_PLUS_UNSIGNED_LIGAND_TARGET_PRIOR_NOT_CAUSAL'
       network_edges[[length(network_edges)+1]] <- row
     }
   }
   # Source-specific ligand activity uses the same receiver gene set; source attribution is via expression.
   cat('NICHE_NET_DONE=',key,' ',test,' ',direction,' genes=',length(gs),' ligands=',length(potential),'\n',sep='')
 }
}
save_rows <- function(lst,name,columns) {
 v <- if(length(lst))do.call(rbind,lst) else as.data.frame(setNames(replicate(length(columns),character(),simplify=FALSE),columns))
 wr(v,name);v
}
ed <- save_rows(lr_rows,'supported_LR_edges.tsv',c('from','to','receiver','filter','threshold','sender'))
ac <- save_rows(activities,'ligand_activities.tsv',c('test_ligand','aupr_corrected','rank','receiver','filter','test','direction'))
tl <- save_rows(target_links,'candidate_ligand_prior.tsv',c('ligand','target','weight','receiver','filter','test','direction'))
ne <- save_rows(network_edges,'hypothesis_network.tsv',c('from','to','target','weight','receiver','filter','test','direction'))
wr(do.call(rbind,set_status),'geneset_status.tsv');wr(do.call(rbind,prior_audit),'candidate_prior_coverage.tsv')
# Joint IL18 receptor chains and selected multimeric receptor checks, per sample; missing chain = NA.
complexes <- list(Il18=c('Il18r1','Il18rap'),Il1b=c('Il1r1','Il1rap'),Il10=c('Il10ra','Il10rb'),Tgfb1=c('Tgfbr1','Tgfbr2'),Il6=c('Il6r','Il6st'),Il12=c('Il12rb1','Il12rb2'))
cr <- list()
for(receiver in receivers)for(filt in c('All','BG1'))for(cx in names(complexes)) {
 m <- a$m;take <- if(receiver=='Core')rep(TRUE,nrow(m)) else m$subtype==receiver
 if(filt=='BG1')take <- take&m$bg<=.01
 chain <- complexes[[cx]];present <- all(chain%in%rownames(a$x))
 for(ss in samples) {
   ii <- which(take&m$sample==ss);n <- length(ii)
   co <- if(present&&n>0)sum(Matrix::colSums(a$x[chain,ii,drop=FALSE]>0)==length(chain)) else NA_integer_
   cr[[length(cr)+1]] <- data.frame(receiver,filter=filt,ligand=cx,sample=ss,group=sub('_.*','',ss),cells=n,eligible20=n>=20,chains=paste(chain,collapse=';'),all_features_present=present,codetected=co,pct=if(n)100*co/n else NA_real_)
 }
}
wr(do.call(rbind,cr),'receptor_complex_codetection.tsv')
# Receiver gene modules are annotated on target links, but a positive prior has no activating/inhibitory sign.
writeLines(c('Primary: independent-library edgeR QL, 3 Saline versus 3 Silica, >=20 cells/library; full original library offsets and TMM. Gene-wise BH within each contrast, no sample-count restriction on BH.',
'All and BG1<=1% are separate. Macrophages are not subjected to the CD4 background filter. No cells used as independent replicates. Core and ActCD4 receive separate NicheNet models.',
'Expression-supported LR: >=10% cells in >=2 eligible Silica libraries for a sender subtype and receiver; >=5% is a separately labeled expression sensitivity. Spp1Mac with sparse controls has case-source support, not a valid balanced disease-DE claim.',
'Background = edgeR-tested AND >=10% expression in >=2 libraries in either condition AND present in the mouse prior. Genesets: FDR<0.05 and |logFC|>=0.5, up/down separately; nominal P<0.05 sets are explicitly exploratory. At least10 genes and500 background genes required. Never substitute a nominal set for a failed primary set.',
'NicheNet v2 mouse priors: Zenodo7074291; downloaded checksums verified. Ranking is corrected AUPR, not a P value or probability of causality. Downregulated gene-set activity does not mean ligand inhibits those genes: the prior is unsigned.',
'Candidate prior ranks use receiver background. Targets outside that background are flagged. Missing prior = unassessable, not no biology. Gene aliases mapped and duplicate counts summed with mapping audit.',
'Output network combines expressed LR edges with ligand-target prior, not a verified receptor-target signaling path. Co-receptor co-detection is an additional audit; dropouts do not prove absence. Transcript scores do not prove mitochondrial injury, protein cleavage, or cell death.',
'Human source-adjusted module correlations are exploratory; subtype composition and leave-one-sample-out checks are sensitivity analyses. Mixed-source human data are not an independent replicated single-cohort disease validation.',
'Official method: https://github.com/saeyslab/nichenetr/blob/master/vignettes/ligand_activity_geneset.md'),file.path(out,'methods.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
cat('NICHE_NET_ANALYSIS_PASS\n')

source(file.path(out,"pairing.R"))

R_NICHE
cat > "$OUT/finish.py" <<'PY_FINISH'
import os,json,tarfile,hashlib
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path(os.environ['OUT']);root=Path(os.environ['ROOT'])
required=['DE_status.tsv','geneset_status.tsv','ligand_activities.tsv','candidate_prior_coverage.tsv','candidate_ligand_prior.tsv','supported_LR_edges.tsv','human_composition_correlations.tsv','human_correlation_LOO.tsv','human_source_models_BH.tsv','receptor_complex_codetection.tsv','paired_ligand_target_correlations.tsv']
for f in required:
    assert (out/'tables'/f).is_file(),f
st=pd.read_csv(out/'tables/geneset_status.tsv',sep='\t')
assert len(st)==16,st
ac=pd.read_csv(out/'tables/ligand_activities.tsv',sep='\t')
assert len(st[st.status=='PASS'])==(ac[['receiver','filter','test','direction']].drop_duplicates().shape[0] if len(ac) else 0)
captions=[]
if len(ac):
    keys=['receiver','test','direction','test_ligand']
    merged=ac[ac['filter']=='All'][keys+['rank','aupr_corrected']].merge(ac[ac['filter']=='BG1'][keys+['rank','aupr_corrected']],on=keys,how='outer',suffixes=('_All','_BG1'))
    merged['top20_both']=(merged['rank_All']<=20)&(merged['rank_BG1']<=20)
    merged.to_csv(out/'tables/ligand_rank_sensitivity.tsv',sep='\t',index=False,na_rep='NA')
else:
    pd.DataFrame(columns=['receiver','test','direction','test_ligand','rank_All','rank_BG1','top20_both']).to_csv(out/'tables/ligand_rank_sensitivity.tsv',sep='\t',index=False)
if len(ac):
    for (receiver,direction),v in ac[(ac['filter']=='All')&(ac.test=='FDR05')].groupby(['receiver','direction']):
        v=v.sort_values(['rank','test_ligand']).head(15).iloc[::-1]
        labels=[x if len(x)<=10 else 'L'+str(i+1) for i,x in enumerate(v.test_ligand)]
        fig,ax=plt.subplots(figsize=(6.8,5.4));ax.barh(labels,v.aupr_corrected,color='#3979a8')
        ax.set_xlabel('AUPR adj');ax.set_title(receiver+' / '+direction);ax.axvline(0,color='gray',linewidth=.5)
        fig.tight_layout();stem=f'NicheNet_{receiver}_{direction}';fig.savefig(out/'figures'/f'{stem}.png',dpi=300);fig.savefig(out/'figures'/f'{stem}.pdf');plt.close(fig)
        captions.append(stem+': FDR-defined receiver gene set; corrected AUPR ligand ranking. '+ '; '.join(f'{a}={b}' for a,b in zip(labels,v.test_ligand)))
(out/'captions.txt').write_text('\n'.join(captions))
summary=dict(status='PASS',stage='mechanism_step3_nichenet',output=str(out),previous=os.environ['PREV'],previous2=os.environ['PREV2'],nichenet_sets_run=int((st.status=='PASS').sum()),nichenet_sets_skipped=int((st.status!='PASS').sum()),note='PASS means workflow completed, not proof of any communication axis or cell death.')
(out/'status.json').write_text(json.dumps(summary,indent=2))
archive=root/'tmp/lung_mechanism_step3_return.tar.gz';part=archive.with_suffix('.tar.gz.part')
files=[p for p in out.rglob('*') if p.is_file() and 'objects' not in p.relative_to(out).parts]
with tarfile.open(part,'w:gz') as tar:
    for p in sorted(files):tar.add(p,arcname=str(p.relative_to(out)))
with tarfile.open(part,'r:gz') as tar:
    for m in tar:
        if m.isfile():
            h=tar.extractfile(m)
            while h.read(1024*1024):pass
part.replace(archive)
print('RETURN='+str(archive),flush=True)
print('RETURN_BYTES='+str(archive.stat().st_size),flush=True)
print('RETURN_SHA256='+hashlib.sha256(archive.read_bytes()).hexdigest(),flush=True)
print('NICHE_NET_SETS_RUN='+str(summary['nichenet_sets_run']),flush=True)
print('LUNG_MECHANISM_STEP3_PASS',flush=True)

PY_FINISH

conda run --no-capture-output -n batch_bench python -m py_compile "$OUT/fetch.py" "$OUT/finish.py"
conda run --no-capture-output -n R_4.5.2 Rscript -e 'for(f in commandArgs(TRUE)) invisible(parse(file=f)); cat("R_PARSE_PASS\n")' "$OUT/analyze.R" "$OUT/pairing.R" "$OUT/human_audit.R"
echo SYNTAX_CHECK=PASS
if [[ "${SYNTAX_ONLY:-0}" == 1 ]]; then exit 0; fi
conda run --no-capture-output -n R_4.5.2 Rscript -e 'p <- c("Matrix","edgeR","statmod","nichenetr"); ok <- vapply(p,requireNamespace,logical(1),quietly=TRUE); print(data.frame(package=p,available=ok)); if(!all(ok)) stop("Required package missing; do not replace with cell-level tests.")'
conda run --no-capture-output -n batch_bench python -u "$OUT/fetch.py"
conda run --no-capture-output -n R_4.5.2 Rscript "$OUT/human_audit.R" "$PREV2" "$OUT"
conda run --no-capture-output -n R_4.5.2 Rscript "$OUT/analyze.R" "$PREV" "$PREV2" "$OUT" "$CACHE"
conda run --no-capture-output -n batch_bench python -u "$OUT/finish.py"
