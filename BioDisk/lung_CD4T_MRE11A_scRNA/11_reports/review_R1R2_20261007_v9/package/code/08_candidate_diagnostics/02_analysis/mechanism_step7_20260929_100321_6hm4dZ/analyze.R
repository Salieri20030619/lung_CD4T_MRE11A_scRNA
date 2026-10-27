args <- commandArgs(TRUE);p1 <- args[1];p2 <- args[2];p6 <- args[3];out <- args[4]
suppressPackageStartupMessages({library(Matrix);library(edgeR)})
rd <- function(p)read.delim(p,check.names=FALSE,stringsAsFactors=FALSE)
wr <- function(x,n)write.table(x,file.path(out,'tables',paste0(n,'.tsv')),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
canon <- function(x){x <- toupper(trimws(x));a <- c(MRE11A='MRE11',ATP5A1='ATP5F1A',ATP5B='ATP5F1B',ATP5C1='ATP5F1C',ATP5G1='ATP5MC1',ATP5O='ATP5PO',PARK2='PRKN');i <- x%in%names(a);x[i] <- unname(a[x[i]]);x}
TFs <- c('STAT1','STAT3','STAT5A','STAT5B','NFKB1','RELA','TBX21','BATF')
targets <- c('IL1R2','PENK','TTC39C','FOXP1','TNFRSF9','LAG3','MRE11')
# FUNCTIONS_START
add_composition <- function(d,cell_metadata){
 stopifnot(!anyDuplicated(d$sample),!anyNA(cell_metadata$sample),!anyNA(cell_metadata$subtype),all(cell_metadata$sample%in%d$sample))
 tt <- table(factor(cell_metadata$sample,levels=d$sample),factor(cell_metadata$subtype))
 # as.matrix(table) can retain class 'table'; rebuild a plain numeric matrix.
 counts <- matrix(as.numeric(tt),nrow=nrow(tt),ncol=ncol(tt),dimnames=dimnames(tt))
 stopifnot(is.matrix(counts),!inherits(counts,'table'),identical(rownames(counts),as.character(d$sample)),all(rowSums(counts)>0),all(rowSums(counts)==d$cells))
 props <- sweep(counts,1,rowSums(counts),'/')
 stopifnot(all(is.finite(props)),all(abs(rowSums(props)-1)<1e-10))
 if(ncol(props)>1){
  vv <- props[,-1,drop=FALSE];colnames(vv) <- paste0('comp_',seq_len(ncol(vv)))
  for(j in seq_len(ncol(vv)))d[[colnames(vv)[j]]] <- as.numeric(vv[,j])
 }
 stopifnot(nrow(d)==nrow(counts),!anyDuplicated(d$sample))
 d
}
check_response <- function(y,d){
 if(length(y)!=nrow(d))stop('RESPONSE_METADATA_DIMENSION_MISMATCH: response=',length(y),' metadata=',nrow(d))
 if(!is.null(names(y))&&!identical(names(y),as.character(d$sample)))stop('RESPONSE_SAMPLE_ORDER_MISMATCH')
 invisible(TRUE)
}
ulm <- function(expr,network,exclude=character(),minsize=10){
 expr <- expr[!rownames(expr)%in%exclude,,drop=FALSE];if(nrow(expr)<500)return(NULL)
 net <- network[network$target%in%rownames(expr),,drop=FALSE];sz <- table(net$tf);tfs <- names(sz)[sz>=minsize];if(!length(tfs))return(NULL)
 ans <- matrix(NA_real_,length(tfs),ncol(expr),dimnames=list(tfs,colnames(expr)));Y <- sweep(expr,2,colMeans(expr),'-');sy <- colSums(Y^2);n <- nrow(expr)
 for(tf in tfs){q <- net[net$tf==tf,];w <- numeric(n);w[match(q$target,rownames(expr))] <- q$mor;w <- w-mean(w);sw <- sum(w*w);if(sw<1e-12)next;r <- as.numeric(crossprod(w,Y))/sqrt(sw*sy);r[sy<1e-12] <- NA_real_;r <- pmax(-1+1e-12,pmin(1-1e-12,r));ans[tf,] <- r*sqrt((n-2)/(1-r*r))}
 list(scores=ans,sizes=as.integer(sz[tfs]),universe=n)
}
make_design <- function(d,arm,rank_depth=FALSE){
 stopifnot(!anyDuplicated(d$sample),!anyNA(d$group),!anyNA(d$series),all(d$group%in%c('Control','CWP')),all(is.finite(d$library)),all(d$library>0))
 rownames(d) <- as.character(d$sample)
 d$group <- factor(d$group,levels=c('Control','CWP'));d$series <- factor(d$series);d$depth <- if(rank_depth)rank(d$library) else log1p(d$library)
 terms <- c('group',if(arm!='WithinG227'&&nlevels(d$series)>1)'series',if(arm%in%c('SourceDep','SourceComp'))'depth',if(arm=='SourceComp')grep('^comp_',names(d),value=TRUE))
 form <- reformulate(terms);mf <- model.frame(form,data=d,na.action=na.fail)
 X <- model.matrix(form,mf)
 stopifnot(nrow(X)==nrow(d),identical(rownames(X),as.character(d$sample)),all(is.finite(X)))
 X
}
fit_effect <- function(y,d,arm){
 check_response(y,d)
 sel <- is.finite(y);if(arm=='WithinG227')sel <- sel&d$series=='G227';d <- d[sel,,drop=FALSE];y <- y[sel]
 nc <- sum(d$group=='Control');nd <- sum(d$group=='CWP');ans <- c(coef=NA_real_,p=NA_real_,df=NA_real_,n_control=nc,n_case=nd)
 if(min(nc,nd)<2)return(ans)
 X <- make_design(d,arm);stopifnot(nrow(X)==length(y));if(qr(X)$rank<ncol(X)||nrow(X)-ncol(X)<3)return(ans)
 f <- lm.fit(X,y);df <- length(y)-f$rank;ss <- sum(f$residuals^2)/df;se <- sqrt(diag(solve(crossprod(X)))*ss);j <- match('groupCWP',colnames(X));if(is.na(j)||se[j]<1e-12)return(ans)
 ans['coef'] <- f$coefficients[j];ans['p'] <- 2*pt(-abs(f$coefficients[j]/se[j]),df);ans['df'] <- df;ans
}
partial <- function(x,y,d,arm){
 check_response(x,d);check_response(y,d)
 ok <- is.finite(x)&is.finite(y);x <- x[ok];y <- y[ok];d <- d[ok,,drop=FALSE];n <- length(x);ans <- c(rho=NA_real_,p=NA_real_,n=n,df=NA_real_)
 if(n<10||length(unique(d$group))!=2||min(table(d$group))<2)return(ans)
 X <- make_design(d,arm,TRUE);stopifnot(nrow(X)==length(x),length(x)==length(y));df <- n-qr(X)$rank-1L;if(qr(X)$rank<ncol(X)||df<3)return(ans)
 u <- qr.resid(qr(X),rank(x));v <- qr.resid(qr(X),rank(y));if(sum(u*u)<1e-10||sum(v*v)<1e-10)return(ans)
 r <- max(-1,min(1,cor(u,v)));ans['rho'] <- r;ans['df'] <- df;ans['p'] <- 2*pt(-abs(r)*sqrt(df/max(1-r*r,1e-15)),df);ans
}
# FUNCTIONS_END
net <- rd(file.path(out,'network.tsv'));net$tf <- canon(net$tf);net$target <- canon(net$target)
net <- net[net$tf%in%TFs&net$confidence%in%c('A','B','C'),];stopifnot(all(is.finite(net$mor)),all(net$mor!=0))
k <- paste(net$tf,net$target);bad <- names(which(vapply(split(sign(net$mor),k),function(z)length(unique(z))>1,logical(1))));wr(net[k%in%bad,],'conflicting_edges');net <- net[!k%in%bad,];net <- net[order(net$confidence),];net <- net[!duplicated(net[,c('tf','target')]),];wr(net,'human_network')
x <- as(readMM(file.path(p2,'input/human_counts.mtx')),'CsparseMatrix');g <- rd(file.path(p2,'input/human_genes.tsv'));m <- rd(file.path(p2,'input/human_cells.tsv'))
stopifnot(nrow(x)==nrow(g),ncol(x)==nrow(m),ncol(x)==6523,!anyDuplicated(m$cell),all(is.finite(x@x)),all(x@x>=0),all(abs(x@x-round(x@x))<1e-6),all(m$library>0))
cb <- toupper(as.character(m$core));stopifnot(all(cb%in%c('TRUE','FALSE')));m$core <- cb=='TRUE'
sy <- canon(g$symbol);keep <- !is.na(sy)&nzchar(sy)&!grepl('^ENS',sy);u <- sort(unique(sy[keep]));pr <- sparseMatrix(i=match(sy[keep],u),j=which(keep),x=1,dims=c(length(u),nrow(x)));x <- pr%*%x;rownames(x) <- u;colnames(x) <- m$cell
stopifnot(all(Matrix::colSums(x)<=m$library+1e-6));wr(data.frame(original=g$symbol,canonical=sy,retained=keep),'gene_mapping')
bgp <- strsplit('SFTPA1 SFTPA2 SFTPB SFTPC SFTPD SCGB1A1 SCGB3A1 SCGB3A2 ABCA3 SLC34A2 NAPSA LYZ LYZ1 LYZ2 S100A8 S100A9 CSF1R LST1 TYROBP FCER1G',' ')[[1]]
b <- u%in%bgp|grepl('^(IGH[VDJMAEG]|IGK[VJC]|IGL[VJC])',u);m$bg <- as.numeric(Matrix::colSums(x[b,,drop=FALSE]))/m$library
ss <- rd(file.path(p2,'tables/human_sample_scores.tsv'));mods <- rd(file.path(p1,'tables/module_definition.tsv'));mods$gene <- canon(mods$gene)
activities <- list();effects <- list();associations <- list();cover <- list();audits <- list();edges <- list();skips <- list();exprs <- list();objects <- list()
for(scope in unique(ss$scope))for(filt in c('All','BG1')){
 sel <- if(scope=='Primary')rep(TRUE,nrow(m)) else if(scope=='Core')m$core else m$subtype==sub('^Subtype:','',scope)
 if(filt=='BG1')sel <- sel&m$bg<=.01
 md <- m[sel,];xx <- x[,sel,drop=FALSE];q <- ss[ss$scope==scope&ss$filter==filt&ss$cells>=20,];q <- q[order(q$sample),];key <- paste(scope,filt,sep='|')
 if(nrow(q)<6||length(unique(q$group))!=2||min(table(q$group))<2){skips[[key]] <- data.frame(scope,filter=filt,reason='INSUFFICIENT_SAMPLE_COUNTS');next}
 ns <- table(factor(md$sample,levels=q$sample));stopifnot(identical(as.integer(ns),as.integer(q$cells)))
 j <- match(md$sample,q$sample);ok <- !is.na(j);md <- md[ok,];xx <- xx[,ok,drop=FALSE];j <- j[ok]
 A <- sparseMatrix(i=seq_len(nrow(md)),j=j,x=1,dims=c(nrow(md),nrow(q)));pb <- as.matrix(xx%*%A);colnames(pb) <- q$sample
 lib <- as.numeric(tapply(md$library,factor(md$sample,levels=q$sample),sum));stopifnot(all(lib>0),all(abs(lib/q$cells-q$library)<1e-6))
 ident <- unique(md[,c('sample','group','series')]);stopifnot(!anyDuplicated(ident$sample));ident <- ident[match(q$sample,ident$sample),];stopifnot(identical(ident$group,q$group),identical(ident$series,q$series))
 # Keep all but a fixed reference subtype as composition covariates for pooled scopes.
 if(scope%in%c('Primary','Core')){
  q <- add_composition(q,md)
 }
 stopifnot(ncol(pb)==nrow(q),identical(colnames(pb),as.character(q$sample)))
 cat('MODEL_INPUT_AUDIT scope=',scope,' filter=',filt,' count_columns=',ncol(pb),' metadata_rows=',nrow(q),'\n',sep='');flush.console()
 wr(q[,c('sample','group','series','cells','library',grep('^comp_',names(q),value=TRUE)),drop=FALSE],paste0('design_',gsub('[^A-Za-z0-9_]','_',key)))
 design <- make_design(q,'Source');if(qr(design)$rank<ncol(design)){skips[[key]] <- data.frame(scope,filter=filt,reason='SOURCE_DESIGN_NOT_IDENTIFIABLE');next}
 y <- DGEList(pb,lib.size=lib);tested <- filterByExpr(y,design=design);if(sum(tested)<500){skips[[key]] <- data.frame(scope,filter=filt,reason='FEWER_THAN_500_GENES');next}
 y <- calcNormFactors(y[tested,,keep.lib.sizes=TRUE]);E <- cpm(y,log=TRUE,prior.count=2);E <- E-rowMeans(E)
 audits[[key]] <- q[,c('scope','filter','sample','group','series','cells','library',grep('^comp_',names(q),value=TRUE)),drop=FALSE]
 for(tf in unique(c(TFs,targets))){i <- match(tf,u);exprs[[length(exprs)+1]] <- data.frame(scope,filter=filt,sample=q$sample,gene=tf,present=!is.na(i),raw_sum=if(!is.na(i))pb[i,] else NA_real_,tested=tf%in%rownames(E))}
 for(level in c('ABC','AB')){
  nn <- net[net$confidence%in%strsplit(level,'')[[1]],];z <- ulm(E,nn,targets);if(is.null(z))next
  V <- z$scores;stopifnot(all(is.finite(V)))
  cover[[length(cover)+1]] <- data.frame(scope,filter=filt,level,exclusion='Drop7',tf=rownames(V),n_targets=z$sizes,universe=z$universe)
  for(tf in rownames(V)){
   v <- V[tf,];activities[[length(activities)+1]] <- data.frame(scope,filter=filt,level,tf,sample=q$sample,group=q$group,series=q$series,activity=as.numeric(v))
   arms <- c('Source','SourceDep','WithinG227',if(scope%in%c('Primary','Core'))'SourceComp')
   for(arm in arms){r <- fit_effect(v,q,arm);loo <- vapply(seq_along(v),function(i)fit_effect(v[-i],q[-i,,drop=FALSE],arm)['coef'],numeric(1));finite <- is.finite(loo)
    effects[[length(effects)+1]] <- data.frame(scope,filter=filt,level,tf,model=arm,coefficient=r['coef'],p_nominal=r['p'],residual_df=r['df'],n_control=r['n_control'],n_case=r['n_case'],LOO_min=if(any(finite))min(loo[finite]) else NA_real_,LOO_max=if(any(finite))max(loo[finite]) else NA_real_,LOO_n=sum(finite),row.names=NULL)
   }
   for(gene in targets)for(arm in c('SourceDep',if(scope%in%c('Primary','Core'))'SourceComp')){
    val <- q[[paste0('expr_',gene)]];stopifnot(length(val)==length(v));r <- partial(v,val,q,arm);if(sum(val>0,na.rm=TRUE)<4)r[c('rho','p')] <- NA_real_
    associations[[length(associations)+1]] <- data.frame(scope,filter=filt,level,tf,model=arm,kind='Candidate',feature=gene,exclusion='Drop7',rho=r['rho'],p_nominal=r['p'],n=r['n'],residual_df=r['df'],row.names=NULL)
   }
  }
  for(mod in unique(mods$module)){
   zz <- ulm(E,nn,union(targets,mods$gene[mods$module==mod]));if(is.null(zz))next
   cover[[length(cover)+1]] <- data.frame(scope,filter=filt,level,exclusion=paste0('Drop7Plus_',mod),tf=rownames(zz$scores),n_targets=zz$sizes,universe=zz$universe)
   for(tf in rownames(zz$scores))for(arm in c('SourceDep',if(scope%in%c('Primary','Core'))'SourceComp')){
    val <- q[[paste0('score_',mod)]];v <- zz$scores[tf,];stopifnot(length(val)==length(v));r <- partial(v,val,q,arm)
    associations[[length(associations)+1]] <- data.frame(scope,filter=filt,level,tf,model=arm,kind='Module',feature=mod,exclusion='Drop7PlusModule',rho=r['rho'],p_nominal=r['p'],n=r['n'],residual_df=r['df'],row.names=NULL)
   }
  }
  ee <- nn[nn$tf%in%rownames(V)&nn$target%in%c(targets,mods$gene),];if(nrow(ee)){ee$scope <- scope;ee$filter <- filt;ee$level <- level;ee$measured <- ee$target%in%rownames(x);ee$tested <- ee$target%in%rownames(E);ee$modules <- vapply(ee$target,function(g)paste(unique(mods$module[mods$gene==g]),collapse=';'),character(1));edges[[length(edges)+1]] <- ee}
  objects[[paste(key,level)]] <- list(activity=V,metadata=q)
 }
 cat('HUMAN_TF_SCOPE_DONE=',key,' samples=',nrow(q),'\n',sep='');flush.console()
}
bind <- function(z)if(length(z))do.call(rbind,z) else data.frame(note=character())
# Covariate columns vary by scope; record common audit columns separately.
basecols <- c('scope','filter','sample','group','series','cells','library');wr(bind(lapply(audits,function(x)x[,basecols])),'human_sample_audit')
wr(bind(activities),'human_TF_activity');wr(bind(cover),'human_TF_coverage');wr(bind(edges),'human_TF_edges');wr(bind(exprs),'human_TF_RNA');wr(bind(skips),'human_skips')
er <- bind(effects);er$FDR <- NA_real_;er$FDR_status <- 'GROUP_N_LT10';for(i in split(seq_len(nrow(er)),interaction(er$scope,er$filter,er$level,er$model,drop=TRUE))){j <- i[er$n_control[i]>=10&er$n_case[i]>=10&is.finite(er$p_nominal[i])];if(length(j)){er$FDR[j] <- p.adjust(er$p_nominal[j],'BH');er$FDR_status[j] <- 'BH'}};wr(er,'human_TF_effects')
ar <- bind(associations);ar$FDR <- NA_real_;for(i in split(seq_len(nrow(ar)),interaction(ar$scope,ar$filter,ar$level,ar$model,drop=TRUE))){j <- i[is.finite(ar$p_nominal[i])];if(length(j))ar$FDR[j] <- p.adjust(ar$p_nominal[j],'BH')};wr(ar,'human_TF_associations')
saveRDS(objects,file.path(out,'objects/human_TF_scores.rds'))
writeLines(c('Human lung official combined cohort previously exported in step2; no new donor added. >=20 cells per scope and donor, existing annotations retained. Human subtypes are not assumed identical to mouse ActCD4.',
'Fixed eight TF panel selected from mouse step6: STAT1,STAT3,STAT5A,STAT5B,NFKB1,RELA,TBX21,BATF. STAT5B is pathway comparator, not a mouse positive. Human-specific DoRothEA ABC and AB; immutable prior embedded. ULM as in step6.',
'Counts summed per subject; full raw library offset and TMM, filterByExpr(~group+series). Gene-centered logCPM, prior.count2. All and BG1 use same original full-count denominators. Drop7 removes candidates from regulons and expression universe; module associations additionally exclude module members.',
'Effect models Source: group+series; SourceDep additionally log1p(mean cell library); SourceComp additionally all but one subtype proportion in pooled scopes; WithinG227 uses mixed-condition source only. Invalid/rank-deficient designs produce NA, not simpler replacement models.',
'Group nominal P retained; BH only when both groups>=10. Correlations require total>=10, >=2 per group, residual df>=3 and >=4 positive subjects for candidate expression; rank residualization on group, series, ranked depth and optional composition; BH across all TF-feature tests within scope/filter/level/model. No cell-level P values.',
'LOO group models remove each donor and refit on frozen activity scores. Regulon activity is inferred and not direct DNA binding. Shared module control genes may remain. Mouse-human comparison concerns direction; activity magnitudes are not pooled.',
'Human macrophage object inventory is path/size metadata only. No human macrophage axis validation or final novelty claim is completed by this stage.'),file.path(out,'methods.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'));cat('HUMAN_TF_ANALYSIS_PASS\n')

