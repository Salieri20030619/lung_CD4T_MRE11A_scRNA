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

