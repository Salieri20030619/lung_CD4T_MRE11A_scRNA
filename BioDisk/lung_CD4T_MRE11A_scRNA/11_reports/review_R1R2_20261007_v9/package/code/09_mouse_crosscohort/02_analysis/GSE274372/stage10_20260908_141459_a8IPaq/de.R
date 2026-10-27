suppressPackageStartupMessages(library(edgeR))
out <- file.path(Sys.getenv('RUN'),'tables')
meta <- read.delim(file.path(out,'sample_summary.tsv'),check.names=FALSE)
features <- read.delim(gzfile(file.path(out,'features.tsv.gz')),row.names=1,check.names=FALSE)
status <- list();target <- list()
for (scope in c('TLike','NoMix','Author')) {
 m <- meta[meta$scope==scope,];m$included <- m$cells>=10 & m$library_UMI>0
 write.table(m,file.path(out,paste0(scope,'sample_inclusion.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
 m <- m[m$included,];m$day <- factor(m$day,levels=c(0,7,28,56))
 n <- table(m$day)
 if(any(n<2)) {status[[scope]]<-data.frame(scope=scope,status='SKIP: fewer than2 samples in a day');next}
 counts <- as.matrix(read.delim(gzfile(file.path(out,paste0(scope,'_counts.tsv.gz'))),row.names=1,check.names=FALSE))[,m$sample,drop=FALSE]
 design <- model.matrix(~0+day,data=m);colnames(design)<-paste0('D',levels(m$day))
 stopifnot(qr(design)$rank==ncol(design),nrow(design)>ncol(design))
 y <- DGEList(counts=counts)
 keep <- filterByExpr(y,design=design)
 audit <- data.frame(gene=rownames(y),symbol=features[rownames(y),'symbol'],tested=keep)
 write.table(audit,file.path(out,paste0(scope,'gene_filter.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
 target[[scope]] <- cbind(scope=scope,audit[audit$symbol %in% 'Mre11a',])
 y <- calcNormFactors(y[keep,,keep.lib.sizes=FALSE]);y<-estimateDisp(y,design,robust=TRUE)
 fit <- glmQLFit(y,design,robust=TRUE)
 # Three day-versus-baseline contrasts jointly test any time effect.
 cm <- cbind(c(-1,1,0,0),c(-1,0,1,0),c(-1,0,0,1))
 tests <- c(list(Omnibus=glmQLFTest(fit,contrast=cm)),setNames(lapply(1:3,function(i)glmQLFTest(fit,contrast=cm[,i])),c('D7_D0','D28_D0','D56_D0')))
 for (nm in names(tests)) {
  z<-topTags(tests[[nm]],n=Inf,sort.by='none')$table
  z<-cbind(gene=rownames(z),symbol=features[rownames(z),'symbol'],z)
  write.table(z,file.path(out,paste0(scope,'_',nm,'.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
 }
 saveRDS(list(y=y,fit=fit,design=design,meta=m),file.path(Sys.getenv('RUN'),paste0(scope,'_edgeR.rds')))
 status[[scope]]<-data.frame(scope=scope,status=paste('PASS samples',nrow(m),'genes',sum(keep)))
}
write.table(do.call(rbind,status),file.path(out,'DE_status.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
if(length(target))write.table(do.call(rbind,target),file.path(out,'Mre_test_eligibility.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'R_session.txt'))
cat('GSE274372_STAGE10_ANALYSIS_PASS\n')
