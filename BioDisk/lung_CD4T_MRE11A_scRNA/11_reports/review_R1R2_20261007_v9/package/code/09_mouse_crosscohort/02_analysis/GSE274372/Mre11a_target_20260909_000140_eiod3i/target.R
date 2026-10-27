suppressPackageStartupMessages(library(edgeR))
root<-Sys.getenv('ROOT');run<-Sys.getenv('RUN');out<-file.path(run,'tables');dir.create(out)
src<-file.path(root,'02_analysis/GSE274372/stage10_20260908_141459_a8IPaq')
f<-read.delim(gzfile(file.path(src,'tables/features.tsv.gz')),row.names=1)
gene<-rownames(f)[which(f$symbol=='Mre11a')];stopifnot(length(gene)==1)
meta<-read.delim(file.path(src,'tables/sample_summary.tsv'))
results<-list();descriptions<-list();fits<-list()
# TLike is the primary candidate scope; Author is a prespecified scope sensitivity.
for(scope in c('TLike','Author')) {
 x<-readRDS(file.path(src,paste0(scope,'_edgeR.rds')))
 raw<-as.matrix(read.delim(gzfile(file.path(src,'tables',paste0(scope,'_counts.tsv.gz'))),row.names=1,check.names=FALSE))
 target<-raw[gene,colnames(x$y),drop=FALSE]
 s<-meta[meta$scope==scope,];s$included<-s$sample %in% colnames(x$y)
 s$detect_pct<-100*s$Mre_positive/s$cells
 s$raw_CPM<-1e6*s$Mre_UMI/s$library_UMI
 descriptions[[scope]]<-s
 stopifnot(sum(target)>0,all(target>=0),all(target==round(target)))
 mat<-rbind(x$y$counts[setdiff(rownames(x$y),gene),,drop=FALSE],target)
 # Keep the original normalization scale. Add only the prespecified target.
 y<-DGEList(counts=mat,lib.size=x$y$samples$lib.size,norm.factors=x$y$samples$norm.factors)
 design<-x$design
 y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE)
 cm<-cbind(c(-1,1,0,0),c(-1,0,1,0),c(-1,0,0,1))
 tests<-c(list(Omnibus=glmQLFTest(fit,contrast=cm)),setNames(lapply(1:3,function(i)glmQLFTest(fit,contrast=cm[,i])),c('D7_D0','D28_D0','D56_D0')))
 for(nm in names(tests)) {
  z<-tests[[nm]]$table[gene,,drop=FALSE]
  lfc<-if('logFC'%in%names(z))z$logFC else NA_real_
  results[[paste(scope,nm)]]<-data.frame(scope=scope,test=nm,gene='Mre11a',samples=ncol(y),target_UMI=sum(target),log2FC=lfc,PValue=z$PValue,nominal_p05=z$PValue<.05)
 }
 fits[[scope]]<-list(y=y,fit=fit,design=design)
}
r<-do.call(rbind,results);rownames(r)<-NULL
write.table(r,file.path(out,'Mre11a_target_tests.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
s<-do.call(rbind,descriptions);rownames(s)<-NULL
write.table(s,file.path(out,'Mre11a_sample_detection.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(c('Prespecified Mre11a target analysis: bypass expression filtering for this gene only; other genes remain the original filtered set for dispersion borrowing. edgeR robust negative-binomial quasi-likelihood model with categorical day. Original Stage10 samples and normalization retained. Primary scope TLike; Author sensitivity. Omnibus tests any day effect, not monotonicity. Day-versus-baseline contrasts are separate exploratory nominal p-values; do not select only the lowest p. Sparse counts and candidate identity limit interpretation. Report all contrasts and sample-level detection. Expression changes alone do not establish Mre11a function or causality.'),file.path(out,'methods.txt'))
saveRDS(fits,file.path(run,'target_models.rds'))
print(r)
cat('MRE11A_TARGET_ANALYSIS_PASS\n')
