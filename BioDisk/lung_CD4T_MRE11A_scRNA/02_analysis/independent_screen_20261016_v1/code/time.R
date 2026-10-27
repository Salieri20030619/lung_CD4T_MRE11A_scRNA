suppressPackageStartupMessages(library(edgeR))
args<-commandArgs(trailingOnly=TRUE);out<-args[1]
p<-'/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/02_analysis/GSE274372/stage10_20260908_141459_a8IPaq/tables'
meta<-read.delim(file.path(p,'sample_summary.tsv'),check.names=FALSE)
m<-meta[meta$scope=='TLike',];stopifnot(nrow(m)==12,all(m$cells>0),all(m$library_UMI>0))
m$included<-TRUE;m$inclusion_reason<-'nonempty measurable sample; no minimum-cell threshold'
m$day<-factor(m$day,levels=c(0,7,28,56));stopifnot(all(table(m$day)==3))
counts<-as.matrix(read.delim(gzfile(file.path(p,'TLike_counts.tsv.gz')),row.names=1,check.names=FALSE))[,m$sample,drop=FALSE]
features<-read.delim(gzfile(file.path(p,'features.tsv.gz')),row.names=1,check.names=FALSE)
stopifnot(all(is.finite(counts)),all(counts>=0),all(abs(counts-round(counts))<1e-8))
stopifnot(all(abs(colSums(counts)-m$library_UMI)<1e-6))
design<-model.matrix(~0+day,data=m);colnames(design)<-paste0('D',levels(m$day))
y<-DGEList(counts=counts,lib.size=m$library_UMI)
keep<-filterByExpr(y,design=design)
write.table(data.frame(gene=rownames(counts),symbol=features[rownames(counts),'symbol'],tested=keep),file.path(out,'time','gene_filter.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
y<-calcNormFactors(y[keep,,keep.lib.sizes=TRUE]);y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE)
cm<-cbind(c(-1,1,0,0),c(-1,0,1,0),c(-1,0,0,1))
tests<-c(list(Omnibus=glmQLFTest(fit,contrast=cm)),setNames(lapply(1:3,function(i)glmQLFTest(fit,contrast=cm[,i])),c('D7_D0','D28_D0','D56_D0')))
tests$D28_D7<-glmQLFTest(fit,contrast=c(0,-1,1,0))
tests$D56_D28<-glmQLFTest(fit,contrast=c(0,0,-1,1))
for(nm in names(tests)){
 z<-topTags(tests[[nm]],n=Inf,sort.by='none')$table
 z<-cbind(gene=rownames(z),symbol=features[rownames(z),'symbol'],z)
 z$family<-paste0('GSE274372_all12_TLike_',nm,'_all_tested_genes')
 write.table(z,file.path(out,'time',paste0(nm,'.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
}
z<-cbind(gene=rownames(y),symbol=features[rownames(y),'symbol'],as.data.frame(cpm(y,log=TRUE,prior.count=2),check.names=FALSE))
write.table(z,file.path(out,'time','sample_TMM_logCPM.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(m,file.path(out,'time','sample_inclusion.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(cbind(sample=m$sample,as.data.frame(design)),file.path(out,'time','design.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(cbind(sample=m$sample,y$samples),file.path(out,'time','normalization.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
saveRDS(list(y=y,fit=fit,design=design,meta=m),file.path(out,'time','all12_TLike_edgeR.rds'))
writeLines(c('Real categorical time, not pseudotime. No contemporaneous late saline controls; age and exposure duration not separable.','All 12 nonempty TLike libraries included; D0 8-cell library retained. TLike is not confirmed pure CD4.','Full original count-library sizes retained after gene filtering; TMM factors recalculated; new results do not overwrite legacy min10 model.','BH separately within complete tested-gene families for omnibus and each contrast.','Nominal selection uses original P<0.05; no fold-change threshold tuned to retain targets.',capture.output(sessionInfo())),file.path(out,'time','methods_and_session.txt'))
cat('ALL12_TIME_MODEL_PASS samples=',nrow(m),'genes=',sum(keep),'\n')
for(nm in names(tests)){z<-topTags(tests[[nm]],n=Inf)$table;cat(nm,'P05=',sum(z$PValue<.05),'BH05=',sum(z$FDR<.05),'\n')}