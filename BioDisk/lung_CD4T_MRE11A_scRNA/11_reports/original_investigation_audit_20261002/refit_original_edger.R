root <- '/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA'
out <- file.path(root,'11_reports/original_investigation_audit_20261002')
dir.create(file.path(out,'results'),recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages(library(edgeR))
src <- file.path(root,'02_analysis/GSE314943/stage5_20260905_222457_a2dApT')
old <- file.path(root,'02_analysis/GSE314943/stage6_20260905_231807_L44XdT')
rd <- function(p)read.delim(p,check.names=FALSE,stringsAsFactors=FALSE)
wr <- function(d,n)write.table(d,file.path(out,'results',n),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
genes <- rd(file.path(src,'tables/genes.tsv'));meta <- rd(file.path(src,'tables/PB_sample_metadata.tsv'))
samples <- c(paste0('Saline_',1:3),paste0('Silica_',1:3));group <- factor(c(rep('Saline',3),rep('Silica',3)),levels=c('Saline','Silica'))
design <- model.matrix(~group);rownames(design)<-samples
background <- c('Sftpc','Sftpa1','Sftpb','Sftpd','Scgb1a1','Scgb3a2','S100a8','S100a9','Lyz2','Igkc','Ighm','Iglc1','Iglc2','Iglc3')
excluded <- genes$symbol %in% background | grepl('^Hb[ab]',genes$symbol)
targets <- c('Il1r2','Penk','Ttc39c','Foxp1','Tnfrsf9','Lag3','Ikzf2','Stat1')
summaries<-list();focal<-list();offsets<-list()
for(scope in c('Core','ActCD4','Naive','Treg')){
 md<-meta[meta$set==scope,];md<-md[match(samples,md$sample),];stopifnot(identical(md$sample,samples),identical(md$group,as.character(group)))
 dat<-rd(gzfile(file.path(src,'tables',paste0('PB_',scope,'_counts.tsv.gz'))));counts<-as.matrix(dat[,-1]);rownames(counts)<-dat[[1]]
 stopifnot(identical(colnames(counts),samples),identical(rownames(counts),genes$gene_id),all(colSums(counts)==md$library_umi),all(counts>=0),all(counts==round(counts)))
 for(arm in c('Full','NoBG')){
  inc<-if(arm=='Full')rep(TRUE,nrow(genes)) else !excluded
  y<-DGEList(counts=counts[inc,,drop=FALSE],group=group,genes=genes[inc,,drop=FALSE]);keep<-filterByExpr(y,design=design)
  y<-y[keep,,keep.lib.sizes=FALSE];y<-calcNormFactors(y,method='TMM');y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE)
  tt<-topTags(glmQLFTest(fit,coef=2),n=Inf,sort.by='PValue')$table
  tag<-paste(scope,arm,sep='_');prior<-rd(file.path(old,'tables',paste0(tag,'_DEG.tsv')))
  j<-merge(tt,prior,by='gene_id',suffixes=c('_refit','_original'));stopifnot(nrow(j)==nrow(tt),nrow(tt)==nrow(prior))
  metrics<-c('logFC','logCPM','F','PValue','FDR');delta<-vapply(metrics,function(v)max(abs(j[[paste0(v,'_refit')]]-j[[paste0(v,'_original')]])),numeric(1))
  stopifnot(all(delta<1e-9))
  summaries[[tag]]<-data.frame(scope=scope,arm=arm,genes_tested=nrow(tt),n=6,n_control=3,n_case=3,cells=sum(md$cells),status='REPRODUCED_FROM_ORIGINAL_COUNTS',t(delta))
  q<-tt[tt$symbol%in%targets,];q$scope<-scope;q$arm<-arm;q$origin<-'Sep05_original_exact_refit';focal[[tag]]<-q
  offsets[[tag]]<-data.frame(scope=scope,arm=arm,mode='original_reset_after_filter',sample=samples,raw_library=colSums(counts),y$samples)
  wr(j[j$symbol_refit%in%targets,],paste0(tag,'_focal_refit_vs_original.tsv'))
  if(arm=='Full'&&scope%in%c('Core','ActCD4')){
   # One-factor diagnostic: identical counts/design/filter implementation, preserve full library offsets.
   y2<-DGEList(counts=counts,group=group,genes=genes);k2<-filterByExpr(y2,design=design);stopifnot(identical(keep,k2))
   y2<-y2[k2,,keep.lib.sizes=TRUE];y2<-calcNormFactors(y2,method='TMM');y2<-estimateDisp(y2,design,robust=TRUE);f2<-glmQLFit(y2,design,robust=TRUE)
   t2<-topTags(glmQLFTest(f2,coef=2),n=Inf,sort.by='PValue')$table;q2<-t2[t2$symbol%in%targets,];q2$scope<-scope;q2$arm<-arm;q2$origin<-'same_counts_full_library_offset_only';focal[[paste0(tag,'_offset')]]<-q2
   offsets[[paste0(tag,'_offset')]]<-data.frame(scope=scope,arm=arm,mode='diagnostic_full_library_preserved',sample=samples,raw_library=colSums(counts),y2$samples)
  }
  cat('EXACT_ORIGINAL_REFIT_PASS',tag,'genes',nrow(tt),'maxPdiff',delta['PValue'],'\n')
 }
}
wr(do.call(rbind,summaries),'original_edgeR_exact_refit_audit.tsv');wr(do.call(rbind,focal),'original_eight_genes_and_offset_diagnostic.tsv');wr(do.call(rbind,offsets),'normalization_offset_comparison.tsv')
capture.output(sessionInfo(),file=file.path(out,'results/refit_R_sessionInfo.txt'))
cat('ALL_ORIGINAL_EDGER_REFITS_PASS\n')
