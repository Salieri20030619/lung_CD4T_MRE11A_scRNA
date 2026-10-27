suppressPackageStartupMessages(library(edgeR))
O<-'/mnt/BioData/bio/wangmingyuan/lung_CD4T_MRE11A_scRNA/26_time_loo_v17_20261019'
dir.create(O,recursive=TRUE,showWarnings=FALSE)
wr<-function(x,n)write.table(as.data.frame(x),file.path(O,n),sep='\t',quote=FALSE,row.names=FALSE,na='NA',fileEncoding='UTF-8')
a<-readRDS('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/death_priority_20261018_v31/objects/time_all12.rds')
genes<-c('HIF1A','IL2RA');rows<-list();nf<-list();warn<-list()
original<-topTags(glmQLFTest(a$fit,contrast=c(0,-1,1,0)),n=Inf,sort.by='none')$table
for(k in 0:ncol(a$y)){
 keep<-if(k==0)seq_len(ncol(a$y)) else setdiff(seq_len(ncol(a$y)),k)
 deleted<-if(k==0)'NONE' else as.character(a$meta$sample[k])
 m<-a$meta[keep,,drop=FALSE];m$day<-factor(m$day,levels=c(0,7,28,56))
 design<-model.matrix(~0+day,m);colnames(design)<-c('D0','D7','D28','D56')
 y<-DGEList(a$y$counts[,keep,drop=FALSE],lib.size=a$y$samples$lib.size[keep])
 result<-withCallingHandlers({
 y<-calcNormFactors(y);y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE)
 topTags(glmQLFTest(fit,contrast=c(0,-1,1,0)),n=Inf,sort.by='none')$table
 },warning=function(w){warn[[length(warn)+1]]<<-data.frame(deleted=deleted,message=conditionMessage(w));invokeRestart('muffleWarning')})
 if(k==0)stopifnot(isTRUE(all.equal(result$logFC,original$logFC,tolerance=1e-8)),isTRUE(all.equal(result$PValue,original$PValue,tolerance=1e-8)))
 for(g in genes){
  z<-result[g,];rows[[length(rows)+1]]<-data.frame(gene=g,deleted=deleted,deleted_day=if(k==0)NA else as.numeric(as.character(a$meta$day[k])),n_total=ncol(y),n_D7=sum(m$day==7),n_D28=sum(m$day==28),gene_family=nrow(y),logFC=z$logFC,P=z$PValue,FDR=z$FDR,full_logFC=original[g,'logFC'],delta_logFC=z$logFC-original[g,'logFC'],same_direction=sign(z$logFC)==sign(original[g,'logFC']),family=paste0('fixed_original_genes|D28-D7|omit_',deleted))
 }
 nf[[length(nf)+1]]<-data.frame(deleted=deleted,sample=as.character(m$sample),day=as.numeric(as.character(m$day)),library=y$samples$lib.size,norm_factor=y$samples$norm.factors,effective_library=y$samples$lib.size*y$samples$norm.factors)
 cat('Completed',deleted,'\n')
}
wr(do.call(rbind,rows),'time_LOO_candidates.tsv');wr(do.call(rbind,nf),'time_LOO_normalization.tsv')
wr(if(length(warn))do.call(rbind,warn) else data.frame(deleted=character(),message=character()),'warnings.tsv')
d<-as.numeric(as.character(a$meta$day))==28
eff<-a$y$samples$lib.size*a$y$samples$norm.factors
prop<-data.frame(sample=a$meta$sample[d],cells=a$meta$cells[d],cell_share=a$meta$cells[d]/sum(a$meta$cells[d]),effective_library=eff[d],effective_library_share=eff[d]/sum(eff[d]),HIF1A_count=a$y$counts['HIF1A',d],HIF1A_count_share=a$y$counts['HIF1A',d]/sum(a$y$counts['HIF1A',d]),IL2RA_count=a$y$counts['IL2RA',d],IL2RA_count_share=a$y$counts['IL2RA',d]/sum(a$y$counts['IL2RA',d]))
wr(prop,'D28_contributions.tsv')
writeLines(rownames(a$y),file.path(O,'fixed_gene_universe.txt'))
writeLines(capture.output(sessionInfo()),file.path(O,'sessionInfo.txt'))
writeLines('SUCCESS: full refit reproduced original FC/P; 12 leave-one-out refits complete; original gene family fixed; TMM and dispersions refitted; exploratory sensitivity only, no reselection',file.path(O,'complete.txt'))
