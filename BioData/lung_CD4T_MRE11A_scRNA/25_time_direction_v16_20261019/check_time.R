suppressPackageStartupMessages(library(edgeR))
P<-'/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA'
V<-file.path(P,'11_reports/death_priority_20261018_v31')
B<-file.path(P,'11_reports/three_cohort_four_axis_20261017_v16/inputs')
O<-'/mnt/BioData/bio/wangmingyuan/lung_CD4T_MRE11A_scRNA/25_time_direction_v16_20261019'
dir.create(O,recursive=TRUE,showWarnings=FALSE)
wr<-function(x,n)write.table(as.data.frame(x),file.path(O,n),sep='\t',quote=FALSE,row.names=FALSE,na='NA',fileEncoding='UTF-8')
a<-readRDS(file.path(V,'objects/time_all12.rds'));raw<-as.matrix(read.delim(gzfile(file.path(B,'GSE274372__TLike.counts.tsv.gz')),row.names=1,check.names=FALSE));md<-read.delim(file.path(B,'GSE274372__TLike.samples.tsv'),check.names=FALSE)
stopifnot(identical(colnames(a$y),as.character(a$meta$sample)),identical(colnames(raw),as.character(md$sample)),identical(as.character(a$meta$sample),as.character(md$sample)),all(as.numeric(as.character(a$meta$day))==md$day))
genes<-c('HIF1A','IL2RA','ARL6IP1');stopifnot(all(a$y$counts[genes,]==raw[genes,]))
lc<-cpm(a$y,log=TRUE,prior.count=.5);lc_old<-read.delim(file.path(V,'tables/time_all12_logCPM.tsv'),check.names=FALSE);rownames(lc_old)<-lc_old$gene;lc_old$gene<-NULL
stopifnot(isTRUE(all.equal(unname(as.matrix(lc_old[genes,])),unname(lc[genes,]),tolerance=1e-8)))
contrast<-c(0,-1,1,0);stopifnot(identical(colnames(a$design),c('D0','D7','D28','D56')))
r<-topTags(glmQLFTest(a$fit,contrast=contrast),n=Inf,sort.by='none')$table
old<-read.delim(file.path(V,'tables/time_all12_D28_vs_D7.tsv'));rownames(old)<-old$gene
stopifnot(isTRUE(all.equal(r[genes,'logFC'],old[genes,'logFC'],tolerance=1e-8)))
stopifnot(isTRUE(all.equal(r[genes,'PValue'],old[genes,'PValue'],tolerance=1e-8)),isTRUE(all.equal(r[genes,'FDR'],old[genes,'FDR'],tolerance=1e-8)))
rows<-list();means<-list();audit<-list()
libs<-a$y$samples$lib.size*a$y$samples$norm.factors
cp<-cpm(a$y,log=FALSE);bet<-a$fit$coefficients[genes,,drop=FALSE]
for(g in genes){
 i<-match(g,rownames(a$y));fitMu<-as.numeric(a$fit$fitted.values[i,]);offset<-if(is.matrix(a$fit$offset))as.numeric(a$fit$offset[i,]) else as.numeric(a$fit$offset)
 for(j in seq_len(ncol(a$y)))rows[[length(rows)+1]]<-data.frame(gene=g,sample=a$meta$sample[j],day=as.numeric(as.character(a$meta$day[j])),cells=a$meta$cells[j],full_library=a$meta$library[j],norm_factor=a$y$samples$norm.factors[j],effective_library=libs[j],count=a$y$counts[g,j],CPM=cp[g,j],logCPM=lc[g,j],fitted_count=fitMu[j],offset=offset[j])
 for(day in c(0,7,28,56)){
  ix<-as.numeric(as.character(a$meta$day))==day;coef<-bet[g,paste0('D',day)]
  means[[length(means)+1]]<-data.frame(gene=g,day=day,observed_mean_logCPM=mean(lc[g,ix]),observed_mean_CPM=mean(cp[g,ix]),model_log2CPM=coef/log(2)+log2(1e6),model_CPM=exp(coef)*1e6,count_sum=sum(a$y$counts[g,ix]),effective_library_sum=sum(libs[ix]))
 }
 z<-do.call(rbind,means);z<-z[z$gene==g,]
 delta<-z$model_log2CPM[z$day==28]-z$model_log2CPM[z$day==7];stopifnot(abs(delta-r[g,'logFC'])<1e-8)
 audit[[g]]<-data.frame(gene=g,contrast='D28-D7',reported_logFC=old[g,'logFC'],reproduced_logFC=r[g,'logFC'],coefficient_logFC=delta,P=old[g,'PValue'],FDR=old[g,'FDR'],observed_mean_logCPM_difference=z$observed_mean_logCPM[z$day==28]-z$observed_mean_logCPM[z$day==7],observed_mean_CPM_difference=z$observed_mean_CPM[z$day==28]-z$observed_mean_CPM[z$day==7],dispersion=if(length(a$fit$dispersion)==1L)a$fit$dispersion else a$fit$dispersion[i],prior_count=a$fit$prior.count)
}
wr(do.call(rbind,rows),'time_candidate_samples.tsv');wr(do.call(rbind,means),'time_candidate_group_estimates.tsv');wr(do.call(rbind,audit),'time_direction_audit.tsv');wr(data.frame(a$meta,a$design,check.names=FALSE),'verified_time_design.tsv')
writeLines(capture.output(sessionInfo()),file.path(O,'sessionInfo.txt'))
writeLines('success: sample IDs, day labels, counts, logCPM, design, contrast and reported coefficients verified; original fit reused, no refit',file.path(O,'complete.txt'))
print(as.data.frame(do.call(rbind,audit)));print(as.data.frame(do.call(rbind,means)))
