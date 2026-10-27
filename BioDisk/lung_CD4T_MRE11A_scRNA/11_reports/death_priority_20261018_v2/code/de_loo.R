suppressPackageStartupMessages(library(edgeR))
P<-'/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA';O<-file.path(P,'11_reports/death_priority_20261018_v2');B<-file.path(P,'11_reports/three_cohort_four_axis_20261017_v16')
fits<-readRDS(file.path(B,'inputs/edgeR_fits.rds'));cand<-read.delim(file.path(O,'tables/frozen_candidate_pairs.tsv'));genes<-unique(cand$gene);rows<-list();dg<-list()
for(sc in c('Core','Naive','ActCD4','Treg')){
 f<-fits[[paste0('joint_stratified__',sc)]]
 for(i in seq_len(nrow(f$meta))){
  m<-f$meta[-i,];m$condition<-factor(m$group,levels=c('Saline','Silica'));m$stratum<-factor(paste(m$dataset,m$day))
  design<-model.matrix(~stratum+condition,m);stopifnot(qr(design)$rank==ncol(design))
  y<-f$y[,-i,keep.lib.sizes=TRUE];y<-calcNormFactors(y);y<-estimateDisp(y,design,robust=TRUE);ff<-glmQLFit(y,design,robust=TRUE);q<-glmQLFTest(ff,coef=which(colnames(design)=='conditionSilica'));d<-topTags(q,n=Inf,sort.by='none')$table;d$gene<-rownames(d);d<-d[d$gene%in%genes,]
  d$scope<-sc;d$omitted_sample<-f$meta$sample[i];d$n<-9;d$full_family_n<-nrow(y);d$family<-paste(sc,'LOO',f$meta$sample[i],sep='|');rows[[length(rows)+1]]<-d
  dg[[length(dg)+1]]<-data.frame(scope=sc,omitted_sample=f$meta$sample[i],rank=qr(design)$rank,residual_df=nrow(design)-ncol(design),n_genes=nrow(y))
 }
 cat('LOO_COMPLETE',sc,'\n')
}
wr<-function(x,n)write.table(as.data.frame(x),file.path(O,'tables',n),sep='\t',quote=FALSE,row.names=FALSE,na='NA',fileEncoding='UTF-8')
wr(do.call(rbind,rows),'candidate_DE_leave_one_donor.tsv');wr(do.call(rbind,dg),'DE_leave_one_donor_diagnostics.tsv')
writeLines('success',file.path(O,'logs/DE_LOO_complete.txt'))
