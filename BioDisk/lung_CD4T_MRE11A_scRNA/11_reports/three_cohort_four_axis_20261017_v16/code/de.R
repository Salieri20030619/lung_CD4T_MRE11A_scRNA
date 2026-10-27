suppressPackageStartupMessages(library(edgeR))
o <- '/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/three_cohort_four_axis_20261017_v16'
rd <- function(ds,sc){s<-paste(ds,sc,sep='__');list(c=as.matrix(read.delim(gzfile(file.path(o,'inputs',paste0(s,'.counts.tsv.gz'))),row.names=1,check.names=FALSE)),m=read.delim(file.path(o,'inputs',paste0(s,'.samples.tsv')),stringsAsFactors=FALSE))}
w <- function(x,n)write.table(as.data.frame(x),file=file.path(o,'tables',n),sep='\t',row.names=FALSE,quote=FALSE,na='NA',fileEncoding='UTF-8')
results<-list();fits<-list();diags<-list()
fit <- function(z,id,sc,time=FALSE,subset=NULL){
 if(!is.null(subset)){z$c<-z$c[,subset,drop=FALSE];z$m<-z$m[subset,,drop=FALSE]}
 stopifnot(identical(colnames(z$c),z$m$sample));m<-z$m
 if(time){m$day<-factor(m$day,levels=c(0,7,28,56));d<-model.matrix(~day,m)} else {m$condition<-factor(m$group,levels=c('Saline','Silica'));m$stratum<-factor(paste(m$dataset,m$day,sep='_'));d<-if(nlevels(m$stratum)>1)model.matrix(~stratum+condition,m) else model.matrix(~condition,m)}
 stopifnot(qr(d)$rank==ncol(d),nrow(d)>ncol(d))
 y<-DGEList(z$c,lib.size=m$library);keep<-filterByExpr(y,design=d);y<-y[keep,,keep.lib.sizes=TRUE];y<-calcNormFactors(y);y<-estimateDisp(y,d,robust=TRUE);f<-glmQLFit(y,d,robust=TRUE)
 cs<-if(time)list(D7_vs_D0=2,D28_vs_D0=3,D56_vs_D0=4,time_omnibus=2:4) else list(Silica_vs_Saline=ncol(d))
 for(n in names(cs)){q<-glmQLFTest(f,coef=cs[[n]]);t<-topTags(q,n=Inf,sort.by='none')$table;t$gene<-rownames(t);t$model<-id;t$scope<-sc;t$contrast<-n;t$family<-paste(id,sc,n,sep='__');t$nominal<-t$PValue<.05;t$FDR_significant<-t$FDR<.05;t$effect_ge_0.5<-if('logFC'%in%names(t))abs(t$logFC)>=.5 else NA;results[[length(results)+1]]<<-t;w(t,paste0(t$family[1],'.tsv'));cat('FIT',id,sc,n,nrow(t),'nominal',sum(t$nominal),'FDR',sum(t$FDR_significant),'\n')}
 lc<-cpm(y,log=TRUE,prior.count=.5);w(data.frame(gene=rownames(lc),lc,check.names=FALSE),paste0(id,'__',sc,'__logCPM.tsv'))
 w(data.frame(m,d,check.names=FALSE),paste0(id,'__',sc,'__design.tsv'));fits[[paste(id,sc,sep='__')]]<<-list(y=y,fit=f,meta=m,design=d)
 diags[[length(diags)+1]]<<-data.frame(model=id,scope=sc,n_samples=nrow(d),design_rank=qr(d)$rank,residual_df=nrow(d)-ncol(d),tested_genes=nrow(y),common_dispersion=y$common.dispersion)
}
for(sc in c('Core','Naive','ActCD4','Treg')){
 a<-rd('GSE314943',sc);b<-rd('GSE183682',sc);fit(a,'314_disease',sc)
 g<-intersect(rownames(a$c),rownames(b$c));z<-list(c=cbind(a$c[g,,drop=FALSE],b$c[g,,drop=FALSE]),m=rbind(a$m,b$m));fit(z,'joint_stratified',sc)
}
z<-rd('GSE274372','TLike');fit(z,'274_time','TLike',TRUE);fit(z,'274_time_min10_sensitivity','TLike',TRUE,z$m$cells>=10)
fit(rd('GSE314943','AllMac'),'314_sender','AllMac')
# Keep full contrast-specific families separate; align unequal omnibus coefficient columns.
cols<-unique(unlist(lapply(results,names)));rr<-lapply(results,function(t){for(k in setdiff(cols,names(t)))t[[k]]<-NA;t[,cols]});all<-do.call(rbind,rr);w(all,'DE_all_models.tsv');w(do.call(rbind,diags),'model_diagnostics.tsv');saveRDS(fits,file.path(o,'inputs','edgeR_fits.rds'));writeLines(capture.output(sessionInfo()),file.path(o,'logs','R_session_DE.txt'));writeLines('success',file.path(o,'logs','DE_complete.txt'))