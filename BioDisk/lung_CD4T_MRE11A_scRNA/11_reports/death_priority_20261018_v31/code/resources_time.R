suppressPackageStartupMessages({library(org.Mm.eg.db);library(GO.db);library(AnnotationDbi);library(edgeR)})
o <- '/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/death_priority_20261018_v31'
b <- '/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/three_cohort_four_axis_20261017_v16'
w <- function(x,n)write.table(as.data.frame(x),file=file.path(o,'tables',n),sep='\t',quote=FALSE,row.names=FALSE,na='NA',fileEncoding='UTF-8')
ids<-c('GO:0006915','GO:0070269','GO:0006914','GO:0000422','GO:0097193')
a<-AnnotationDbi::select(org.Mm.eg.db,keys=ids,keytype='GOALL',columns=c('SYMBOL','ENTREZID','EVIDENCEALL'));t<-AnnotationDbi::select(GO.db,keys=ids,keytype='GOID',columns='TERM');a<-merge(a,t,by.x='GOALL',by.y='GOID');a$gene<-toupper(a$SYMBOL);w(a,'GO_priority_members.tsv');w(AnnotationDbi::metadata(org.Mm.eg.db),'GO_resource_metadata.tsv')
readz<-function(){s<-'GSE274372__TLike';list(c=as.matrix(read.delim(gzfile(file.path(b,'inputs',paste0(s,'.counts.tsv.gz'))),row.names=1,check.names=FALSE)),m=read.delim(file.path(b,'inputs',paste0(s,'.samples.tsv')),stringsAsFactors=FALSE))}
for(v in c('all12','min10')){
 z<-readz();if(v=='min10'){ii<-z$m$cells>=10;z$c<-z$c[,ii,drop=FALSE];z$m<-z$m[ii,,drop=FALSE]}
 stopifnot(identical(colnames(z$c),z$m$sample));m<-z$m;m$day<-factor(m$day,levels=c(0,7,28,56));design<-model.matrix(~0+day,m);colnames(design)<-paste0('D',c(0,7,28,56));y<-DGEList(z$c,lib.size=m$library);keep<-filterByExpr(y,design=design);y<-y[keep,,keep.lib.sizes=TRUE];y<-calcNormFactors(y);y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE)
 cs<-list(D28_vs_D7=c(0,-1,1,0),D7_vs_D0=c(-1,1,0,0),D28_vs_D0=c(-1,0,1,0),D56_vs_D0=c(-1,0,0,1),D56_vs_D28=c(0,0,-1,1),time_omnibus=cbind(c(-1,1,0,0),c(-1,0,1,0),c(-1,0,0,1)))
 for(n in names(cs)){q<-glmQLFTest(fit,contrast=cs[[n]]);r<-topTags(q,n=Inf,sort.by='none')$table;r$gene<-rownames(r);r$contrast<-n;r$variant<-v;r$family<-paste(v,n,sep='|');w(r,paste0('time_',v,'_',n,'.tsv'));cat('TIME',v,n,nrow(r),sum(r$PValue<.05),'nominal\n',sep=' ')}
 lc<-cpm(y,log=TRUE,prior.count=.5);w(data.frame(gene=rownames(lc),lc,check.names=FALSE),paste0('time_',v,'_logCPM.tsv'));w(data.frame(m,design,check.names=FALSE),paste0('time_',v,'_design.tsv'));saveRDS(list(y=y,fit=fit,meta=m,design=design),file.path(o,'objects',paste0('time_',v,'.rds')))
}
writeLines(capture.output(sessionInfo()),file.path(o,'logs','R_session.txt'));writeLines('success',file.path(o,'logs','resources_time_complete.txt'))