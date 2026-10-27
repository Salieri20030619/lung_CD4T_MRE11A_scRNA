suppressPackageStartupMessages({library(scTenifoldKnk);library(Matrix);library(jsonlite)})
O<-'/mnt/BioData/bio/wangmingyuan/lung_CD4T_MRE11A_scRNA/12_scTenifoldKnk_20261018';jobs<-fromJSON(file.path(O,'jobs.json'),simplifyVector=FALSE);rr<-list()
for(j in jobs){
 d<-file.path(O,'runs',j$id);WT<-as.matrix(readRDS(file.path(d,'WT.rds')))
 for(g in unlist(j$targets)){
  KO<-WT;KO[g,]<-0;set.seed(j$seed);MA<-scTenifoldNet::manifoldAlignment(WT,KO,d=2,nCores=4)
  saveRDS(MA,file.path(d,paste0(g,'_manifold.rds')))
  ng<-nrow(WT);distances<-sqrt(rowSums((MA[seq_len(ng),,drop=FALSE]-MA[ng+seq_len(ng),,drop=FALSE])^2));names(distances)<-sub('^X_','',rownames(MA)[seq_len(ng)])
  z<-read.delim(file.path(d,paste0(g,'.tsv')));dev<-max(abs(distances[z$gene]-z$distance));stopifnot(dev<1e-12)
  noise<-sqrt(.Machine$double.eps)*max(abs(MA));pp<-pchisq(z$FC,df=1,lower.tail=FALSE);pp[z$distance<=noise]<-1
  stopifnot(max(abs(pp-z$p.value))<1e-10)
  rr[[length(rr)+1]]<-data.frame(id=j$id,target=g,noiseLevel=noise,noise_genes=sum(z$distance<=noise),distance_reproduction_max_abs=dev)
 }
}
write.table(do.call(rbind,rr),file.path(O,'tables/numerical_noise_audit.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
cat('MANIFOLD REPRODUCED',length(rr),'\n')