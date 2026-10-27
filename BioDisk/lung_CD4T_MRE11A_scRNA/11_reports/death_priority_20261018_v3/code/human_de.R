suppressPackageStartupMessages(library(edgeR))
P<-'/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA';O<-file.path(P,'11_reports/death_priority_20261018_v3');V<-file.path(P,'11_reports/death_priority_20261018_v2')
wr<-function(x,n)write.table(as.data.frame(x),file.path(O,'tables',n),sep='\t',quote=FALSE,row.names=FALSE,na='NA',fileEncoding='UTF-8')
x<-as.matrix(read.delim(gzfile(file.path(O,'tables/human_G227_counts.tsv.gz')),row.names=1,check.names=FALSE));m<-read.delim(file.path(O,'tables/human_G227_samples.tsv'))
all<-list();qc<-list()
for(sc in unique(m$subtype)){
 z<-m[m$subtype==sc,];z$group<-factor(z$group,levels=c('Control','CWP'));design<-model.matrix(~group,z);y<-DGEList(x[,z$column],lib.size=z$library);keep<-filterByExpr(y,design=design);y<-y[keep,,keep.lib.sizes=TRUE];y<-calcNormFactors(y);y<-estimateDisp(y,design,robust=TRUE);f<-glmQLFit(y,design,robust=TRUE);tt<-topTags(glmQLFTest(f,coef=2),n=Inf,sort.by='none')$table;tt$gene<-rownames(tt);tt$scope<-sc;tt$family<-paste0('G227|',sc,'|all_genes');wr(tt,paste0('human_G227_',sc,'_DE.tsv'));all[[sc]]<-tt[tt$gene%in%c('HIF1A','IL2RA','ARL6IP1','CD44','FOXP1','PRKCQ','TRPS1','PYCARD','GCLM','NFKB2'),]
 qc[[sc]]<-data.frame(scope=sc,n_control=sum(z$group=='Control'),n_case=sum(z$group=='CWP'),min_cells=min(z$cells),tested=nrow(tt))
}
wr(do.call(rbind,all),'human_G227_candidate_DE.tsv');wr(do.call(rbind,qc),'human_G227_DE_QC.tsv');writeLines('success',file.path(O,'logs/human_DE_complete.txt'))
