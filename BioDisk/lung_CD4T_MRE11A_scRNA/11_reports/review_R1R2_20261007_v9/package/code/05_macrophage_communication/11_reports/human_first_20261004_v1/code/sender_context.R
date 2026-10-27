P<-'/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA';O<-file.path(P,'11_reports/human_first_20261004_v1')
writeLines('Supplementary sender disease context, specified before its computation. All reviewed MacAM, MonoMac and union MacWide; all genes, edgeR robust QL Base and Source on fixed >=20-cell donors, >=2/group, full library offsets/TMM, BH all tested genes per sender/model. This context does not select or rerank receiver NicheNet endpoints. Recorded after receiver computation; not a prospectively registered independent hypothesis.',file.path(O,'provenance/SENDER_CONTEXT_SPEC.txt'))
suppressPackageStartupMessages({library(edgeR);library(Matrix)})
ag<-readRDS(file.path(O,'results/stage2_aggregates.rds'));rows<-list();statuses<-list()
for(sen in names(ag$sender))for(model in c('Base','Source')){
 z<-ag$sender[[sen]];take<-z$sm$eligible20;md<-z$sm[take,];md$group<-factor(md$group,levels=c('Control','CWP'));md$series<-factor(md$series);X<-if(model=='Source')model.matrix(~group+series,md)else model.matrix(~group,md)
 st<-data.frame(sender=sen,model,n_control=sum(md$group=='Control'),n_case=sum(md$group=='CWP'),samples=paste(md$sample,collapse=';'),status='INSUFFICIENT')
 if(min(st$n_control,st$n_case)>=2&&qr(X)$rank==ncol(X)&&nrow(X)-ncol(X)>=2){y<-DGEList(z$pb[,take],lib.size=md$library);k<-filterByExpr(y,design=X);y<-y[k,,keep.lib.sizes=TRUE];y<-calcNormFactors(y);y<-estimateDisp(y,X,robust=TRUE);fit<-glmQLFit(y,X,robust=TRUE);tt<-topTags(glmQLFTest(fit,coef=2),n=Inf,sort.by='none')$table;rows[[length(rows)+1]]<-data.frame(gene=rownames(tt),tt,sender=sen,model,check.names=FALSE);st$status<-'FIT'}
 statuses[[length(statuses)+1]]<-st;cat('SENDER_CONTEXT',sen,model,'\n');flush.console()
}
con<-gzfile(file.path(O,'results/stage2_sender_disease_DE.tsv.gz'),'wt');write.table(do.call(rbind,rows),con,sep='\t',quote=FALSE,row.names=FALSE);close(con);write.table(do.call(rbind,statuses),file.path(O,'results/stage2_sender_DE_status.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
