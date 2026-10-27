suppressPackageStartupMessages(library(fgsea))
P<-'/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA';O<-file.path(P,'11_reports/death_priority_20261018_v3');V<-file.path(P,'11_reports/death_priority_20261018_v2')
wr<-function(x,n)write.table(as.data.frame(x),file.path(O,'tables',n),sep='\t',quote=FALSE,row.names=FALSE,na='NA',fileEncoding='UTF-8')
edges<-read.delim(gzfile(file.path(V,'tables/joint_network_edge_stability.tsv.gz')));cand<-read.delim(file.path(V,'tables/frozen_candidate_pairs.tsv'));mods<-read.delim(file.path(V,'tables/module_members.tsv'));cov<-read.delim(file.path(V,'tables/joint_network_TF_coverage.tsv'))
excluded<-unique(c(cand$gene,mods$gene));all<-list();cover<-list();signatures<-list()
for(arm in c('All','BG1')){
 ee<-edges[edges$arm==arm & edges$seed_support>=2 & !edges$target%in%excluded & edges$tf!=edges$target,]
 ee<-ee[ee$tf%in%cov$tf[cov$arm==arm & cov$status=='SCORED'],]
 sets<-lapply(split(ee$target,ee$tf),unique);sets<-sets[lengths(sets)>=5]
 for(tf in names(sets))signatures[[length(signatures)+1]]<-data.frame(arm=arm,tf=tf,gene=sets[[tf]])
 for(id in c('GSE121883_spleen','GSE121883_thymus','GSE324038_HumanPanT_Hypoxia','GSE253090_Hypoxia','GSE253090_Normoxia')){
  if(startsWith(id,'GSE253090')){
   x<-read.delim(file.path(O,'tables',paste0(id,'_author_DE_mapped.tsv')));x<-x[!is.na(x$gene)&is.finite(x$stat),];x<-x[order(-x$baseMean,x$ems),];x<-x[!duplicated(x$gene),];r<-setNames(x$stat,x$gene)
  }else{
   x<-read.delim(file.path(O,'tables',paste0(id,'_DE.tsv')));r<-setNames(sign(x$logFC)*sqrt(x$F),x$gene)
  }
  r<-sort(r,decreasing=TRUE);ss<-lapply(sets,intersect,y=names(r));cc<-data.frame(arm=arm,contrast=id,tf=names(ss),original_targets=lengths(sets),tested_targets=lengths(ss),eligible=lengths(ss)>=5);cover[[length(cover)+1]]<-cc
  ss<-ss[lengths(ss)>=5]
  set.seed(20261018);f<-as.data.frame(fgseaMultilevel(ss,r,minSize=5,maxSize=length(r),eps=0,nPermSimple=10000,nproc=2));f$leadingEdge<-vapply(f$leadingEdge,paste,character(1),collapse=';');f$arm<-arm;f$contrast<-id;f$family<-paste(id,arm,'joint_regulon_GSEA',sep='|');all[[length(all)+1]]<-f
  cat('DONE',arm,id,nrow(f),'\n')
 }
}
wr(do.call(rbind,all),'external_joint_regulon_GSEA.tsv');wr(do.call(rbind,cover),'external_joint_regulon_coverage.tsv');wr(do.call(rbind,signatures),'external_joint_regulon_members.tsv')
writeLines('success',file.path(O,'logs/external_regulon_complete.txt'))
