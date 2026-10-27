suppressPackageStartupMessages(library(fgsea))
P<-'/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA';O<-file.path(P,'11_reports/death_priority_20261018_v3');V<-file.path(P,'11_reports/death_priority_20261018_v2')
wr<-function(x,n)write.table(as.data.frame(x),file.path(O,'tables',n),sep='\t',quote=FALSE,row.names=FALSE,na='NA',fileEncoding='UTF-8')
if(file.exists(file.path(O,'logs/GSEA_precision_complete.txt')))quit(status=0)
mm<-read.delim(file.path(V,'tables/module_members.tsv'));kept<-mm[tolower(as.character(mm$retained))=='true',];mods<-lapply(split(kept$gene,kept$module),unique)
kegg<-readRDS(file.path(V,'inputs/pathways_symbol.rds'))$KEGG
rm<-read.delim(file.path(O,'tables/external_joint_regulon_members.tsv'))
getr<-function(id){
 if(startsWith(id,'GSE253090')){
  x<-read.delim(file.path(O,'tables',paste0(id,'_author_DE_mapped.tsv')));x<-x[!is.na(x$gene)&is.finite(x$stat),];x<-x[order(-x$baseMean,x$ems),];x<-x[!duplicated(x$gene),];r<-setNames(x$stat,x$gene)
 }else{x<-read.delim(file.path(O,'tables',paste0(id,'_DE.tsv')));r<-setNames(sign(x$logFC)*sqrt(x$F),x$gene)}
 sort(r,decreasing=TRUE)
}
files<-c('external_module_GSEA.tsv','mouse_HIF1A_perturbation_GSEA.tsv','external_joint_regulon_GSEA.tsv','GSE121883_spleen_KEGG_GSEA.tsv','GSE121883_thymus_KEGG_GSEA.tsv')
plan<-list()
for(name in files){
 path<-file.path(O,'tables',name);backup<-sub('.tsv$','_before_precision.tsv',path)
 if(!file.exists(backup))stopifnot(file.copy(path,backup))
 x<-read.delim(backup);q<-x[is.finite(x$pval)&x$pval>=.04&x$pval<=.06,];q$file<-name;plan[[name]]<-q[,c('file','family','contrast','pathway','pval')]
}
wr(do.call(rbind,plan),'GSEA_precision_frozen_plan.tsv')
writeLines('One fixed refinement for ALL new GSEA entries with first-pass P in [0.04,0.06]; nPermSimple100000, sampleSize1001, eps0, seed20261018. No rerun based on significance. Replace estimates and recalculate complete original family BH.',file.path(O,'logs/GSEA_precision_protocol.txt'))
audit<-list()
for(name in files){
 path<-file.path(O,'tables',name);backup<-sub('.tsv$','_before_precision.tsv',path);x<-read.delim(backup)
 for(fam in unique(x$family)){
  ii<-which(x$family==fam & is.finite(x$pval)&x$pval>=.04&x$pval<=.06)
  if(!length(ii))next
  id<-x$contrast[ii[1]];r<-getr(id);minsize<-3;maxsize<-500
  if(name=='external_joint_regulon_GSEA.tsv'){
   z<-rm[rm$arm==x$arm[ii[1]],];sets<-lapply(split(z$gene,z$tf),unique);minsize<-5;maxsize<-length(r)
  }else if(grepl('KEGG',name,fixed=TRUE)||(name=='mouse_HIF1A_perturbation_GSEA.tsv'&&x$database[ii[1]]=='KEGG')){
   sets<-kegg;minsize<-15
  }else sets<-mods
  set.seed(20261018);u<-as.data.frame(fgseaMultilevel(sets[x$pathway[ii]],r,minSize=minsize,maxSize=maxsize,eps=0,nPermSimple=100000,sampleSize=1001,nproc=2))
  stopifnot(setequal(u$pathway,x$pathway[ii]),!anyNA(u$pval))
  jj<-match(u$pathway,x$pathway);# path names repeat across contexts; restrict the match to this family
  jj<-ii[match(u$pathway,x$pathway[ii])]
  stopifnot(identical(as.integer(u$size),as.integer(x$size[jj])))
  audit[[length(audit)+1]]<-data.frame(file=name,family=fam,pathway=u$pathway,P_before=x$pval[jj],P_after=u$pval,NES_before=x$NES[jj],NES_after=u$NES,log2err_after=u$log2err)
  u$leadingEdge<-vapply(u$leadingEdge,paste,character(1),collapse=';')
  for(col in c('pval','log2err','ES','NES','size','leadingEdge'))x[jj,col]<-u[[col]]
  cat('REFINED',name,fam,length(ii),'\n')
 }
 for(fam in unique(x$family)){ii<-which(x$family==fam);x$padj[ii]<-p.adjust(x$pval[ii],'BH')}
 wr(x,name)
}
wr(do.call(rbind,audit),'GSEA_precision_audit.tsv');writeLines('success',file.path(O,'logs/GSEA_precision_complete.txt'))
