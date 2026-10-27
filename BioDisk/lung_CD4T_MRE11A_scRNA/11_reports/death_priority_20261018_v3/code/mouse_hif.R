suppressPackageStartupMessages({library(fgsea);library(AnnotationDbi);library(org.Mm.eg.db)})
P<-'/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA';O<-file.path(P,'11_reports/death_priority_20261018_v3');V<-file.path(P,'11_reports/death_priority_20261018_v2')
wr<-function(x,n)write.table(as.data.frame(x),file.path(O,'tables',n),sep='\t',quote=FALSE,row.names=FALSE,na='NA',fileEncoding='UTF-8')
mm<-read.delim(file.path(V,'tables/module_members.tsv'));mm<-mm[tolower(as.character(mm$retained))=='true',];mods<-lapply(split(mm$gene,mm$module),unique)
pw<-readRDS(file.path(V,'inputs/pathways_symbol.rds'))$KEGG
res<-list();audit<-list();kg<-list()
for(ct in c('Hypoxia','Normoxia')){
 x<-read.delim(file.path(O,'external/GSE253090',paste0(ct,'_author_all_DE.tsv')),check.names=FALSE);x<-unique(x[,setdiff(names(x),'SYMBOL')]);stopifnot(!anyDuplicated(x$ems))
 mp<-AnnotationDbi::select(org.Mm.eg.db,keys=intersect(x$ems,keys(org.Mm.eg.db,keytype='ENSEMBL')),keytype='ENSEMBL',columns='SYMBOL');mp<-unique(mp[!is.na(mp$SYMBOL),]);nn<-table(mp$ENSEMBL);mp<-mp[mp$ENSEMBL%in%names(nn[nn==1]),];x$gene<-toupper(mp$SYMBOL[match(x$ems,mp$ENSEMBL)]);wr(x,paste0('GSE253090_',ct,'_author_DE_mapped.tsv'))
 y<-x[!is.na(x$gene)&is.finite(x$stat),];y<-y[order(-y$baseMean,y$ems),];y<-y[!duplicated(y$gene),];r<-setNames(y$stat,y$gene);r<-sort(r,decreasing=TRUE)
 for(db in c('Modules','KEGG')){
  sets<-if(db=='Modules')mods else pw
  set.seed(20261018);f<-as.data.frame(fgseaMultilevel(sets,r,minSize=if(db=='Modules')3 else 15,maxSize=500,eps=0,nPermSimple=10000,nproc=2))
  f$leadingEdge<-vapply(f$leadingEdge,paste,character(1),collapse=';');f$contrast<-paste0('GSE253090_',ct);f$database<-db;f$family<-paste(f$contrast,db,sep='|');res[[length(res)+1]]<-f
 }
 kg[[ct]]<-x[x$gene%in%c('HIF1A','IL2RA','ARL6IP1','PYCARD','GCLM','CASP3','CASP8','CASP9','CASP1','GSDMD','GSDME','BAX','BAK1','BCL2','BCL2L1','MCL1','CYCS','CYC1','HK2','LDHA','SLC2A1'),];kg[[ct]]$context<-ct
 audit[[ct]]<-data.frame(context=ct,unique_ensembl=nrow(x),ranked_unique_symbols=length(r),source_P_zero=sum(x$pvalue==0,na.rm=TRUE),data_type='author_DESeq2_all_gene_statistics; new downstream_GSEA')
}
wr(do.call(rbind,res),'mouse_HIF1A_perturbation_GSEA.tsv');wr(do.call(rbind,audit),'mouse_HIF1A_perturbation_audit.tsv');wr(do.call(rbind,kg),'mouse_HIF1A_key_genes.tsv')
# Refine only bounded/low-precision existing enrichment results and recalculate original complete families.
for(ct in c('thymus','spleen')){
 id<-paste0('GSE121883_',ct);f<-read.delim(file.path(O,'tables',paste0(id,'_KEGG_GSEA.tsv')));bad<-is.na(f$log2err)|f$pval<=1e-10
 if(any(bad)){
  d<-read.delim(file.path(O,'tables',paste0(id,'_DE.tsv')));r<-sort(setNames(sign(d$logFC)*sqrt(d$F),d$gene),decreasing=TRUE)
  set.seed(20261018);u<-as.data.frame(fgseaMultilevel(pw[f$pathway[bad]],r,minSize=15,maxSize=500,eps=0,nPermSimple=20000,nproc=2));u$leadingEdge<-vapply(u$leadingEdge,paste,character(1),collapse=';')
  ix<-match(u$pathway,f$pathway);for(col in c('pval','log2err','ES','NES','size','leadingEdge'))f[ix,col]<-u[[col]]
  f$padj<-p.adjust(f$pval,'BH');wr(f,paste0(id,'_KEGG_GSEA.tsv'))
 }
}
writeLines('success',file.path(O,'logs/mouse_HIF1A_complete.txt'))
