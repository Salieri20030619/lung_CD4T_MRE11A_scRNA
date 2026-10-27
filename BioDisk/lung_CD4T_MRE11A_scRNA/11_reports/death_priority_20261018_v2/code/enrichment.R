suppressPackageStartupMessages({library(AnnotationDbi);library(org.Mm.eg.db);library(GO.db);library(fgsea)})
P<-'/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA';O<-file.path(P,'11_reports/death_priority_20261018_v2');old<-file.path(P,'02_analysis/GSE314943/stage7_20260905_234200_aiQYL0')
wr<-function(x,n)write.table(as.data.frame(x),file.path(O,'tables',n),sep='\t',quote=FALSE,row.names=FALSE,na='NA',fileEncoding='UTF-8')
canon<-function(x){x<-toupper(x);a<-c(MRE11A='MRE11',ATP5A1='ATP5F1A',ATP5B='ATP5F1B',ATP5C1='ATP5F1C',ATP5G1='ATP5MC1',ATP5O='ATP5PO',PARK2='PRKN');ii<-x%in%names(a);x[ii]<-a[x[ii]];x}
mp<-AnnotationDbi::select(org.Mm.eg.db,keys=keys(org.Mm.eg.db,keytype='SYMBOL'),keytype='SYMBOL',columns=c('ENTREZID','GOALL','ONTOLOGYALL'));mp$gene<-canon(mp$SYMBOL)
go<-unique(mp[!is.na(mp$GOALL)&mp$ONTOLOGYALL=='BP',c('gene','GOALL')]);gosets<-lapply(split(go$gene,go$GOALL),unique)
kn<-read.delim(file.path(old,'database/KEGG_names.tsv'),header=FALSE);kl<-read.delim(file.path(old,'database/KEGG_links.tsv'),header=FALSE);kl[[1]]<-sub('^mmu:','',kl[[1]]);kl[[2]]<-sub('^path:','',kl[[2]])
ids<-unique(mp[,c('gene','ENTREZID')]);ids<-ids[!is.na(ids$ENTREZID),];ks<-merge(kl,ids,by.x='V1',by.y='ENTREZID');ksets<-lapply(split(ks$gene,ks$V2),unique)
gn<-AnnotationDbi::select(GO.db,keys=names(gosets),keytype='GOID',columns='TERM')
desc<-c(setNames(gn$TERM,gn$GOID),setNames(kn[[2]],sub('^path:','',kn[[1]])))
pws<-list(GO_BP=gosets,KEGG=ksets);saveRDS(pws,file.path(O,'inputs/pathways_symbol.rds'))
wr(data.frame(pathway=names(desc),description=unname(desc)),'pathway_names.tsv');wr(unique(ids),'symbol_entrez_mapping.tsv')
res<-list();ora<-list();audit<-list()
for(sc in c('Core','Naive','ActCD4','Treg')){
 d<-read.delim(file.path(O,'inputs',paste0('joint_DE_',sc,'.tsv')),check.names=FALSE)
 stopifnot(!anyDuplicated(d$gene));d$rank<-sign(d$logFC)*sqrt(pmax(d$F,0));d<-d[is.finite(d$rank),];d<-d[order(-d$rank,d$gene),]
 wr(d[,c('gene','rank','logFC','PValue','FDR')],paste0('joint_rank_',sc,'.tsv'))
 for(db in names(pws)){
  universe<-intersect(d$gene,unique(unlist(pws[[db]])));pw<-lapply(pws[[db]],intersect,y=d$gene);pw<-pw[lengths(pw)>=15&lengths(pw)<=500]
  set.seed(20261018);f<-as.data.frame(fgseaMultilevel(pw,setNames(d$rank,d$gene),minSize=15,maxSize=500,eps=1e-10,nproc=2))
  f$scope<-sc;f$database<-db;f$family<-paste(sc,db,'joint_GSEA',sep='|');f$description<-unname(desc[f$pathway]);f$leadingEdge<-vapply(f$leadingEdge,paste,character(1),collapse=';');res[[length(res)+1]]<-f
  opw<-lapply(pws[[db]],intersect,y=universe);opw<-opw[lengths(opw)>=15&lengths(opw)<=500]
  for(direction in c('up','down')){
   selected<-intersect(d$gene[d$PValue<.05 & if(direction=='up')d$logFC>0 else d$logFC<0],universe)
   t<-data.frame(pathway=names(opw),size=lengths(opw),hits=vapply(opw,function(g)sum(g%in%selected),integer(1)))
   t$selected_n<-length(selected);t$universe_n<-length(universe);t$P<-phyper(t$hits-1,t$size,length(universe)-t$size,length(selected),lower.tail=FALSE);t$FDR<-p.adjust(t$P,'BH')
   t$genes<-vapply(opw,function(g)paste(intersect(g,selected),collapse=';'),character(1));t$scope<-sc;t$database<-db;t$direction<-direction;t$family<-paste(sc,db,direction,'joint_ORA',sep='|');t$description<-unname(desc[t$pathway]);ora[[length(ora)+1]]<-t
  }
  audit[[length(audit)+1]]<-data.frame(scope=sc,database=db,n_genes=nrow(d),universe=length(universe),GSEA_tests=nrow(f),ties=sum(duplicated(d$rank)))
  cat('COMPLETE',sc,db,nrow(f),'\n')
 }
}
wr(do.call(rbind,res),'joint_GSEA_all.tsv');wr(do.call(rbind,ora),'joint_ORA_all.tsv');wr(do.call(rbind,audit),'enrichment_audit.tsv')
writeLines(capture.output(sessionInfo()),file.path(O,'logs/enrichment_R_session.txt'));writeLines('success',file.path(O,'logs/enrichment_complete.txt'))
