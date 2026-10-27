P <- '/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA'
O <- file.path(P,'11_reports/manuscript_graphics_20261006_v5')
suppressPackageStartupMessages({library(SeuratObject);library(DOSE);library(enrichplot);library(jsonlite)})
wr <- function(x,n)write.table(x,file.path(O,'tables',n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
rd <- function(f)read.delim(f,check.names=FALSE,stringsAsFactors=FALSE)
cat('BEGIN_GSEA\n');flush.console()
z <- readRDS(file.path(P,'07_DEG/CD4_pseudobulk_v1/pathway_v1/GSE288887_CD4_pathway_results.rds'))
ref <- rd(file.path(P,'10_tables/CD4_final/pathway_v1/GO_KEGG_GSEA.tsv'))
qa <- list()
for(nm in names(z$GSEA)){
 g <- z$GSEA[[nm]]$object;rr<-g@result;tag<-strsplit(nm,'_')[[1]];db<-if(tag[2]=='GO')'GO_BP'else'KEGG';r<-ref[ref$analysis==tag[1]&ref$database==db,];r<-r[match(rr$ID,r$ID),]
 stopifnot(inherits(g,'gseaResult'),!anyNA(r$ID),nrow(rr)==nrow(r),all(is.finite(g@geneList)),!anyDuplicated(names(g@geneList)))
 delta<-sapply(c('NES','enrichmentScore','pvalue','p.adjust'),function(k)max(abs(rr[[k]]-r[[k]]),na.rm=TRUE));stopifnot(all(delta<1e-12))
 cat('GSEA',nm,'pathways',nrow(rr),'rank',length(g@geneList),'sets',length(g@geneSets),'deltas',delta,'\n')
 print(rr[order(-abs(rr$NES),rr$ID)[1:5],c('ID','Description','NES','pvalue','p.adjust')]);flush.console()
 qa[[nm]]<-list(class=class(g),pathways=nrow(rr),ranked_genes=length(g@geneList),gene_sets=length(g@geneSets),stat_delta=as.list(delta),params=g@params)
 saveRDS(g,file.path(O,'tables',paste0('Human_',nm,'_frozen_gsea.rds')))
}
write_json(qa,file.path(O,'provenance/GSEA_input_audit.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
rm(z,g);gc()
cat('BEGIN_CD4\n');flush.console()
cd <- readRDS(file.path(P,'06_CD4T_analysis/CD4_final/GSE288887_CD4_final.rds'))
cat('CD4_META',paste(colnames(cd@meta.data),collapse='|'),'\n');cat('CD4_REDUCTIONS',names(cd@reductions),'\n');print(table(cd@meta.data$CD4_subtype,useNA='ifany'));flush.console()
emb<-cd@reductions$umap@cell.embeddings;md<-cd@meta.data;stopifnot(identical(rownames(emb),rownames(md)))
md$cell<-rownames(md);md$UMAP_1<-emb[,1];md$UMAP_2<-emb[,2]
wr(md,'Human_CD4_UMAP_SERVER_ONLY.tsv');saveRDS(list(metadata=md,embedding=emb),file.path(O,'tables/Human_CD4_embedding_frozen.rds'))
rm(cd,md,emb);gc()
cat('BEGIN_BROAD_LUNG\n');flush.console()
b <- readRDS(file.path(P,'04_preprocessing/GSE288887_broad_clustered_scored.rds'))
cat('BROAD_META',paste(colnames(b@meta.data),collapse='|'),'\n');cat('BROAD_REDUCTIONS',names(b@reductions),'\n');print(table(b@meta.data$lineage_prelim,useNA='ifany'));flush.console()
emb<-b@reductions$umap@cell.embeddings;md<-b@meta.data;stopifnot(identical(rownames(emb),rownames(md)))
md$cell<-rownames(md);md$UMAP_1<-emb[,1];md$UMAP_2<-emb[,2]
wr(md,'Human_lung_UMAP_SERVER_ONLY.tsv');saveRDS(list(metadata=md,embedding=emb),file.path(O,'tables/Human_lung_embedding_frozen.rds'))
cat('INPUT_EXPORT_DONE',nrow(md),'\n');writeLines('complete',file.path(O,'provenance/inputs.done'));flush.console()
