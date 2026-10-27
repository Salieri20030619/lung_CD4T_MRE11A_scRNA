args <- commandArgs(TRUE)
out <- args[1]
suppressPackageStartupMessages({library(Seurat);library(Matrix)})
root <- '/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA'
a <- readRDS(file.path(root,'06_CD4T_analysis/TNK_subcluster_v1/GSE288887_TNK_subclustered.rds'))
x <- LayerData(a,assay='RNA',layer='counts')
md <- a[[]][colnames(x),,drop=FALSE]
gs <- c('IL15','IL15RA','TNF','IL18','CXCL16','IL27','EBI3','IFNB1','IFNG')
symbols <- toupper(rownames(x))
v <- sapply(gs,function(g) if(any(symbols==g)) as.numeric(Matrix::colSums(x[symbols==g,,drop=FALSE])) else rep(NA_real_,ncol(x)))
d <- data.frame(cell=colnames(x),sample=toupper(as.character(md$subject_id)),group=as.character(md$group),library=Matrix::colSums(x),v,check.names=FALSE)
stopifnot(!anyDuplicated(d$cell),all(d$library>0),all(d$sample!=''))
write.table(d,file.path(out,'inputs/human_TNK_selected_counts.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
write.table(data.frame(gene=gs,present=gs%in%symbols),file.path(out,'inputs/human_TNK_gene_coverage.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'provenance/tnk_export_session.txt'))
cat('TNK_SELECTED_EXPORT_OK',nrow(d),'cells\n')
