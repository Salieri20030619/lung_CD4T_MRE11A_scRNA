args <- commandArgs(TRUE); p8 <- args[1]; out <- args[2]
suppressPackageStartupMessages(library(Matrix))
ck <- readRDS(file.path(p8,'objects/marker_audit_checkpoint.rds'))
md <- ck$metadata; col <- 'RNA_snn_res.0.4'
stopifnot(col%in%names(md),!anyDuplicated(md$cell))
sel <- as.character(md[[col]])%in%c('3','6','16','22')
md <- md[sel,,drop=FALSE]; wanted <- md$cell
stopifnot(length(wanted)>0,all(is.finite(ck$library)),file.exists(ck$source))
oldlib <- ck$library[match(wanted,ck$metadata$cell)]
cat('HUMAN_SELECTED_CELLS=',length(wanted),'\n',sep='');flush.console()
suppressPackageStartupMessages(library(Seurat))
a <- readRDS(ck$source)
stopifnot(inherits(a,'Seurat'),'RNA'%in%Assays(a),all(wanted%in%colnames(a)))
layers <- if(inherits(a[['RNA']],'Assay5'))Layers(a[['RNA']],search='^counts') else 'counts'
pieces <- list(); seen <- integer(length(wanted)); audit <- list()
for(ly in layers){
 x <- if(inherits(a[['RNA']],'Assay5'))LayerData(a,assay='RNA',layer=ly) else GetAssayData(a,assay='RNA',slot='counts')
 keep <- which(colnames(x)%in%wanted)
 if(!length(keep))next
 x <- as(x[,keep,drop=FALSE],'CsparseMatrix'); ix <- match(colnames(x),wanted)
 stopifnot(!anyNA(ix),all(seen[ix]==0),all(is.finite(x@x)),all(x@x>=0),all(abs(x@x-round(x@x))<1e-6))
 stopifnot(isTRUE(all.equal(as.numeric(Matrix::colSums(x)),oldlib[ix],tolerance=1e-8)))
 genes <- toupper(rownames(x));stopifnot(!anyNA(genes),all(nzchar(genes)))
 rownames(x) <- genes;pieces[[ly]] <- x;seen[ix] <- seen[ix]+1L
 audit[[ly]] <- data.frame(layer=ly,cells=ncol(x),features=nrow(x))
}
stopifnot(all(seen==1));rm(a);gc(FALSE)
genes <- sort(unique(unlist(lapply(pieces,rownames),use.names=FALSE)))
# Aggregate duplicate symbols and place every cell in its original checkpoint order.
ii <- jj <- vv <- list()
for(k in seq_along(pieces)){
 x <- pieces[[k]];s <- summary(x)
 ii[[k]] <- match(rownames(x)[s$i],genes);jj[[k]] <- match(colnames(x)[s$j],wanted);vv[[k]] <- s$x
}
x <- sparseMatrix(i=unlist(ii,use.names=FALSE),j=unlist(jj,use.names=FALSE),x=unlist(vv,use.names=FALSE),dims=c(length(genes),length(wanted)),dimnames=list(genes,wanted))
stopifnot(isTRUE(all.equal(as.numeric(Matrix::colSums(x)),oldlib,tolerance=1e-8)))
md$parent <- as.character(md[[col]]);md$library <- oldlib
wr <- function(d,p)write.table(d,p,sep='\t',quote=FALSE,row.names=FALSE,na='NA')
wr(md,file.path(out,'input/human_mac_cells.tsv'))
wr(data.frame(gene=genes),file.path(out,'input/human_mac_genes.tsv'))
writeMM(x,file.path(out,'input/human_mac_counts.mtx'))
wr(do.call(rbind,audit),file.path(out,'tables/human_layer_audit.tsv'))
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
cat('ANALYSIS_LOGIC_PASS HUMAN_RAW_COUNTS_AND_MAPPING\n')

