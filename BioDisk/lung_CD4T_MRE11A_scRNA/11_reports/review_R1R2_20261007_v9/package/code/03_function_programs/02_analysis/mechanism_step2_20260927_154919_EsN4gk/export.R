args <- commandArgs(TRUE)
root <- args[1]; out <- args[2]
suppressPackageStartupMessages({library(Seurat); library(Matrix)})
p <- file.path(root,'06_CD4T_analysis/CD4_final/GSE288887_CD4_final.rds')
stopifnot(file.exists(p))
a <- readRDS(p); DefaultAssay(a) <- 'RNA'
if(inherits(a[['RNA']], 'Assay5')) {
  lyr <- Layers(a[['RNA']], search='^counts')
  if(length(lyr)!=1L || lyr!='counts') a <- JoinLayers(a,assay='RNA')
}
x <- GetAssayData(a,assay='RNA',layer='counts')
x <- as(x,'CsparseMatrix'); md <- a[[]][colnames(x),,drop=FALSE]
stopifnot(identical(rownames(md),colnames(x)),!anyDuplicated(colnames(x)),!anyDuplicated(rownames(x)))
stopifnot(all(is.finite(x@x)),all(x@x>=0),all(abs(x@x-round(x@x))<1e-6),all(Matrix::colSums(x)>0))
pick <- function(candidates) {
  found <- intersect(candidates,colnames(md))
  if(!length(found)) stop('Missing metadata: ',paste(candidates,collapse=','))
  found[1]
}
ss <- as.character(md[[pick(c('subject_id','subject','Subject','donor_id','patient_id'))]])
ss <- toupper(ss); hits <- regexpr('(CTL|CWP)[0-9]+',ss)
stopifnot(all(hits>0)); ss <- regmatches(ss,hits)
group <- tolower(as.character(md[[pick(c('group','Group','condition','Condition'))]]))
group <- ifelse(grepl('cwp',group),'CWP',ifelse(grepl('control|ctl',group),'Control',NA_character_))
ty <- as.character(md[[pick(c('cd4_label','CD4_anno','cd4_anno','Subtype','subtype'))]])
g227 <- c('CWP85','CWP86','CTL092','CTL094','CTL098','CTL103','CTL105','CTL106','CTL108','CTL111','CTL117')
g135 <- c('CTL068','CTL070');g288 <- c('CWP89','CWP94','CWP109','CWP127','CWP144','CWP149','CWP152')
series <- ifelse(ss%in%g227,'G227',ifelse(ss%in%g135,'G135',ifelse(ss%in%g288,'G288',NA_character_)))
stopifnot(!anyNA(series),!anyNA(group),!anyNA(ty),all(nzchar(ty)))
stopifnot(all((grepl('^CTL',ss)&group=='Control')|(grepl('^CWP',ss)&group=='CWP')))
m <- data.frame(cell=colnames(x),sample=ss,group=group,series=series,subtype=ty,core=ty!='CD4state',library=Matrix::colSums(x))
u <- unique(m[,c('sample','group','series')]);stopifnot(!anyDuplicated(u$sample))
write.table(m,file.path(out,'input/human_cells.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(data.frame(symbol=rownames(x)),file.path(out,'input/human_genes.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
writeMM(x,file.path(out,'input/human_counts.mtx'))
write.table(u,file.path(out,'tables/human_subject_design.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(as.data.frame(table(m$subtype,m$sample)),file.path(out,'tables/human_subtype_counts.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(data.frame(package=c('Seurat','Matrix','nichenetr','multinichenetr'),available=vapply(c('Seurat','Matrix','nichenetr','multinichenetr'),requireNamespace,logical(1),quietly=TRUE)),file.path(out,'tables/package_readiness.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(c(p,paste('cells',ncol(x),'genes',nrow(x))),file.path(out,'input/human_source.txt'))
cat('HUMAN_EXPORT_PASS\n')

