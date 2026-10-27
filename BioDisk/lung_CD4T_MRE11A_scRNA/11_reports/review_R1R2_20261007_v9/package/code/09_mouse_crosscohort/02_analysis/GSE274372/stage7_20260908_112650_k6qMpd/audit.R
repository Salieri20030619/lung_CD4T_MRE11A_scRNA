root <- Sys.getenv('ROOT'); out <- file.path(Sys.getenv('RUN'),'tables')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
for (pkg in c('Biobase','SummarizedExperiment','monocle')) {
 if (requireNamespace(pkg, quietly=TRUE)) cat('AVAILABLE_PACKAGE=',pkg,'\n',sep='')
}
f <- file.path(root,'00_raw/GSE274372/GSE274372_cds_object.rds.gz')
stopifnot(file.exists(f))
f <- file.path(Sys.getenv('RUN'),'author_payload.bin')
x <- tryCatch(readRDS(f), error=function(e) NULL)
if (is.null(x)) {
 env <- new.env(parent=baseenv())
 loaded <- tryCatch(load(f,envir=env), error=function(e) stop(paste('RDS/RData recognition failed:',conditionMessage(e))))
 inventory <- vapply(loaded,function(n)paste(class(env[[n]]),collapse=','),character(1))
 write.table(data.frame(object=loaded,class=inventory),file.path(out,'author_loaded_objects.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
 usable <- loaded[vapply(loaded,function(n) {
  z <- env[[n]]
  isS4(z) && any(c('phenoData','colData') %in% names(attributes(z)))
 },logical(1))]
 if(length(usable)!=1L) stop('Expected exactly one CDS-like object; inspect author_loaded_objects.tsv')
 x <- env[[usable]]
 cat('AUTHOR_FORMAT=RData; OBJECT=',usable,'\n',sep='')
} else cat('AUTHOR_FORMAT=RDS\n')
info <- c(paste('class:',paste(class(x),collapse=',')), if(isS4(x)) paste('stored_attributes:',paste(names(attributes(x)),collapse=',')))
writeLines(info,file.path(out,'author_structure.txt'))
stored <- names(attributes(x))
if (isS4(x) && 'phenoData' %in% stored) {
 p <- attr(x,'phenoData',exact=TRUE)
 md <- attr(p,'data',exact=TRUE)
} else if (isS4(x) && 'colData' %in% stored) {
 md <- attr(x,'colData',exact=TRUE)
} else stop('No stored phenoData or colData attribute found')
if(is.null(md)) stop('Stored metadata is NULL')
# Conversion dispatch is on the metadata class, never on the monocle3 CDS.
md <- as.data.frame(md)
stopifnot(nrow(md)>0, !is.null(rownames(md)))
# Flatten only metadata columns, never the expression matrix.
for (k in names(md)) if (is.list(md[[k]])) md[[k]] <- vapply(md[[k]],paste,collapse=';',FUN.VALUE=character(1))
md <- cbind(author_cell_id=rownames(md),md)
con <- gzfile(file.path(out,'author_metadata.tsv.gz'),'wt')
write.table(md,con,sep='\t',quote=TRUE,row.names=FALSE,na='');close(con)
summary <- data.frame(column=names(md),unique_values=vapply(md,function(z)length(unique(z)),integer(1)))
write.table(summary,file.path(out,'author_columns.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
sink(file.path(out,'author_categories.txt'))
for(k in names(md)) if(length(unique(md[[k]]))<=100) {cat('\nCOLUMN:',k,'\n');print(table(md[[k]],useNA='ifany'))}
sink()
cat('AUTHOR_METADATA_PASS\n')
