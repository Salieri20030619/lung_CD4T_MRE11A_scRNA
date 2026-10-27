args <- commandArgs(TRUE)
stopifnot(length(args)==1)
root <- normalizePath(args[1],mustWork=TRUE)
files <- list.files(file.path(root,'code'),pattern='\\.[Rr]$',recursive=TRUE,full.names=TRUE)
rows <- lapply(files,function(f){
 message <- ''
 ok <- tryCatch({parse(file=f);TRUE},error=function(e){message <<- conditionMessage(e);FALSE})
 data.frame(file=substring(f,nchar(root)+2),check='R_parse_only',status=if(ok)'PASS' else 'FAIL',message=message)
})
out <- do.call(rbind,rows)
write.table(out,file.path(root,'provenance/r_parse_checks.tsv'),sep='\t',quote=TRUE,row.names=FALSE)
cat('R_PARSE_FILES',nrow(out),'FAILURES',sum(out$status!='PASS'),'\n')
if(any(out$status!='PASS'))print(out[out$status!='PASS',])
