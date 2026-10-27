#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA
mkdir -p "$ROOT/02_analysis/GSE274372" "$ROOT/tmp"
RUN=$(mktemp -d "$ROOT/02_analysis/GSE274372/stage7_$(date +%Y%m%d_%H%M%S)_XXXXXX")
export ROOT RUN OPENBLAS_NUM_THREADS=1
trap 'rc=$?; echo "GSE274372_STAGE7_FAIL EXIT_CODE=$rc"; exit "$rc"' ERR
cat > "$RUN/audit.R" <<'RS'
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
RS
cat > "$RUN/audit.py" <<'PY'
import os
from pathlib import Path
import numpy as np,pandas as pd,scanpy as sc
root=Path(os.environ['ROOT']);run=Path(os.environ['RUN']);out=run/'tables'
a=sc.read_h5ad(root/'02_analysis/GSE274372/stage5_20260907_224503_cUAf52/objects/T_refined.h5ad')
md=pd.read_csv(out/'author_metadata.tsv.gz',sep='\t',dtype=str,keep_default_na=False)
markers=['Cd4','Cd8a','Cd8b1','Zbtb7b','Runx3','Cd3d','Cd3e','Cd3g','Bcl11b','Themis','Skap1','Tcf7','Lef1','Ccr7','Sell','Il7r','Foxp3','Il2ra','Ctla4','Ikzf2','Icos','Nkg7','Ccl5','Prf1','Mre11a']
sym=a.var['symbol'].astype('string').fillna('').to_numpy()
e=a.obs.copy();e.index.name='cell_id'
for g in markers:
 idx=np.flatnonzero(sym==g)
 e[g]=np.asarray(a.layers['counts'][:,idx].sum(axis=1)).ravel() if len(idx) else np.nan
pd.DataFrame({'gene':markers,'feature_matches':[int((sym==g).sum()) for g in markers]}).to_csv(out/'marker_availability.tsv',sep='\t',index=False)
e.to_csv(out/'T_marker_evidence.tsv.gz',sep='\t')
e.loc[e['K20R08'].astype(str).eq('0')].to_csv(out/'Naive_marker_evidence.tsv.gz',sep='\t')
# Audit exact keys only. Repeated sci-RNA barcodes cannot establish sample identity.
rows=[]
for col in md.columns:
 vals=md[col].astype(str)
 for key,target in [('cell_id',pd.Series(a.obs_names)),('barcode',a.obs['barcode'].astype(str))]:
  matched=target.isin(set(vals))
  rows.append([col,key,int(matched.sum()),bool(vals.is_unique),bool(target.is_unique)])
pd.DataFrame(rows,columns=['author_column','local_key','matched_cells','author_unique','local_unique']).to_csv(out/'exact_key_audit.tsv',sep='\t',index=False)
ids=md['author_cell_id']; local=pd.Index(a.obs_names)
if ids.is_unique and local.isin(ids).any():
 joined=e.join(md.set_index('author_cell_id'),how='left',rsuffix='_author')
 joined.to_csv(out/'exact_id_join.tsv.gz',sep='\t')
(out/'notes.txt').write_text('Author labels are independent reference evidence, not ground truth. No fuzzy barcode matching or automatic relabeling. Barcode-only matches across samples are not sufficient. No expression imputation, new clustering, or time trend testing performed.\n')
print('GSE274372_STAGE7_ANALYSIS_PASS')
PY
cat > "$RUN/unpack.py" <<'PYCODE'
import os,gzip,shutil
from pathlib import Path
root=Path(os.environ['ROOT']);run=Path(os.environ['RUN'])
src=root/'00_raw/GSE274372/GSE274372_cds_object.rds.gz'
if not src.is_file(): raise FileNotFoundError(src)
current=src
for level in range(5):
 with current.open('rb') as f: magic=f.read(12)
 print('AUTHOR_LAYER=',level,' MAGIC_HEX=',magic.hex(),flush=True)
 if magic[:2]!=b'\x1f\x8b': break
 if level==4: raise RuntimeError('More than four gzip layers')
 dest=run/f'author_unpacked_{level}.bin'
 with gzip.open(current,'rb') as fi, dest.open('wb') as fo:
  shutil.copyfileobj(fi,fo,length=1024*1024)
 current=dest
# Only a run-local copy is passed to R. R handles any remaining supported compression.
shutil.copyfile(current,run/'author_payload.bin')
print('AUTHOR_UNPACK_PASS',flush=True)
PYCODE
bash -n "$0"
conda run -n R_4.5.2 Rscript -e 'parse(file=file.path(Sys.getenv("RUN"),"audit.R"));cat("R_SYNTAX_PASS\n")' >/dev/null
conda run -n batch_bench python -m py_compile "$RUN/audit.py" "$RUN/unpack.py"
echo SYNTAX_CHECK=PASS
echo "OUTPUT=$RUN"
conda run --no-capture-output -n batch_bench python -u "$RUN/unpack.py"
conda run --no-capture-output -n R_4.5.2 Rscript "$RUN/audit.R"
conda run --no-capture-output -n batch_bench python -u "$RUN/audit.py"
tar -czf "$ROOT/tmp/GSE274372_stage7_return.tar.gz" -C "$RUN" tables audit.R audit.py unpack.py
echo GSE274372_STAGE7_PASS
