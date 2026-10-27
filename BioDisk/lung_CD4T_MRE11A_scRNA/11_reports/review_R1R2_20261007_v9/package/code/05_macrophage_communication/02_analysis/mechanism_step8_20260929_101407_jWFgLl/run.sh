#!/usr/bin/env bash
set -Eeuo pipefail
ROOT="${ROOT:-/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA}"
PREV2="${PREV2:-$ROOT/02_analysis/mechanism_step2_20260927_154919_EsN4gk}"
PREV7="${PREV7:-$ROOT/02_analysis/mechanism_step7_20260929_100321_6hm4dZ}"
HUMAN_OBJECT="${HUMAN_OBJECT:-$ROOT/04_preprocessing/GSE288887_broad_clustered_scored.rds}"
export ROOT PREV2 PREV7 HUMAN_OBJECT
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-2}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
bash -n "$0"
mkdir -p "$ROOT/tmp" "$ROOT/logs" "$ROOT/02_analysis"
exec 9>"$ROOT/tmp/lung_mechanism_step8.lock"
flock -n 9 || { echo MECHANISM_STEP8_ALREADY_RUNNING; exit 1; }
OUT=$(mktemp -d "$ROOT/02_analysis/mechanism_step8_$(date +%Y%m%d_%H%M%S)_XXXXXX")
export OUT
trap 'rc=$?; if ((rc != 0)); then printf "LUNG_MECHANISM_STEP8_FAIL EXIT_CODE=%s OUTPUT=%s\n" "$rc" "$OUT"; fi' EXIT
mkdir -p "$OUT/tables" "$OUT/objects"
cp -- "$0" "$OUT/run.sh"
printf 'OUTPUT=%s\n' "$OUT"
cat > "$OUT/preflight.py" <<'PY_PREFLIGHT'
import os,json
from pathlib import Path
import pandas as pd
root=Path(os.environ['ROOT']);p7=Path(os.environ['PREV7']);p2=Path(os.environ['PREV2']);out=Path(os.environ['OUT']);obj=Path(os.environ['HUMAN_OBJECT'])
for p,stage in [(p7,'mechanism_step7_human_TF'),(p2,'mechanism_step2_human_and_driver_audit')]:
 s=json.loads((p/'status.json').read_text());assert s['status']=='PASS' and s['stage']==stage
inventory=pd.read_csv(p7/'tables/human_mac_object_inventory.tsv',sep='\t');row=inventory[inventory.path==str(obj)]
assert len(row)==1 and obj.is_file() and obj.stat().st_size==int(row.bytes.iloc[0]),'Selected object absent from verified inventory or size changed'
assert (p2/'input/human_cells.tsv').is_file()
(out/'selected_input.json').write_text(json.dumps({'path':str(obj),'bytes':obj.stat().st_size,'inventory':str(p7/'tables/human_mac_object_inventory.tsv'),'identity_status':'TO_BE_AUDITED'},indent=2))
print('INPUT_INVENTORY_CHECK=PASS',flush=True)
mem=Path('/proc/meminfo')
if mem.exists():
 for line in mem.read_text().splitlines():
  if line.startswith('MemAvailable:'):print(line,flush=True)

PY_PREFLIGHT
cat > "$OUT/audit.R" <<'R_AUDIT'
args <- commandArgs(TRUE);root <- args[1];p2 <- args[2];out <- args[3];input <- args[4]
suppressPackageStartupMessages({library(Seurat);library(Matrix)})
rd <- function(p)read.delim(p,check.names=FALSE,stringsAsFactors=FALSE)
wr <- function(d,n)write.table(d,file.path(out,'tables',paste0(n,'.tsv')),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
bind <- function(z)if(length(z))do.call(rbind,z) else data.frame(note=character())
cat('READ_OBJECT=',input,'\n',sep='');flush.console();a <- readRDS(input);stopifnot(inherits(a,'Seurat'),'RNA'%in%Assays(a))
md <- a[[]];cells <- rownames(md);stopifnot(!anyDuplicated(cells),length(cells)>6523)
anchor <- rd(file.path(p2,'input/human_cells.tsv'));stopifnot(!anyDuplicated(anchor$cell),all(c('sample','group','series')%in%names(anchor)))
ix <- match(cells,anchor$cell);shared <- !is.na(ix);cat('SHARED_CD4_ANCHOR_CELLS=',sum(shared),'\n',sep='');stopifnot(sum(shared)>0)
# Column inventory includes labels only, not free-text personal metadata.
cols <- list();label_cols <- character();sample_cols <- character()
for(n in names(md)){
 v <- md[[n]];if(!is.atomic(v)&&!is.factor(v))next
 u <- unique(as.character(v));u <- u[!is.na(u)&nzchar(u)]
 relevant <- grepl('type|anno|label|cluster|broad|ident|sample|subject|donor|patient|group|condition|series|batch',n,ignore.case=TRUE)
 cols[[length(cols)+1]] <- data.frame(column=n,class=paste(class(v),collapse=';'),unique_n=length(u),missing_n=sum(is.na(v)),examples=if(relevant&&length(u)<=200)paste(head(u,30),collapse=';') else '')
 if(grepl('subject|donor|patient|sample|orig.ident',n,ignore.case=TRUE)&&length(u)>1&&length(u)<=250)sample_cols <- c(sample_cols,n)
 if(grepl('type|anno|label|cluster|broad|snn_res',n,ignore.case=TRUE)&&!grepl('sample|subject|donor|patient|group|condition|series|batch',n,ignore.case=TRUE)&&length(u)>1&&length(u)<=150)label_cols <- c(label_cols,n)
}
wr(bind(cols),'metadata_columns');stopifnot(length(label_cols)>0)
# Rank annotation columns for marker audit only; no cell type is assigned by this ranking.
priority <- ifelse(grepl('type|anno|label',label_cols,ignore.case=TRUE),0,1);label_cols <- label_cols[order(priority,label_cols)]
label_cols <- head(label_cols,12);wr(data.frame(column=label_cols),'annotation_columns_audited')
# Learn one-key-to-one-subject maps only from exact frozen CD4 barcode matches.
map_stats <- list();maps <- list();mapping_edges <- list()
for(n in sample_cols){
 vals <- as.character(md[[n]]);valid <- shared&!is.na(vals)&nzchar(vals);dat <- data.frame(key=vals[valid],sample=anchor$sample[ix[valid]])
 keys <- split(dat$sample,dat$key);good <- keys[vapply(keys,function(v)length(unique(v))==1&&length(v)>=3,logical(1))]
 map <- if(length(good))setNames(vapply(good,function(v)v[1],character(1)),names(good)) else character()
 pred <- unname(map[vals]);bad_anchor <- sum(shared&!is.na(pred)&pred!=anchor$sample[ix],na.rm=TRUE)
 map_stats[[n]] <- data.frame(column=n,keys=length(keys),usable_keys=length(good),ambiguous_keys=sum(vapply(keys,function(v)length(unique(v))>1,logical(1))),mapped_cells=sum(!is.na(pred)),anchor_disagreements=bad_anchor)
 maps[[n]] <- pred
 if(length(good))mapping_edges[[n]] <- data.frame(column=n,key=names(good),sample=unname(map),anchor_cells=vapply(good,length,integer(1)))
}
ms <- bind(map_stats);wr(ms,'sample_mapping_candidates');wr(bind(mapping_edges),'sample_mapping_edges')
chosen <- NA_character_;donor <- rep(NA_character_,length(cells))
if(length(map_stats)){
 eligible <- ms[ms$anchor_disagreements==0&ms$ambiguous_keys==0&ms$usable_keys>=2,]
 if(nrow(eligible)){eligible <- eligible[order(-eligible$mapped_cells,eligible$column),];chosen <- eligible$column[1];donor <- maps[[chosen]]}
}
# Anchored cells keep their known identities even if the wider mapping is unavailable.
donor[shared] <- anchor$sample[ix[shared]]
subject <- unique(anchor[,c('sample','group','series')]);stopifnot(!anyDuplicated(subject$sample));di <- match(donor,subject$sample)
wr(data.frame(mapping_column=chosen,total_cells=length(cells),anchor_cells=sum(shared),mapped_cells=sum(!is.na(donor)),unknown_cells=sum(is.na(donor))),'sample_mapping_status')
# Read raw RNA layers separately; no re-clustering or layer integration.
panel <- strsplit('CD3D CD3E TRAC CD4 CD8A NKG7 GNLY MS4A1 CD79A EPCAM KRT19 SFTPC COL1A1 DCN LST1 TYROBP FCER1G CD68 MSR1 MERTK C1QA C1QB C1QC APOE MARCO FABP4 PPARG FCN1 S100A8 S100A9 VCAN CD1C FCER1A CLEC9A IL15 IL15RA TNF IL18 CXCL16 IL2RB IL2RG IL1R2 IFNG IFNGR1 IFNGR2 STAT1 STAT3 BATF',' ')[[1]]
mat <- matrix(0,length(panel),length(cells),dimnames=list(panel,cells));present <- setNames(rep(FALSE,length(panel)),panel);lib <- rep(NA_real_,length(cells));seen <- integer(length(cells));layerlog <- list()
layer_names <- if(inherits(a[['RNA']],'Assay5'))Layers(a[['RNA']],search='^counts') else 'counts';stopifnot(length(layer_names)>0)
for(ly in layer_names){
 x <- if(inherits(a[['RNA']],'Assay5'))LayerData(a,assay='RNA',layer=ly) else GetAssayData(a,assay='RNA',slot='counts')
 x <- as(x,'CsparseMatrix');idx <- match(colnames(x),cells);stopifnot(!anyNA(idx),!anyDuplicated(idx),all(seen[idx]==0),all(is.finite(x@x)),all(x@x>=0),all(abs(x@x-round(x@x))<1e-6))
 genes <- toupper(rownames(x));stopifnot(!is.null(rownames(x)));valid <- genes%in%panel
 if(any(valid)){
  pr <- sparseMatrix(i=match(genes[valid],panel),j=which(valid),x=1,dims=c(length(panel),nrow(x)))
  mat[,idx] <- as.matrix(pr%*%x);present[unique(genes[valid])] <- TRUE
 }
 lib[idx] <- as.numeric(Matrix::colSums(x));seen[idx] <- seen[idx]+1L
 layerlog[[ly]] <- data.frame(layer=ly,genes=nrow(x),cells=ncol(x));rm(x);gc(FALSE)
}
stopifnot(all(seen==1),all(is.finite(lib)),all(lib>0),all(colSums(mat)<=lib+1e-6))
wr(bind(layerlog),'RNA_layers');wr(data.frame(gene=panel,present=as.logical(present)),'panel_coverage')
cat('RAW_COUNT_AND_CELL_ALIGNMENT_PASS\n');flush.console()
# Missing panel genes stay NA in summaries; measured zero stays zero.
norm <- sweep(mat,2,lib,'/')*1e4;det <- mat>0;overall <- list();bysample <- list();counts <- list();anchor_tab <- list()
for(col in label_cols){
 labs <- as.character(md[[col]]);labs[is.na(labs)|!nzchar(labs)] <- 'Unlabeled';labels <- sort(unique(labs));counts[[col]] <- data.frame(column=col,label=labels,cells=as.integer(table(factor(labs,levels=labels))))
 for(label in labels){
  sel <- labs==label;n <- sum(sel);v <- rowMeans(norm[,sel,drop=FALSE]);pc <- rowMeans(det[,sel,drop=FALSE])*100;v[!present] <- pc[!present] <- NA_real_
  overall[[length(overall)+1]] <- data.frame(column=col,label,cells=n,gene=panel,mean_CP10k=v,pct=pc)
  for(ss in sort(unique(donor[sel&!is.na(donor)]))){
   ii <- sel&!is.na(donor)&donor==ss;nn <- sum(ii);j <- match(ss,subject$sample)
   v <- rowMeans(norm[,ii,drop=FALSE]);pc <- rowMeans(det[,ii,drop=FALSE])*100;raw <- rowSums(mat[,ii,drop=FALSE]);v[!present] <- pc[!present] <- raw[!present] <- NA_real_
   both <- if(all(present[c('IL15','IL15RA')]))mean(mat['IL15',ii]>0&mat['IL15RA',ii]>0)*100 else NA_real_
   bysample[[length(bysample)+1]] <- data.frame(column=col,label,sample=ss,group=subject$group[j],series=subject$series[j],cells=nn,eligible20=nn>=20,gene=panel,raw_sum=raw,mean_CP10k=v,pct=pc,IL15_IL15RA_pct=both)
  }
 }
 tab <- as.data.frame(table(label=labs[shared],CD4_subtype=anchor$subtype[ix[shared]]));tab$column <- col;anchor_tab[[col]] <- tab
 cat('ANNOTATION_AUDIT_DONE=',col,' labels=',length(labels),'\n',sep='');flush.console()
}
wr(bind(counts),'annotation_counts');wr(bind(overall),'annotation_marker_summary');wr(bind(bysample),'annotation_sample_expression');wr(bind(anchor_tab),'CD4_annotation_crosscheck')
# Retain a small server checkpoint for subsequent reviewed macrophage selection.
saveRDS(list(counts=mat,library=lib,present=present,metadata=data.frame(cell=cells,sample=donor,group=subject$group[di],series=subject$series[di],md[,label_cols,drop=FALSE],check.names=FALSE),source=input,mapping_column=chosen),file.path(out,'objects/marker_audit_checkpoint.rds'))
writeLines(c(input,paste('total_cells',length(cells)),paste('shared_CD4_cells',sum(shared)),paste('mapping_column',chosen)),file.path(out,'input_source.txt'))
writeLines(c('Human lung sender-input audit only. Read previously inventoried broad annotated Seurat object; retain all cells and original labels. No annotation is accepted solely from its name.',
'Subject mapping learned from exact shared CD4 barcodes to frozen step2 identities. Only unambiguous metadata keys with >=3 anchors used; mapping candidates/edges exported; unknown cells remain unassigned. Sample IDs, sources and disease groups never guessed from filename.',
'Raw RNA counts read layer-wise, each cell exactly once. Original full library totals are denominators. Panel genes absent from feature list reported NA; measured zeros retained.',
'Expression and IL15/IL15RA same-cell co-detection exported for every audited annotation label and mapped subject. Eligibility >=20 cells; no disease-DE, communication score or mechanistic validation claimed. Macrophage scope must be reviewed using markers and labels before the next statistical comparison.'),file.path(out,'methods.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'));cat('HUMAN_MAC_INPUT_AUDIT_PASS\n')

R_AUDIT
cat > "$OUT/finish.py" <<'PY_FINISH'
import os,json,hashlib,tarfile
from pathlib import Path
import pandas as pd
out=Path(os.environ['OUT']);root=Path(os.environ['ROOT'])
a=pd.read_csv(out/'tables/annotation_sample_expression.tsv',sep='\t');m=pd.read_csv(out/'tables/sample_mapping_status.tsv',sep='\t');coverage=pd.read_csv(out/'tables/panel_coverage.tsv',sep='\t')
assert not a.duplicated(['column','label','sample','gene']).any()
assert a.pct.dropna().between(0,100).all()
assert a.IL15_IL15RA_pct.dropna().between(0,100).all()
assert (m.mapped_cells+m.unknown_cells==m.total_cells).all()
status=dict(status='PASS',stage='mechanism_step8_human_mac_input_audit',output=str(out),input=os.environ['HUMAN_OBJECT'],annotation_columns=int(a['column'].nunique()),mapped_subjects=int(a['sample'].nunique()),total_cells=int(m.total_cells.iloc[0]),mapped_cells=int(m.mapped_cells.iloc[0]),panel_genes_present=int(coverage.present.sum()),next='review_macrophage_identity_and_matched_subjects_before_sender_statistics')
(out/'status.json').write_text(json.dumps(status,indent=2))
archive=root/'tmp/lung_mechanism_step8_return.tar.gz';part=archive.with_suffix('.tar.gz.part')
with tarfile.open(part,'w:gz') as t:
 for p in sorted(out.rglob('*')):
  if p.is_file() and not {'objects','__pycache__'}&set(p.relative_to(out).parts):t.add(p,arcname=str(p.relative_to(out)))
with tarfile.open(part) as t:
 for m in t:
  if m.isfile():
   f=t.extractfile(m)
   while f.read(1024*1024):pass
part.replace(archive)
print('RETURN='+str(archive),flush=True)
print('RETURN_BYTES='+str(archive.stat().st_size),flush=True)
print('RETURN_SHA256='+hashlib.sha256(archive.read_bytes()).hexdigest(),flush=True)
print('LUNG_MECHANISM_STEP8_PASS',flush=True)

PY_FINISH
conda run --no-capture-output -n batch_bench python -m py_compile "$OUT/preflight.py" "$OUT/finish.py"
conda run --no-capture-output -n R_4.5.2 Rscript -e 'invisible(parse(file=commandArgs(TRUE)[1]));cat("R_PARSE_PASS\n")' "$OUT/audit.R"
echo SYNTAX_CHECK=PASS
conda run --no-capture-output -n batch_bench python -u "$OUT/preflight.py"
conda run --no-capture-output -n R_4.5.2 Rscript "$OUT/audit.R" "$ROOT" "$PREV2" "$OUT" "$HUMAN_OBJECT"
conda run --no-capture-output -n batch_bench python -u "$OUT/finish.py"
