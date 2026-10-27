#!/usr/bin/env bash
set -Eeuo pipefail
ROOT="${ROOT:-/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA}"
PREV1="${PREV1:-$ROOT/02_analysis/mechanism_step1_20260927_152218_s4V5eF}"
PREV3="${PREV3:-$ROOT/02_analysis/mechanism_step3_20260928_103057_lIRL38}"
PREV4="${PREV4:-$ROOT/02_analysis/mechanism_step4_20260928_150337_J58beC}"
CACHE="${CACHE:-$ROOT/00_raw/nichenet_v2_mouse}"
export ROOT PREV1 PREV3 PREV4 CACHE
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-2}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
bash -n "$0"
mkdir -p "$ROOT/tmp" "$ROOT/logs" "$ROOT/02_analysis"
exec 9>"$ROOT/tmp/lung_mechanism_step5.lock"
flock -n 9 || { echo MECHANISM_STEP5_ALREADY_RUNNING; exit 1; }
OUT=$(mktemp -d "$ROOT/02_analysis/mechanism_step5_$(date +%Y%m%d_%H%M%S)_XXXXXX")
export OUT
trap 'rc=$?; if ((rc != 0)); then printf "LUNG_MECHANISM_STEP5_FAIL EXIT_CODE=%s OUTPUT=%s\n" "$rc" "$OUT"; fi' EXIT
mkdir -p "$OUT/tables"
cp -- "$0" "$OUT/run.sh"
printf 'OUTPUT=%s\n' "$OUT"
cat > "$OUT/preflight.py" <<'PY_CHECK'
import os,json,hashlib
from pathlib import Path
for var,stage in [('PREV1','mechanism_step1'),('PREV3','mechanism_step3_nichenet'),('PREV4','mechanism_step4_receptor_audit')]:
 p=Path(os.environ[var]);s=json.loads((p/'status.json').read_text());assert s['status']=='PASS' and s['stage']==stage
p3=Path(os.environ['PREV3']);cache=Path(os.environ['CACHE'])
manifest=json.loads((p3/'resource_manifest.json').read_text())
record=next(v for v in manifest if v['file']=='ligand_target_matrix_nsga2r_final_mouse.rds')
p=cache/record['file'];assert p.is_file() and p.stat().st_size==record['bytes']
algorithm,expected=record['checksum'].split(':',1);h=hashlib.new(algorithm)
with p.open('rb') as f:
 for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
assert h.hexdigest()==expected,'NicheNet cached prior checksum mismatch'
print('INPUT_AND_PRIOR_CHECK=PASS',flush=True)

PY_CHECK
cat > "$OUT/analyze.R" <<'R_ANALYZE'
args <- commandArgs(TRUE);p1 <- args[1];p3 <- args[2];out <- args[3];cache <- args[4]
suppressPackageStartupMessages({library(Matrix);library(nichenetr)})
set.seed(314943)
rd <- function(p)read.delim(p,check.names=FALSE,stringsAsFactors=FALSE)
wr <- function(x,n)write.table(x,file.path(out,'tables',n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
canon <- function(x){x <- toupper(x);x[x=='MRE11A'] <- 'MRE11';x}
targets <- c('IL1R2','PENK','TTC39C','FOXP1','TNFRSF9','LAG3','MRE11')
g <- rd(file.path(p1,'comm_input/CD4_genes.tsv'));m <- rd(file.path(p1,'comm_input/CD4_cells.tsv'))
g3 <- rd(file.path(p3,'tables/CD4_gene_mapping.tsv'))
stopifnot(identical(g$canonical,g3$canonical),!anyNA(g3$mapped),!anyDuplicated(m$cell),nrow(m)==2311)
con <- gzfile(file.path(p1,'comm_input/CD4_counts.mtx.gz'),'rt');x <- as(readMM(con),'CsparseMatrix');close(con)
stopifnot(nrow(x)==nrow(g),ncol(x)==nrow(m),all(x@x>=0),all(is.finite(x@x)),all(abs(x@x-round(x@x))<1e-6),all(m$library>0))
bgp <- strsplit('SFTPA1 SFTPA2 SFTPB SFTPC SFTPD SCGB1A1 SCGB3A1 SCGB3A2 ABCA3 SLC34A2 NAPSA LYZ LYZ1 LYZ2 S100A8 S100A9 CSF1R LST1 TYROBP FCER1G',' ')[[1]]
bgi <- g$canonical%in%bgp|grepl('^(IGH[VDJMAEG]|IGK[VJC]|IGL[VJC])',g$canonical)
m$bg <- as.numeric(Matrix::colSums(x[bgi,,drop=FALSE]))/m$library
u <- sort(unique(g3$mapped));pr <- sparseMatrix(i=match(g3$mapped,u),j=seq_len(nrow(g3)),x=1,dims=c(length(u),nrow(g3)))
x <- pr%*%x;rownames(x) <- u;colnames(x) <- m$cell
prior <- file.path(cache,'ligand_target_matrix_nsga2r_final_mouse.rds');stopifnot(file.exists(prior))
lt <- readRDS(prior);stopifnot(!anyDuplicated(rownames(lt)),!anyDuplicated(colnames(lt)))
lr <- rd(file.path(p3,'tables/supported_LR_edges.tsv'))
old <- rd(file.path(p3,'tables/ligand_activities.tsv'));old_status <- rd(file.path(p3,'tables/geneset_status.tsv'))
activities <- list();statuses <- list();coverage <- list();sets <- list();reproduction <- list()
make_sets <- function(de,background,test,direction,exclude=character()) {
 bg <- setdiff(background,exclude)
 p <- if(test=='FDR05')de$FDR else de$PValue
 keep <- is.finite(p)&p<.05&(if(direction=='Up')de$logFC>=.5 else de$logFC<= -.5)
 list(bg=bg,gs=intersect(de$gene[keep],bg))
}
for(receiver in c('Core','ActCD4'))for(filt in c('All','BG1')) {
 take <- if(receiver=='Core')rep(TRUE,nrow(m)) else m$subtype==receiver
 if(filt=='BG1')take <- take&m$bg<=.01
 mm <- m[take,,drop=FALSE];xx <- x[,take,drop=FALSE]
 ss <- sort(unique(m$sample));n <- as.numeric(table(factor(mm$sample,levels=ss)))
 agg <- sparseMatrix(i=seq_len(nrow(mm)),j=match(mm$sample,ss),x=1,dims=c(nrow(mm),length(ss)))
 pct <- sweep(as.matrix((xx>0)%*%agg),2,pmax(n,1),'/')
 stopifnot(all(n>=20),sum(grepl('^Saline_',ss))==3,sum(grepl('^Silica_',ss))==3)
 de <- rd(file.path(p3,'tables',paste0('DE_CD4_',receiver,'_',filt,'.tsv')))
 potential <- sort(intersect(unique(lr$from[lr$receiver==receiver & lr$filter==filt & lr$threshold=='Strict10']),colnames(lt)))
 for(cut in c(.1,.05)) {
   gate <- if(cut==.1)'Pct10' else 'Pct5'
   ec <- rownames(pct)[rowSums(pct[,grepl('^Saline_',ss),drop=FALSE]>=cut)>=2]
   ed <- rownames(pct)[rowSums(pct[,grepl('^Silica_',ss),drop=FALSE]>=cut)>=2]
   background <- Reduce(intersect,list(de$gene,union(ec,ed),rownames(lt)))
   coverage[[length(coverage)+1]] <- data.frame(receiver,filter=filt,gate,target=targets,in_background=targets%in%canon(background),DE_tested=targets%in%canon(de$gene))
   for(test in c('FDR05','Nominal05'))for(direction in c('Up','Down'))for(omit in c('Full','Drop7')) {
     exclude <- if(omit=='Drop7')rownames(lt)[canon(rownames(lt))%in%targets] else character()
     z <- make_sets(de,background,test,direction,exclude);gs <- z$gs;bg <- z$bg
     st <- data.frame(receiver,filter=filt,gate,test,direction,omit,background=length(bg),geneset=length(gs),ligands=length(potential),status='PASS')
     if(length(bg)<500)st$status <- 'SKIP_SMALL_BACKGROUND'
     if(length(gs)<10)st$status <- 'SKIP_GENESET_LT10'
     if(length(setdiff(bg,gs))<10)st$status <- 'SKIP_NEGATIVE_GENES_LT10'
     if(length(potential)<2)st$status <- 'SKIP_LIGANDS_LT2'
     statuses[[length(statuses)+1]] <- st
     if(length(gs))sets[[length(sets)+1]] <- data.frame(receiver,filter=filt,gate,test,direction,omit,gene=gs)
     if(gate=='Pct10'&&omit=='Full') {
       os <- old_status[old_status$receiver==receiver & old_status$filter==filt & old_status$test==test & old_status$direction==direction,,drop=FALSE]
       stopifnot(nrow(os)==1,st$background==os$background,st$geneset==os$geneset,st$ligands==os$ligands,st$status==os$status)
     }
     if(st$status!='PASS')next
     # Reuse verified original 10% scores rather than recomputing them.
     if(gate=='Pct10'&&omit=='Full') {
       act <- old[old$receiver==receiver & old$filter==filt & old$test==test & old$direction==direction,,drop=FALSE]
       stopifnot(nrow(act)==length(potential),setequal(act$test_ligand,potential))
       act <- act[,c('test_ligand','auroc','aupr','aupr_corrected','pearson','rank'),drop=FALSE]
       reproduction[[length(reproduction)+1]] <- data.frame(receiver,filter=filt,test,direction,status='GENESET_BACKGROUND_LIGAND_EXACT_MATCH_ORIGINAL_SCORES_REUSED')
     } else {
       act <- as.data.frame(nichenetr::predict_ligand_activities(geneset=gs,background_expressed_genes=bg,ligand_target_matrix=lt,potential_ligands=potential))
       stopifnot(all(c('test_ligand','aupr_corrected')%in%names(act)),any(is.finite(act$aupr_corrected)))
       act$rank <- rank(-act$aupr_corrected,ties.method='min')
       act <- act[,c('test_ligand','auroc','aupr','aupr_corrected','pearson','rank'),drop=FALSE]
     }
     act$receiver <- receiver;act$filter <- filt;act$gate <- gate;act$test <- test;act$direction <- direction;act$omit <- omit;act$geneset_n <- length(gs);act$background_n <- length(bg)
     activities[[length(activities)+1]] <- act
     cat('SENSITIVITY_DONE=',receiver,' ',filt,' ',gate,' ',test,' ',direction,' ',omit,' genes=',length(gs),'\n',sep='')
   }
 }
}
wr(do.call(rbind,activities),'ligand_activities.tsv');wr(do.call(rbind,statuses),'geneset_status.tsv');wr(do.call(rbind,coverage),'candidate_gate_coverage.tsv');wr(do.call(rbind,sets),'geneset_members.tsv');wr(do.call(rbind,reproduction),'baseline_reproduction.tsv')
writeLines(c('Step5: predefined 5% and10% target-expression gates, keeping >=2 libraries in either group, independent-library DE and sender/receptor eligibility unchanged. All and BG1 are separate.',
'Full uses all genes passing the frozen rule. Drop7 removes IL1R2/PENK/TTC39C/FOXP1/TNFRSF9/LAG3/MRE11 from BOTH the gene set and background; evaluates whether ligand ranking depends on the selected candidates. It is within-dataset sensitivity, not independent validation.',
'Only >=10 genes and>=500 background genes are scored. FDR and nominal P sets stay separately labeled. Failed FDR sets are not replaced by nominal sets. This is a prespecified follow-up prompted by observed low detection, not a new primary test.',
'Original 10% Full gene sets, backgrounds, ligand counts and statuses are checked against step3; their scores are reused unchanged. All other viable combinations are recomputed with the same NicheNet mouse prior.',
'Ranks are relative prioritization, not interaction P values. No new DE or causal sign is inferred. The original six-library design limits precision of sample correlations.',
'Output flags IL15/TNF plus all other ligands; no favorable threshold selection. No new database download.'),file.path(out,'methods.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
cat('ANALYSIS_LOGIC_PASS FROZEN_DE_AND_PRIOR_BASELINE_MATCHED\n')

R_ANALYZE
cat > "$OUT/finish.py" <<'PY_FINISH'
import os,json,tarfile,hashlib
from pathlib import Path
import pandas as pd
out=Path(os.environ['OUT']);root=Path(os.environ['ROOT'])
a=pd.read_csv(out/'tables/ligand_activities.tsv',sep='\t');s=pd.read_csv(out/'tables/geneset_status.tsv',sep='\t')
assert len(s)==64
cols=['receiver','filter','gate','test','direction','omit']
assert a[cols].drop_duplicates().shape[0]==(s.status=='PASS').sum()
assert not a.duplicated(cols+['test_ligand']).any()
assert (a['rank']>=1).all()
a[a.test_ligand.isin(['Il15','Tnf','Cd86','Tgfb1','Il1b','Tnfsf9','Il18','Cxcl16'])].sort_values(cols+['rank']).to_csv(out/'tables/focused_ligand_ranks.tsv',sep='\t',index=False)
key=['receiver','filter','gate','test','direction','test_ligand']
x=a[a.omit=='Full'][key+['rank','aupr_corrected']].merge(a[a.omit=='Drop7'][key+['rank','aupr_corrected']],on=key,how='outer',suffixes=('_Full','_Drop7'),validate='one_to_one')
x['top20_both']=(x.rank_Full<=20)&(x.rank_Drop7<=20)
x.to_csv(out/'tables/candidate_removal_sensitivity.tsv',sep='\t',index=False,na_rep='NA')
status=dict(status='PASS',stage='mechanism_step5_target_gate_sensitivity',output=str(out),previous=os.environ['PREV3'],sets_run=int((s.status=='PASS').sum()),sets_skipped=int((s.status!='PASS').sum()))
(out/'status.json').write_text(json.dumps(status,indent=2))
archive=root/'tmp/lung_mechanism_step5_return.tar.gz';part=archive.with_suffix('.tar.gz.part')
with tarfile.open(part,'w:gz') as t:
    for p in sorted(out.rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts:t.add(p,arcname=str(p.relative_to(out)))
with tarfile.open(part) as t:
    for m in t:
        if m.isfile():
            f=t.extractfile(m)
            while f.read(1024*1024):pass
part.replace(archive)
print('RETURN='+str(archive),flush=True)
print('RETURN_BYTES='+str(archive.stat().st_size),flush=True)
print('RETURN_SHA256='+hashlib.sha256(archive.read_bytes()).hexdigest(),flush=True)
print('LUNG_MECHANISM_STEP5_PASS',flush=True)

PY_FINISH
conda run --no-capture-output -n batch_bench python -m py_compile "$OUT/preflight.py" "$OUT/finish.py"
conda run --no-capture-output -n R_4.5.2 Rscript -e 'invisible(parse(file=commandArgs(TRUE)[1]));cat("R_PARSE_PASS\n")' "$OUT/analyze.R"
echo SYNTAX_CHECK=PASS
if [[ "${SYNTAX_ONLY:-0}" == 1 ]]; then exit 0; fi
conda run --no-capture-output -n batch_bench python -u "$OUT/preflight.py"
conda run --no-capture-output -n R_4.5.2 Rscript "$OUT/analyze.R" "$PREV1" "$PREV3" "$OUT" "$CACHE"
conda run --no-capture-output -n batch_bench python -u "$OUT/finish.py"
