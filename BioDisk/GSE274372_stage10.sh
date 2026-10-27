#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA
RUN=$(mktemp -d "$ROOT/02_analysis/GSE274372/stage10_$(date +%Y%m%d_%H%M%S)_XXXXXX")
export ROOT RUN OPENBLAS_NUM_THREADS=1
trap 'rc=$?; echo "GSE274372_STAGE10_FAIL EXIT_CODE=$rc"; exit "$rc"' ERR
cat > "$RUN/export.py" <<'PY'
import os
from pathlib import Path
import scanpy as sc,numpy as np,pandas as pd
root=Path(os.environ['ROOT']);run=Path(os.environ['RUN']);out=run/'tables';out.mkdir()
a=sc.read_h5ad(root/'02_analysis/GSE274372/stage9_20260908_140721_Yj56oq/objects/CD4_author_reclustered.h5ad')
assert a.n_obs==888
x=a.layers['counts'];assert np.all(x.data>=0) and np.allclose(x.data,np.round(x.data))
sy=a.var['symbol'].astype('string').fillna('');ix=np.flatnonzero(sy.eq('Mre11a'));assert len(ix)==1
cl=a.obs.R08.astype(str)
scopes={'TLike':cl.isin(['0','1']),'NoMix':~cl.isin(['4','5','6']),'Author':np.ones(a.n_obs,dtype=bool)}
rows=[]
for name,selection in scopes.items():
 cols={}
 for sample in sorted(a.obs['sample'].astype(str).unique()):
  mask=np.asarray(selection)&a.obs['sample'].astype(str).eq(sample).to_numpy()
  v=np.asarray(x[mask].sum(axis=0)).ravel().astype('int64');cols[sample]=v
  target=np.asarray(x[mask][:,ix].sum(axis=1)).ravel();n=int(mask.sum())
  rows.append([name,sample,int(a.obs.loc[a.obs['sample'].astype(str).eq(sample),'day'].iloc[0]),n,int((target>0).sum()),int(target.sum()),int(v.sum())])
 pd.DataFrame(cols,index=a.var_names).to_csv(out/f'{name}_counts.tsv.gz',sep='\t')
pd.DataFrame(rows,columns=['scope','sample','day','cells','Mre_positive','Mre_UMI','library_UMI']).to_csv(out/'sample_summary.tsv',sep='\t',index=False)
a.var[['symbol']].to_csv(out/'features.tsv.gz',sep='\t')
a.obs.to_csv(out/'cell_metadata.tsv.gz',sep='\t')
(out/'scope_notes.txt').write_text('TLike: R08 clusters0/1, activated-regulatory and naive T-like author CD4 candidates, not definitively pure CD4. NoMix: excludes clusters4/5/6; retains uncertain clusters2/3. Author: all888. No identities overwritten. Independent samples are analysis units. Minimum10 cells/sample is an operational inclusion rule, not a universal validity threshold. Per-day biological replicates remain small. No Mre11a-positive gating.\n')
print('EXPORT_PASS',flush=True)
PY
cat > "$RUN/de.R" <<'RS'
suppressPackageStartupMessages(library(edgeR))
out <- file.path(Sys.getenv('RUN'),'tables')
meta <- read.delim(file.path(out,'sample_summary.tsv'),check.names=FALSE)
features <- read.delim(gzfile(file.path(out,'features.tsv.gz')),row.names=1,check.names=FALSE)
status <- list();target <- list()
for (scope in c('TLike','NoMix','Author')) {
 m <- meta[meta$scope==scope,];m$included <- m$cells>=10 & m$library_UMI>0
 write.table(m,file.path(out,paste0(scope,'sample_inclusion.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
 m <- m[m$included,];m$day <- factor(m$day,levels=c(0,7,28,56))
 n <- table(m$day)
 if(any(n<2)) {status[[scope]]<-data.frame(scope=scope,status='SKIP: fewer than2 samples in a day');next}
 counts <- as.matrix(read.delim(gzfile(file.path(out,paste0(scope,'_counts.tsv.gz'))),row.names=1,check.names=FALSE))[,m$sample,drop=FALSE]
 design <- model.matrix(~0+day,data=m);colnames(design)<-paste0('D',levels(m$day))
 stopifnot(qr(design)$rank==ncol(design),nrow(design)>ncol(design))
 y <- DGEList(counts=counts)
 keep <- filterByExpr(y,design=design)
 audit <- data.frame(gene=rownames(y),symbol=features[rownames(y),'symbol'],tested=keep)
 write.table(audit,file.path(out,paste0(scope,'gene_filter.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
 target[[scope]] <- cbind(scope=scope,audit[audit$symbol %in% 'Mre11a',])
 y <- calcNormFactors(y[keep,,keep.lib.sizes=FALSE]);y<-estimateDisp(y,design,robust=TRUE)
 fit <- glmQLFit(y,design,robust=TRUE)
 # Three day-versus-baseline contrasts jointly test any time effect.
 cm <- cbind(c(-1,1,0,0),c(-1,0,1,0),c(-1,0,0,1))
 tests <- c(list(Omnibus=glmQLFTest(fit,contrast=cm)),setNames(lapply(1:3,function(i)glmQLFTest(fit,contrast=cm[,i])),c('D7_D0','D28_D0','D56_D0')))
 for (nm in names(tests)) {
  z<-topTags(tests[[nm]],n=Inf,sort.by='none')$table
  z<-cbind(gene=rownames(z),symbol=features[rownames(z),'symbol'],z)
  write.table(z,file.path(out,paste0(scope,'_',nm,'.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
 }
 saveRDS(list(y=y,fit=fit,design=design,meta=m),file.path(Sys.getenv('RUN'),paste0(scope,'_edgeR.rds')))
 status[[scope]]<-data.frame(scope=scope,status=paste('PASS samples',nrow(m),'genes',sum(keep)))
}
write.table(do.call(rbind,status),file.path(out,'DE_status.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
if(length(target))write.table(do.call(rbind,target),file.path(out,'Mre_test_eligibility.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'R_session.txt'))
cat('GSE274372_STAGE10_ANALYSIS_PASS\n')
RS
bash -n "$0"
conda run -n batch_bench python -m py_compile "$RUN/export.py"
conda run -n R_4.5.2 Rscript -e 'invisible(parse(file=file.path(Sys.getenv("RUN"),"de.R")))'
echo SYNTAX_CHECK=PASS
echo "OUTPUT=$RUN"
conda run --no-capture-output -n batch_bench python -u "$RUN/export.py"
conda run --no-capture-output -n R_4.5.2 Rscript "$RUN/de.R"
tar -czf "$ROOT/tmp/GSE274372_stage10_return.tar.gz" -C "$RUN" tables export.py de.R
echo GSE274372_STAGE10_PASS
