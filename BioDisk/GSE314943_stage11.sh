#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA
OUT=$(mktemp -d "$ROOT/02_analysis/GSE314943/stage11_$(date +%Y%m%d_%H%M%S)_XXXXXX")
trap 'rc=$?; if ((rc!=0)); then echo "GSE314943_STAGE11_FAIL EXIT_CODE=$rc"; fi' EXIT
cat > "$OUT/export.py" <<'PYTHON'
import sys
from pathlib import Path
import pandas as pd,numpy as np,scanpy as sc
from scipy.io import mmwrite
root=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/02_analysis/GSE314943')
out=Path(sys.argv[1]);(out/'input').mkdir();(out/'objects').mkdir();(out/'tables').mkdir()
a=sc.read_h5ad(root/'stage9_20260906_080538_wN3sjO/objects/K15_trajectory.h5ad')
pt=pd.read_csv(root/'stage9_20260906_080538_wN3sjO/tables/pseudotime_all_runs.tsv.gz',sep='\t',index_col=0)
assert a.n_obs==1734 and a.obs_names.isin(pt.index).all()
m=a.obs.copy();m.insert(0,'cell',a.obs_names)
m['PT']=pt.loc[a.obs_names,'K15_R1'].values;m['AltPT']=pt.loc[a.obs_names,'K15_R2'].values
m.to_csv(out/'input/cells.tsv',sep='\t',index=False)
mmwrite(out/'input/counts.mtx',a.layers['counts'].T)
pd.DataFrame({'gene_id':a.var_names,'symbol':a.var.symbol.values}).to_csv(out/'input/genes.tsv',sep='\t',index=False)
mac=sc.read_h5ad(root/'stage10_20260906_094533_ZUUfJn/objects/myeloid_candidates.h5ad')
mp={'1':'AMlike','2':'Spp1Mac','5':'C1qMac','11':'Folr2Mac'}
mac.obs['CoreMac']=mac.obs.M06.astype(str).isin(mp)
mac.obs.to_csv(out/'tables/myeloid_membership.tsv.gz',sep='\t')
mac=mac[mac.obs.CoreMac].copy();mac.obs['MacType']=mac.obs.M06.astype(str).map(mp).astype('category')
assert mac.n_obs==1796
pd.crosstab(mac.obs.MacType,mac.obs['sample']).to_csv(out/'tables/macrophage_sample.tsv',sep='\t')
mac.write_h5ad(out/'objects/macrophages_core.h5ad',compression='gzip')
print('EXPORT_PASS',flush=True)
PYTHON
cat > "$OUT/analysis.R" <<'RSCRIPT'
out<-commandArgs(TRUE)[1]
req<-c('tradeSeq','SingleCellExperiment','SummarizedExperiment','Matrix','BiocParallel','ggplot2')
miss<-req[!vapply(req,requireNamespace,logical(1),quietly=TRUE)]
if(length(miss)){
 options(repos=c(CRAN='https://cloud.r-project.org'),timeout=600)
 if(!requireNamespace('BiocManager',quietly=TRUE))install.packages('BiocManager')
 BiocManager::install(miss,ask=FALSE,update=FALSE)
}
stopifnot(all(vapply(req,requireNamespace,logical(1),quietly=TRUE)))
suppressPackageStartupMessages({library(tradeSeq);library(Matrix);library(SingleCellExperiment);library(ggplot2)})
dir.create(file.path(out,'figures'),showWarnings=FALSE)
wr<-function(t,n)write.table(t,file.path(out,'tables',n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
m<-read.delim(file.path(out,'input/cells.tsv'),stringsAsFactors=FALSE)
g<-read.delim(file.path(out,'input/genes.tsv'),stringsAsFactors=FALSE)
x<-as(readMM(file.path(out,'input/counts.mtx')),'CsparseMatrix');rownames(x)<-g$gene_id;colnames(x)<-m$cell
samples<-c(paste0('Saline_',1:3),paste0('Silica_',1:3))
m$sample<-factor(m$sample,levels=samples);m$group<-factor(m$group,levels=c('Saline','Silica'))
stopifnot(ncol(x)==1734,all(is.finite(m$PT)),all(colSums(x)>0))
off<-log(Matrix::colSums(x))
det<-sapply(samples,function(s)Matrix::rowSums(x[,m$sample==s,drop=FALSE]>0))
keep<-Matrix::rowSums(x)>=(50) & Matrix::rowSums(x>0)>=30 & rowSums(det>=5)>=3 & rowSums(det[,1:3,drop=FALSE]>=5)>=1 & rowSums(det[,4:6,drop=FALSE]>=5)>=1
wr(data.frame(g,total_umi=Matrix::rowSums(x),detected_cells=Matrix::rowSums(x>0),tested=keep),'gene_filter.tsv')
x<-x[keep,,drop=FALSE];stopifnot(nrow(x)>100)
U<-model.matrix(~sample,data=m);stopifnot(qr(U)$rank==ncol(U))
pt<-matrix(m$PT,ncol=1,dimnames=list(m$cell,'L1'));w<-matrix(1,nrow(m),1,dimnames=dimnames(pt))
bp<-BiocParallel::MulticoreParam(workers=4,RNGseed=42,progressbar=FALSE)
# Prespecified small knot range; two independent random gene subsets.
ks<-3:6;ics<-list()
for(seed in c(42,43)){
 set.seed(seed);message('EVALUATE_K=',seed)
 z<-evaluateK(counts=x,pseudotime=pt,cellWeights=w,U=U,offset=off,k=ks,nGenes=min(100,nrow(x)),plot=FALSE,verbose=FALSE,parallel=TRUE,BPPARAM=bp)
 stopifnot(is.matrix(z),ncol(z)==length(ks));ics[[as.character(seed)]]<-z
 wr(data.frame(gene_id=rownames(z),z,check.names=FALSE),paste0('evaluateK_',seed,'.tsv'))
}
ic<-do.call(rbind,ics);ok<-apply(ic,1,function(v)all(is.finite(v)));stopifnot(sum(ok)>=50)
delta<-sweep(ic[ok,,drop=FALSE],1,apply(ic[ok,,drop=FALSE],1,min),'-')
meanDelta<-colMeans(delta);nk<-min(ks[meanDelta<=min(meanDelta)+2])
wr(data.frame(knots=ks,mean_delta_AIC=meanDelta,selected=ks==nk),'knot_choice.tsv')
message('SELECTED_K=',nk)
fit<-function(xx,mm,pp,oo,conditions=FALSE){
 uu<-model.matrix(~sample,data=transform(mm,sample=droplevels(sample)))
 if(conditions){
  # Within-group sample contrasts: sum to zero across biological samples,
  # avoiding redundant sample/group indicators alongside condition smoothers.
  blocks<-lapply(levels(mm$group),function(gr){
   sm<-levels(droplevels(mm$sample[mm$group==gr]));h<-contr.helmert(length(sm))
   z<-matrix(0,nrow(mm),ncol(h));i<-which(mm$group==gr);z[i,]<-h[match(as.character(mm$sample[i]),sm),,drop=FALSE];z
  })
  uu<-cbind(Intercept=1,do.call(cbind,blocks))
 }
 stopifnot(qr(uu)$rank==ncol(uu))
 set.seed(42)
 do.call(tradeSeq::fitGAM,c(list(counts=xx,pseudotime=matrix(pp,ncol=1),cellWeights=matrix(1,nrow(mm),1),U=uu,offset=oo,nknots=nk,verbose=FALSE,parallel=TRUE,BPPARAM=bp),if(conditions)list(conditions=droplevels(mm$group)) else list()))
}
report<-function(sce,fun,name){
 z<-as.data.frame(fun(sce));z$gene_id<-rownames(z);z$symbol<-g$symbol[match(z$gene_id,g$gene_id)]
 cv<-rowData(sce)$tradeSeq$converged
 stopifnot(length(cv)==nrow(sce));z$converged<-as.logical(cv[match(z$gene_id,rownames(sce))])
 z$pvalue[is.na(z$converged)|!z$converged]<-NA_real_;z$FDR<-p.adjust(z$pvalue,'BH')
 z<-z[order(z$FDR,z$pvalue,na.last=TRUE),];wr(z,paste0(name,'.tsv'));z
}
message('FIT_COMMON_START')
sce<-fit(x,m,m$PT,off)
saveRDS(sce,file.path(out,'objects/common_GAM.rds'))
a<-report(sce,function(z)associationTest(z),'associationTest')
b<-report(sce,function(z)startVsEndTest(z),'startVsEndTest')
# Condition tests only use the shared 5-95% sample support range.
ranges<-t(vapply(samples,function(s)quantile(m$PT[m$sample==s],c(.05,.95)),numeric(2)))
lo<-max(ranges[,1]);hi<-min(ranges[,2]);wr(data.frame(sample=samples,lower=ranges[,1],upper=ranges[,2]),'condition_support.tsv')
use<-m$PT>=lo & m$PT<=hi
conditionOK<-hi>lo & all(table(factor(m$sample[use],levels=samples))>=30)
ct<-NULL;cm<-NULL;status<-data.frame(analysis='conditionTest',status=if(conditionOK)'READY' else 'INSUFFICIENT_SHARED_SUPPORT')
if(conditionOK){
 mm<-m[use,];xx<-x[,use,drop=FALSE];oo<-off[use]
 ck<-Matrix::rowSums(xx>0)>=30 & Matrix::rowSums(xx)>=50
 xx<-xx[ck,,drop=FALSE];stopifnot(nrow(xx)>100)
 message('FIT_CONDITION_START')
 cm<-fit(xx,mm,mm$PT,oo,TRUE);saveRDS(cm,file.path(out,'objects/condition_GAM.rds'))
 ct<-report(cm,function(z)conditionTest(z),'conditionTest');status$status<-'COMPLETED'
}
wr(status,'condition_status.tsv')
# Candidate validation: predeclared target plus top ten converged significant
# genes per test, without re-labelling these selected-gene P values as a new FDR screen.
gettop<-function(z){if(is.null(z))return(character());head(z$gene_id[is.finite(z$FDR)&z$FDR<.05],10)}
candidates<-unique(c(intersect('ENSMUSG00000031928',rownames(sce)),gettop(a),gettop(b),gettop(ct)))
wr(data.frame(gene_id=candidates,symbol=g$symbol[match(candidates,g$gene_id)]),'candidate_genes.tsv')
if(length(candidates)){
 pred<-as.data.frame(predictSmooth(sce,gene=candidates,nPoints=100,tidy=TRUE));wr(pred,'common_fitted_curves.tsv')
 pred$label<-sprintf('G%02d',match(pred$gene,candidates));pred$scaled<-ave(pred$yhat,pred$gene,FUN=function(v)if(sd(v)>0)as.numeric(scale(v))else rep(0,length(v)))
 labels<-data.frame(label=sprintf('G%02d',seq_along(candidates)),gene_id=candidates,symbol=g$symbol[match(candidates,g$gene_id)]);wr(labels,'curve_labels.tsv')
 p<-ggplot(pred,aes(time,scaled))+geom_line(color='#0072B2')+facet_wrap(~label,ncol=4)+labs(x='DPT',y='Scaled fit')+theme_classic(base_size=11)
 for(ext in c('png','pdf'))ggsave(file.path(out,'figures',paste0('common_curves.',ext)),p,width=10,height=max(4,ceiling(length(candidates)/4)*2),dpi=350,bg='white')
 writeLines(c('Sample-adjusted common NB-GAM curves, standardized within each gene for display. Genes are Mre11a if eligible and top ten FDR<0.05 genes per primary test. Curves use tradeSeq prediction defaults; y is not an absolute expression measurement.',paste(labels$label,labels$symbol,sep=' = ')),file.path(out,'figures/common_curves_caption.txt'))
 # Alternate root and leave-one-sample-out models: same selected candidates and knots.
 checks<-list();configs<-c('AltRoot',paste0('Drop_',samples))
 for(cfg in configs){
  take<-if(cfg=='AltRoot')rep(TRUE,nrow(m))else as.character(m$sample)!=sub('^Drop_','',cfg)
  mx<-m[take,];mx$sample<-droplevels(mx$sample)
  pp<-if(cfg=='AltRoot')mx$AltPT else mx$PT
  message('CANDIDATE_CHECK=',cfg)
  cs<-fit(x[candidates,take,drop=FALSE],mx,pp,off[take]);saveRDS(cs,file.path(out,'objects',paste0(cfg,'.rds')))
  za<-report(cs,function(z)associationTest(z),paste0(cfg,'_association'))
  ze<-report(cs,function(z)startVsEndTest(z),paste0(cfg,'_start_end'))
  pr<-as.data.frame(predictSmooth(cs,gene=candidates,nPoints=50,tidy=TRUE))
  wr(pr,paste0(cfg,'_curves.tsv'))
  for(id in candidates){v<-pr[pr$gene==id,];v<-v[order(v$time),];checks[[length(checks)+1]]<-data.frame(run=cfg,gene_id=id,association_P=za$pvalue[match(id,za$gene_id)],start_end_P=ze$pvalue[match(id,ze$gene_id)],end_minus_start=tail(v$yhat,1)-v$yhat[1])}
 }
 wr(do.call(rbind,checks),'candidate_sensitivity.tsv')
 if(!is.null(cm)){
  cc<-intersect(candidates,rownames(cm))
  if(length(cc)){
   cp<-as.data.frame(predictSmooth(cm,gene=cc,nPoints=100,tidy=TRUE));wr(cp,'condition_fitted_curves.tsv')
   cp$label<-labels$label[match(cp$gene,labels$gene_id)]
   cp$scaled<-ave(cp$yhat,cp$gene,FUN=function(v)if(sd(v)>0)as.numeric(scale(v))else rep(0,length(v)))
   pp<-ggplot(cp,aes(time,scaled,color=condition))+geom_line()+facet_wrap(~label,ncol=4)+scale_color_manual(values=c(Saline='#0072B2',Silica='#D55E00'))+labs(x='DPT',y='Scaled fit')+theme_classic(base_size=11)
   for(ext in c('png','pdf'))ggsave(file.path(out,'figures',paste0('condition_curves.',ext)),pp,width=10,height=max(4,ceiling(length(cc)/4)*2),dpi=350,bg='white')
   writeLines(c('Condition-specific curves in the shared sample-supported DPT interval; both curves standardized jointly within each gene. No separate condition-wise rescaling. Within-group sample fixed contrasts included; cell-level model evidence, not a donor-level test.',paste(labels$label,labels$symbol,sep=' = ')),file.path(out,'figures/condition_curves_caption.txt'))
  }
 }
}
sumtest<-function(z,n){if(is.null(z))return(data.frame(test=n,tested=0,FDR05=0));data.frame(test=n,tested=sum(is.finite(z$pvalue)),FDR05=sum(z$FDR<.05,na.rm=TRUE))}
wr(rbind(sumtest(a,'Association'),sumtest(b,'StartEnd'),sumtest(ct,'Condition')),'tradeSeq_summary.tsv')
# Tests have different columns, so export each target test independently.
for(n in c('Association','StartEnd','Condition')){
 z<-switch(n,Association=a,StartEnd=b,Condition=ct)
 if(!is.null(z))wr(z[z$gene_id=='ENSMUSG00000031928',,drop=FALSE],paste0('Mre11a_',n,'.tsv'))
}
writeLines(c('Reference K15_R1, alternate root K15_R2; 1734 Naive/ActCD4 cells. Raw counts; offset=log full per-cell library UMI before gene filtering.',
 'Expression filter: >=50 UMIs, >=30 detected cells, >=5 detected cells in >=3 samples including >=1 sample per group. Mre11a is subject to the same filter.',
 'evaluateK: k3:6 on two seeded random subsets of100 eligible genes, choose smallest K within2 of minimum mean per-gene delta AIC.',
 'Common model U=intercept plus sample fixed effects. Conditions model uses shared sample support (intersection of per-sample 5-95% DPT), >=30 cells/sample, intercept plus within-group sample contrasts without redundant group indicators. Condition-specific curves compare matched absolute DPT, not separate group rescaling.',
 'BH separately per primary test; nonconverged fits have NA P/FDR. Primary tests are cell-level model evidence; sample fixed effects do not create donor-level independent errors.',
 'Candidate checks: same K, alternate root and leave-one-sample-out common models; selected-gene P/FDR are descriptive sensitivity checks, not independent confirmation. ConditionTest has not been donor-bootstrap validated.',
 'Macrophage core retained separately: M1 AMlike, M2 Spp1Mac, M5 C1qMac, M11 Folr2Mac. Subtype labels express marker-supported states, not developmental origin. M2 has control counts9,7,0; no complete 3+3 subtype communication inference.',
 'State ordering at one endpoint is not observed time or a causal test of Mre11a.'),file.path(out,'tables/methods.txt'))
capture.output(sessionInfo(),file=file.path(out,'tables/sessionInfo.txt'))
message('GSE314943_STAGE11_ANALYSIS_PASS')
RSCRIPT
conda run -n batch_bench python -m py_compile "$OUT/export.py"
conda run -n R_4.5.2 Rscript -e 'invisible(parse(file=commandArgs(TRUE)[1]));cat("R_PARSE_PASS\n")' "$OUT/analysis.R" > "$OUT/syntax.log" 2>&1
echo 'SYNTAX_CHECK=PASS'
echo "OUTPUT=$OUT"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
conda run --no-capture-output -n batch_bench python -u "$OUT/export.py" "$OUT"
conda run --no-capture-output -n R_4.5.2 Rscript "$OUT/analysis.R" "$OUT"
tar -C "$OUT" -czf "$ROOT/tmp/GSE314943_stage11_return.tar.gz.part" tables figures analysis.R export.py syntax.log
mv "$ROOT/tmp/GSE314943_stage11_return.tar.gz.part" "$ROOT/tmp/GSE314943_stage11_return.tar.gz"
echo 'GSE314943_STAGE11_PASS'
