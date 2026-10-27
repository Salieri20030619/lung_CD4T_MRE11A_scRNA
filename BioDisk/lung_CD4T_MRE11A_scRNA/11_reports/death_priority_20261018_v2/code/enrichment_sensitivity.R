suppressPackageStartupMessages({library(edgeR);library(fgsea)})
P<-'/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA';O<-file.path(P,'11_reports/death_priority_20261018_v2');B<-file.path(P,'11_reports/three_cohort_four_axis_20261017_v16')
wr<-function(x,n)write.table(as.data.frame(x),file.path(O,'tables',n),sep='\t',quote=FALSE,row.names=FALSE,na='NA',fileEncoding='UTF-8')
pws<-readRDS(file.path(O,'inputs/pathways_symbol.rds'));nm<-read.delim(file.path(O,'tables/pathway_names.tsv'));desc<-setNames(nm$description,nm$pathway)
# Refine only imprecise or bounded first-pass probabilities; keep full families.
f<-read.delim(file.path(O,'tables/joint_GSEA_all.tsv'));file.copy(file.path(O,'tables/joint_GSEA_all.tsv'),file.path(O,'inputs/joint_GSEA_firstpass.tsv'),overwrite=FALSE)
for(sc in unique(f$scope))for(db in unique(f$database)){
 ii<-which(f$scope==sc&f$database==db);bad<-ii[is.na(f$log2err[ii])|f$pval[ii]<=1e-10]
 if(!length(bad))next
 d<-read.delim(file.path(O,'tables',paste0('joint_rank_',sc,'.tsv')));stat<-setNames(d$rank,d$gene)
 set.seed(20261018);z<-as.data.frame(fgseaMultilevel(pws[[db]][f$pathway[bad]],stat,minSize=15,maxSize=500,eps=0,nPermSimple=10000,nproc=2))
 for(j in seq_len(nrow(z))){
  k<-bad[match(z$pathway[j],f$pathway[bad])]
  for(col in c('pval','log2err','ES','NES','size'))f[k,col]<-z[j,col]
  f$leadingEdge[k]<-paste(z$leadingEdge[[j]],collapse=';')
 }
 f$padj[ii]<-p.adjust(f$pval[ii],'BH')
}
f$numeric_status<-ifelse(is.na(f$log2err),'precision_warning','OK');wr(f,'joint_GSEA_all.tsv')
fit<-readRDS(file.path(B,'inputs/edgeR_fits.rds'))[['joint_stratified__Core']]
m<-fit$meta;s<-read.delim(file.path(O,'tables/sample_scores.tsv'));s<-s[s$variant=='Raw'&s$scope%in%c('Naive','ActCD4','Treg')&s$dataset!='GSE274372',]
for(sc in c('ActCD4','Treg')){ss<-s[s$scope==sc,];tot<-aggregate(cells~dataset+sample,s,sum);ss<-merge(ss,tot,by=c('dataset','sample'),suffixes=c('','_total'));m[[paste0('fraction_',sc)]]<-(ss$cells/ss$cells_total)[match(paste(m$dataset,m$sample),paste(ss$dataset,ss$sample))]}
m$condition<-factor(m$group,levels=c('Saline','Silica'));m$stratum<-factor(paste(m$dataset,m$day))
design<-model.matrix(~stratum+fraction_ActCD4+fraction_Treg+condition,m)
stopifnot(qr(design)$rank==ncol(design));y<-estimateDisp(fit$y,design,robust=TRUE);ff<-glmQLFit(y,design,robust=TRUE);q<-glmQLFTest(ff,coef=which(colnames(design)=='conditionSilica'));d<-topTags(q,n=Inf,sort.by='none')$table;d$gene<-rownames(d);d$family<-'Core_composition_adjusted_disease';wr(d,'Core_composition_adjusted_DE.tsv')
wr(data.frame(m,design,check.names=FALSE),'Core_composition_adjusted_design.tsv');wr(data.frame(n=nrow(design),rank=qr(design)$rank,residual_df=nrow(design)-ncol(design),condition_number=kappa(design),genes=nrow(d)),'Core_composition_adjusted_diagnostics.tsv')
d$rank<-sign(d$logFC)*sqrt(pmax(d$F,0));d<-d[order(-d$rank,d$gene),];wr(d,'Core_composition_adjusted_rank.tsv')
rr<-list()
for(db in names(pws)){
 set.seed(20261018);z<-as.data.frame(fgseaMultilevel(pws[[db]],setNames(d$rank,d$gene),minSize=15,maxSize=500,eps=0,nPermSimple=10000,nproc=2));z$scope<-'Core';z$database<-db;z$family<-paste('Core_composition_adjusted',db,sep='|');z$description<-unname(desc[z$pathway]);z$leadingEdge<-vapply(z$leadingEdge,paste,character(1),collapse=';');rr[[db]]<-z
}
wr(do.call(rbind,rr),'Core_composition_adjusted_GSEA.tsv');writeLines('success',file.path(O,'logs/composition_sensitivity_complete.txt'))
cat('COMPOSITION_SENSITIVITY_COMPLETE\n')
