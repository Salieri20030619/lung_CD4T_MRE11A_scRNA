suppressPackageStartupMessages({library(edgeR);library(AnnotationDbi);library(org.Mm.eg.db);library(fgsea)})
P<-'/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA';O<-file.path(P,'11_reports/death_priority_20261018_v3');V<-file.path(P,'11_reports/death_priority_20261018_v2')
wr<-function(x,n)write.table(as.data.frame(x),file.path(O,'tables',n),sep='\t',quote=FALSE,row.names=FALSE,na='NA',fileEncoding='UTF-8')
canon<-function(x){x<-toupper(x);a<-c(MRE11A='MRE11',ATP5A1='ATP5F1A',ATP5B='ATP5F1B',ATP5C1='ATP5F1C',ATP5G1='ATP5MC1',ATP5O='ATP5PO',PARK2='PRKN');i<-x%in%names(a);x[i]<-a[x[i]];x}
mm<-read.delim(file.path(V,'tables/module_members.tsv'));print(names(mm))
modules<-lapply(split(mm$gene[tolower(as.character(mm$retained))=='true'],mm$module[tolower(as.character(mm$retained))=='true']),unique)
qc<-list();ms<-list();me<-list();gs<-list();selected<-list()
for(acc in c('GSE121883','GSE324038')){
 d<-file.path(O,'external',acc);x<-as.matrix(read.delim(gzfile(file.path(d,'analysis_counts.tsv.gz')),row.names=1,check.names=FALSE))
 m<-read.delim(file.path(d,'analysis_samples.tsv'),check.names=FALSE);stopifnot(identical(colnames(x),m$column))
 db<-if(acc=='GSE121883')org.Mm.eg.db else {suppressPackageStartupMessages(library(org.Hs.eg.db));org.Hs.eg.db}
 ids<-sub('\\..*$','',rownames(x));valid<-intersect(ids,keys(db,keytype='ENSEMBL'));map<-AnnotationDbi::select(db,keys=valid,keytype='ENSEMBL',columns='SYMBOL');map<-unique(map[!is.na(map$SYMBOL),])
 n<-table(map$ENSEMBL);map<-map[map$ENSEMBL%in%names(n[n==1]),];symbol<-canon(map$SYMBOL[match(ids,map$ENSEMBL)]);symbol[is.na(symbol)]<-paste0('UNMAPPED_',ids[is.na(symbol)])
 wr(data.frame(ensembl=rownames(x),symbol=symbol),paste0(acc,'_mapping.tsv'));x<-rowsum(x,symbol,reorder=FALSE)
 for(ct in unique(m$context)){
  z<-m[m$context==ct,];xx<-x[,z$column,drop=FALSE];z$group<-factor(z$group,levels=c('WT','KO'))
  if(acc=='GSE324038'){z$donor<-factor(z$donor);design<-model.matrix(~donor+group,z)}else design<-model.matrix(~group,z)
  id<-paste(acc,ct,sep='_');y<-DGEList(xx);keep<-filterByExpr(y,design=design);y<-y[keep,,keep.lib.sizes=FALSE];y<-calcNormFactors(y);y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE);test<-glmQLFTest(fit,coef=ncol(design));tt<-topTags(test,n=Inf,sort.by='none')$table;tt$gene<-rownames(tt);tt$contrast<-id;tt$family<-paste0(id,'|all_tested_genes');wr(tt,paste0(id,'_DE.tsv'))
  lc<-cpm(y,log=TRUE,prior.count=2);wr(data.frame(gene=rownames(lc),lc,check.names=FALSE),paste0(id,'_logCPM.tsv'))
  zz<-z;zz$contrast<-id;zz$library_size<-y$samples$lib.size;zz$TMM<-y$samples$norm.factors;wr(zz,paste0(id,'_sample_QC.tsv'))
  pcs<-prcomp(t(lc),center=TRUE,scale.=FALSE)$x[,1:2];wr(data.frame(column=rownames(pcs),pcs,group=z$group,donor=z$donor),paste0(id,'_PCA.tsv'))
  for(mod in names(modules)){
   mem<-intersect(modules[[mod]],rownames(lc));if(length(mem)<3)next
   score<-colMeans(lc[mem,,drop=FALSE]);lmf<-lm.fit(design,score);resdf<-nrow(design)-qr(design)$rank;se<-sqrt(sum(lmf$residuals^2)/resdf*solve(crossprod(design))[ncol(design),ncol(design)]);effect<-tail(lmf$coefficients,1);pv<-2*pt(-abs(effect/se),resdf)
   me[[length(me)+1]]<-data.frame(contrast=id,program=mod,effect=effect,P=pv,n_libraries=nrow(z),n_biological=if(acc=='GSE324038')length(unique(z$donor)) else nrow(z),df=resdf,members=length(mem),family=paste0(id,'|module_effects'),member_genes=paste(mem,collapse=';'))
   ms[[length(ms)+1]]<-data.frame(contrast=id,column=z$column,donor=z$donor,group=z$group,program=mod,score=score)
  }
  ranks<-setNames(sign(tt$logFC)*sqrt(tt$F),tt$gene);ranks<-sort(ranks,decreasing=TRUE)
  sets<-lapply(modules,intersect,y=names(ranks));sets<-sets[lengths(sets)>=3]
  set.seed(20261018);fg<-as.data.frame(fgseaMultilevel(sets,ranks,minSize=3,maxSize=500,eps=1e-10,nproc=2));fg$leadingEdge<-vapply(fg$leadingEdge,paste,character(1),collapse=';');fg$contrast<-id;fg$family<-paste0(id,'|targeted_module_GSEA');gs[[length(gs)+1]]<-fg
  if(acc=='GSE121883'){
   k<-readRDS(file.path(V,'inputs/pathways_symbol.rds'))$KEGG
   set.seed(20261018);f<-as.data.frame(fgseaMultilevel(k,ranks,minSize=15,maxSize=500,eps=1e-10,nproc=2));f$leadingEdge<-vapply(f$leadingEdge,paste,character(1),collapse=';');f$contrast<-id;f$family<-paste0(id,'|KEGG_GSEA');wr(f,paste0(id,'_KEGG_GSEA.tsv'))
  }
  qc[[length(qc)+1]]<-data.frame(contrast=id,n_libraries=nrow(z),n_biological=if(acc=='GSE324038')length(unique(z$donor)) else nrow(z),tested=nrow(tt),nominal=sum(tt$PValue<.05),FDR05=sum(tt$FDR<.05),design_rank=qr(design)$rank,residual_df=nrow(z)-qr(design)$rank)
  selected[[length(selected)+1]]<-tt[tt$gene%in%c('HIF1A','IL2RA','ARL6IP1','NFKB1','NFKB2','PYCARD','GCLM','CASP3','CASP8','CASP9','CASP1','GSDMD','GSDME','BAX','BAK1','BCL2','BCL2L1','MCL1','CYCS','CYC1'),]
  cat('COMPLETE',id,'\n')
 }
}
eff<-do.call(rbind,me);eff$FDR<-ave(eff$P,eff$family,FUN=function(p)p.adjust(p,'BH'));wr(eff,'external_module_effects.tsv');wr(do.call(rbind,ms),'external_module_scores.tsv');wr(do.call(rbind,gs),'external_module_GSEA.tsv');wr(do.call(rbind,qc),'external_fit_QC.tsv');wr(do.call(rbind,selected),'external_key_genes.tsv')
writeLines(capture.output(sessionInfo()),file.path(O,'logs/external_R_session.txt'));writeLines('success',file.path(O,'logs/external_DE_complete.txt'))
