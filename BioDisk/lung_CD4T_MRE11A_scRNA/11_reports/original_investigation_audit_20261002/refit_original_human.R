suppressPackageStartupMessages({library(Matrix);library(edgeR)})
P<-'/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA';O<-file.path(P,'11_reports/original_investigation_audit_20261002/results');I<-file.path(P,'02_analysis/mechanism_step2_20260927_154919_EsN4gk/input')
rd<-function(p)read.delim(p,check.names=FALSE,stringsAsFactors=FALSE);wr<-function(d,n)write.table(d,file.path(O,n),sep='\t',row.names=FALSE,quote=FALSE)
x<-as(readMM(file.path(I,'human_counts.mtx')),'CsparseMatrix');m<-rd(file.path(I,'human_cells.tsv'));g<-rd(file.path(I,'human_genes.tsv'));rownames(x)<-g$symbol;colnames(x)<-m$cell
old<-rd(file.path(P,'10_tables/CD4_final/pathway_v1/source_adjusted_full_DEG.tsv'));audit<-list();focal<-list();rosters<-list();bridges<-list();samplevalues<-list()
for(scope in c('Primary','Core')){
 take<-if(scope=='Primary')rep(TRUE,nrow(m))else m$subtype!='CD4state';z<-m[take,];ok<-names(which(table(z$sample)>=20));z<-z[z$sample%in%ok,];sm<-unique(z[,c('sample','group','series')]);sm<-sm[order(sm$sample),];names(sm)<-c('Subject','Group','Series')
 pb<-do.call(cbind,lapply(sm$Subject,function(ss)Matrix::rowSums(x[,z$cell[z$sample==ss],drop=FALSE])));rownames(pb)<-g$symbol;colnames(pb)<-sm$Subject
 sm$Group<-factor(sm$Group,levels=c('Control','CWP'));sm$Series<-factor(sm$Series,levels=c('G227','G288','G135'));design<-model.matrix(~Series+Group,data=sm);rownames(design)<-sm$Subject
 y<-DGEList(counts=pb,samples=sm);keep<-filterByExpr(y,design=design);y<-y[keep,,keep.lib.sizes=FALSE];y<-calcNormFactors(y);y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE);tt<-topTags(glmQLFTest(fit,coef='GroupCWP'),n=Inf,sort.by='none')$table;tt$gene<-rownames(tt)
 j<-merge(tt,old[old$analysis==scope,],by='gene',suffixes=c('_refit','_original'));stopifnot(nrow(j)==nrow(tt),nrow(tt)==sum(old$analysis==scope))
 metrics<-c('logFC','logCPM','F','PValue','FDR');delta<-vapply(metrics,function(v)max(abs(j[[paste0(v,'_refit')]]-j[[paste0(v,'_original')]])),numeric(1));stopifnot(all(delta<1e-8))
 audit[[scope]]<-data.frame(scope=scope,n=nrow(sm),cells=nrow(z),genes=nrow(tt),status='EXACT_ORIGINAL_SOURCE_ADJUSTED_REFIT',t(delta));q<-j[j$gene%in%c('IKZF2','STAT1'),];q$scope<-scope;focal[[scope]]<-q
 sm$cells<-as.numeric(table(factor(z$sample,levels=sm$Subject)));sm$scope<-scope;rosters[[scope]]<-sm
 if(scope=='Core'){
  sc<-rd(file.path(P,'02_analysis/scenic1_20260929_103513_xRoEnr/Human/All/cells.tsv'));stopifnot(setequal(z$cell,sc$cell),setequal(sm$Subject,sc$sample));cat('HUMAN_CORE_EXACT_SCENIC_CELL_SET',nrow(z),nrow(sm),'\n')
  lcpm<-cpm(y,log=TRUE,prior.count=2)
  for(gene in c('IKZF2','STAT1')){
   rawmean<-vapply(sm$Subject,function(ss){ii<-which(m$cell%in%z$cell[z$sample==ss]);mean(as.numeric(x[gene,ii])/m$library[ii]*1e4)},numeric(1))
   rawlib<-vapply(sm$Subject,function(ss)sum(z$library[z$sample==ss]),numeric(1));weighted<-as.numeric(pb[gene,])/rawlib*1e4
   samplevalues[[gene]]<-data.frame(gene=gene,sm,raw_gene_count=as.numeric(pb[gene,]),raw_library=rawlib,mean_cell_CP10k=rawmean,weighted_pseudobulk_CP10k=weighted,TMM_logCPM=as.numeric(lcpm[gene,]))
   values<-list(log1p_mean_cell_CP10k=log1p(rawmean),log1p_weighted_pseudobulk_CP10k=log1p(weighted),same_edgeR_TMM_logCPM=as.numeric(lcpm[gene,]))
   for(method in names(values)){
    tmp<-sm;tmp$Y<-values[[method]];lf<-lm(Y~Series+Group,data=tmp);co<-summary(lf)$coefficients['GroupCWP',]
    bridges[[paste(gene,method)]]<-data.frame(gene=gene,scope=scope,method=method,model='OLS_Source',n=nrow(sm),coefficient=co[1],se=co[2],PValue=co[4],df=df.residual(lf))
   }
  }
 }
 cat('HUMAN_REFIT_PASS',scope,'Pdiff',delta['PValue'],'\n')
}
wr(do.call(rbind,audit),'human_original_edgeR_exact_refit_audit.tsv');wr(do.call(rbind,focal),'human_original_TF_refit_comparison.tsv');wr(do.call(rbind,rosters),'human_original_edgeR_roster.tsv')
wr(do.call(rbind,bridges),'human_same_roster_estimand_bridge.tsv');wr(do.call(rbind,samplevalues),'human_same_roster_RNA_sample_values.tsv')
