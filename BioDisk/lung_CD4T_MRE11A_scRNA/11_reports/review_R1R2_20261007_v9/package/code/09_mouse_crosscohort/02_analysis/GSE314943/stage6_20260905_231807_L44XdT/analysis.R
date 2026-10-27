args <- commandArgs(trailingOnly=TRUE)
src <- args[1]; out <- args[2]
required <- c('edgeR','statmod','ggplot2','ggrepel')
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) stop('Missing R packages: ',paste(missing,collapse=', '))
suppressPackageStartupMessages({library(edgeR);library(ggplot2);library(ggrepel)})
set.seed(42)
for(d in c('tables','figures','models')) dir.create(file.path(out,d),showWarnings=FALSE)
readtab <- function(f) read.delim(file.path(src,'tables',f),check.names=FALSE,stringsAsFactors=FALSE)
write <- function(x,f) write.table(x,file.path(out,'tables',f),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
meta <- readtab('PB_sample_metadata.tsv')
genes <- readtab('genes.tsv')
samples <- c(paste0('Saline_',1:3),paste0('Silica_',1:3))
group <- factor(c(rep('Saline',3),rep('Silica',3)),levels=c('Saline','Silica'))
design <- model.matrix(~group)
rownames(design) <- samples
stopifnot(qr(design)$rank==2L, !anyDuplicated(genes$gene_id))
write(data.frame(sample=samples,group=group,GSM=c('GSM9473144','GSM9529135','GSM9529136','GSM9529140','GSM9529141','GSM9529142')),'design.tsv')
# Predeclared background sensitivity; no count decontamination or Harmony counts.
background <- c('Sftpc','Sftpa1','Sftpb','Sftpd','Scgb1a1','Scgb3a2','S100a8','S100a9','Lyz2','Igkc','Ighm','Iglc1','Iglc2','Iglc3')
excluded <- genes$symbol %in% background | grepl('^Hb[ab]',genes$symbol)
write(genes[excluded,,drop=FALSE],'NoBG_excluded.tsv')
summary <- list(); target <- list(); resultlist <- list()
plot_save <- function(p,stem,caption,w=7,h=5){
  ggsave(file.path(out,'figures',paste0(stem,'.png')),p,width=w,height=h,dpi=350,bg='white')
  ggsave(file.path(out,'figures',paste0(stem,'.pdf')),p,width=w,height=h,bg='white')
  writeLines(caption,file.path(out,'figures',paste0(stem,'_caption.txt')))
}
for(setname in c('Core','Naive','ActCD4','Treg')){
  message('SET_START=',setname)
  m <- meta[meta$set==setname,,drop=FALSE];m <- m[match(samples,m$sample),,drop=FALSE]
  stopifnot(identical(m$sample,samples),identical(m$group,as.character(group)))
  if(setname!='Treg') stopifnot(all(m$cells>=20))
  stopifnot(all(m$cells>0))
  status <- if(setname=='Treg') 'exploratory_low_cells' else 'primary'
  con <- gzfile(file.path(src,'tables',paste0('PB_',setname,'_counts.tsv.gz')),'rt')
  dat <- read.delim(con,check.names=FALSE,stringsAsFactors=FALSE);close(con)
  ids <- dat[[1]]; counts <- as.matrix(dat[,-1,drop=FALSE]);rownames(counts)<-ids
  stopifnot(identical(colnames(counts),samples),identical(ids,genes$gene_id),
            all(is.finite(counts)),all(counts>=0),all(counts==round(counts)),
            all(colSums(counts)==m$library_umi))
  for(arm in c('Full','NoBG')){
    tag <- paste(setname,arm,sep='_')
    inc <- if(arm=='Full') rep(TRUE,nrow(genes)) else !excluded
    y <- DGEList(counts=counts[inc,,drop=FALSE],group=group,genes=genes[inc,,drop=FALSE])
    keep <- filterByExpr(y,design=design)
    write(data.frame(gene_id=y$genes$gene_id,symbol=y$genes$symbol,tested=keep),paste0(tag,'_filter.tsv'))
    stopifnot(sum(keep)>100)
    y <- y[keep,,keep.lib.sizes=FALSE]
    y <- calcNormFactors(y,method='TMM')
    y <- estimateDisp(y,design,robust=TRUE)
    fit <- glmQLFit(y,design,robust=TRUE)
    test <- glmQLFTest(fit,coef=2)
    tt <- topTags(test,n=Inf,sort.by='PValue')$table
    tt$background <- tt$gene_id %in% genes$gene_id[excluded]
    tt$direction <- ifelse(tt$logFC>0,'Up','Down')
    tt$set <- setname;tt$arm <- arm;tt$status <- status
    write(tt,paste0(tag,'_DEG.tsv'))
    rank <- data.frame(gene_id=tt$gene_id,symbol=tt$symbol,
                       rank=sign(tt$logFC)*sqrt(pmax(tt$F,0)),FDR=tt$FDR)
    rank <- rank[order(-rank$rank),];write(rank,paste0(tag,'_rank.tsv'))
    lcpm <- cpm(y,log=TRUE,prior.count=2)
    write(data.frame(gene_id=rownames(lcpm),lcpm,check.names=FALSE),paste0(tag,'_logCPM.tsv'))
    write(data.frame(sample=samples,cells=m$cells,y$samples,check.names=FALSE),paste0(tag,'_normalization.tsv'))
    saveRDS(list(y=y,design=design,fit=fit,test=test),file.path(out,'models',paste0(tag,'.rds')))
    summary[[tag]] <- data.frame(set=setname,arm=arm,status=status,min_cells=min(m$cells),
       tested=nrow(tt),FDR05=sum(tt$FDR<.05),Up=sum(tt$FDR<.05 & tt$logFC>0),Down=sum(tt$FDR<.05 & tt$logFC<0))
    z <- tt[tt$gene_id=='ENSMUSG00000031928',,drop=FALSE]
    if(nrow(z)) target[[tag]] <- data.frame(set=setname,arm=arm,status=status,tested=TRUE,
                 logFC=z$logFC,PValue=z$PValue,FDR=z$FDR)
    else target[[tag]] <- data.frame(set=setname,arm=arm,status=status,tested=FALSE,logFC=NA,PValue=NA,FDR=NA)
    resultlist[[tag]] <- tt
    # Each point in PCA is one sample; top 2000 variable tested genes.
    vv <- apply(lcpm,1,var); sel <- order(vv,decreasing=TRUE)[seq_len(min(2000,sum(vv>0)))]
    pc <- prcomp(t(lcpm[sel,,drop=FALSE]),center=TRUE,scale.=FALSE)
    vpc <- 100*pc$sdev^2/sum(pc$sdev^2)
    pd <- data.frame(PC1=pc$x[,1],PC2=pc$x[,2],group=group,label=c('Sal1','Sal2','Sal3','Si1','Si2','Si3'))
    p <- ggplot(pd,aes(PC1,PC2,color=group,label=label))+geom_point(size=3)+
      geom_text_repel(seed=42,max.overlaps=Inf)+scale_color_manual(values=c(Saline='#0072B2',Silica='#D55E00'))+
      labs(x=sprintf('PC1 (%.1f%%)',vpc[1]),y=sprintf('PC2 (%.1f%%)',vpc[2]))+theme_classic(base_size=12)
    plot_save(p,paste0(tag,'_PCA'),paste('Sample PCA:',tag,'. Sal1-3=Saline_1-3; Si1-3=Silica_1-3. PCA uses the 2000 most variable tested genes (or all if fewer), TMM logCPM; no sample correction.'))
    tt$Signif <- ifelse(tt$FDR<.05,tt$direction,'NS')
    tt$Label <- ''; top <- seq_len(min(10,nrow(tt)))
    lab <- tt$symbol[top];empty <- is.na(lab)|lab=='';lab[empty]<-tt$gene_id[top][empty]
    full <- lab;long <- nchar(lab)>10;lab[long] <- sprintf('G%02d',top[long]);tt$Label[top]<-lab
    write(data.frame(label=lab,gene_id=tt$gene_id[top],symbol=full,FDR=tt$FDR[top]),paste0(tag,'_Top10_labels.tsv'))
    p <- ggplot(tt,aes(logFC,-log10(pmax(FDR,1e-300)),color=Signif))+geom_point(size=1,alpha=.7)+
      geom_hline(yintercept=-log10(.05),linetype=2,color='grey50')+
      geom_text_repel(aes(label=Label),seed=42,max.overlaps=Inf,size=3,min.segment.length=0)+
      scale_color_manual(values=c(Up='#D55E00',Down='#0072B2',NS='grey75'))+
      labs(x='log2FC',y='-log10 FDR')+theme_classic(base_size=12)
    plot_save(p,paste0(tag,'_volcano'),paste('Silica versus Saline:',tag,'. Positive log2FC denotes higher expression in Silica. Labels identify the ten genes with smallest P values (FDR ordered with P-value tie breaking); labels do not imply significance. NS=FDR>=0.05. Long gene labels:',paste(paste(lab,full,sep='='),collapse='; ')))
    message('MODEL_PASS=',tag)
  }
}
write(do.call(rbind,summary),'DEG_summary.tsv');write(do.call(rbind,target),'Mre11a_DEG.tsv')
# Concordance for matched tested genes, preserving separate model FDR values.
for(setname in c('Core','Naive','ActCD4','Treg')){
 a <- resultlist[[paste0(setname,'_Full')]];b <- resultlist[[paste0(setname,'_NoBG')]]
 j <- merge(a[,c('gene_id','symbol','logFC','FDR')],b[,c('gene_id','logFC','FDR')],by='gene_id',suffixes=c('_Full','_NoBG'))
 j$same_direction <- sign(j$logFC_Full)==sign(j$logFC_NoBG)
 j$both_FDR05 <- j$FDR_Full<.05 & j$FDR_NoBG<.05
 write(j,paste0(setname,'_sensitivity.tsv'))
}
writeLines(c('Contrast: Silica versus Saline; six sample-level replicates; design ~group.',
 'GEO GSE314943 overall design states n=3 per group. Sample GSM identities are recorded in design.tsv. No additional batch covariates are inferred.',
 'Core, Naive and ActCD4 retain all 3+3 samples and at least 20 cells per sample. Treg retains all six samples, including Saline_3 with 13 cells, as a separately flagged exploratory fit.',
 'Full: all genes enter expression filtering. NoBG: a predeclared background panel and hemoglobin genes are removed before expression filtering and TMM normalization. NoBG is normalization/feature sensitivity, not ambient RNA decontamination.',
 'edgeR: filterByExpr, TMM, robust dispersion, robust quasi-likelihood F test. FDR is BH within each set/arm across tested genes. Cross-model FDR is not pooled.',
 'Core DEG includes CD4 subtype composition effects. Subtype DEG supports assessment of within-state changes. State composition is not regressed out of the six-sample Core model.',
 'GSEA ranking is sign(logFC)*sqrt(QLF statistic) over all tested genes; no significant-gene preselection. Gene IDs remain primary; symbol deduplication is deferred to pathway mapping.',
 'Mre11a is tested only if it passes the same expression filter as other genes; failure is reported without forcing significance.'),file.path(out,'tables/methods.txt'))
capture.output(sessionInfo(),file=file.path(out,'tables/sessionInfo.txt'))
message('GSE314943_STAGE6_ANALYSIS_PASS')
