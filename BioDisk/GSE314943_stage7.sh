#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA
SRC="$ROOT/02_analysis/GSE314943/stage6_20260905_231807_L44XdT"
OUT=$(mktemp -d "$ROOT/02_analysis/GSE314943/stage7_$(date +%Y%m%d_%H%M%S)_XXXXXX")
trap 'rc=$?; if ((rc!=0)); then echo "GSE314943_STAGE7_FAIL EXIT_CODE=$rc"; fi' EXIT
cat > "$OUT/analysis.R" <<'RSCRIPT'
args<-commandArgs(TRUE);src<-args[1];out<-args[2]
req<-c('AnnotationDbi','org.Mm.eg.db','GO.db','fgsea','ggplot2')
miss<-req[!vapply(req,requireNamespace,logical(1),quietly=TRUE)]
if(length(miss)){
 message('INSTALL_MISSING=',paste(miss,collapse=','))
 options(repos=c(CRAN='https://cloud.r-project.org'),timeout=600)
 if(!requireNamespace('BiocManager',quietly=TRUE))install.packages('BiocManager')
 BiocManager::install(miss,ask=FALSE,update=FALSE)
 miss<-req[!vapply(req,requireNamespace,logical(1),quietly=TRUE)]
 if(length(miss))stop('Missing R packages after installation: ',paste(miss,collapse=', '))
}
suppressPackageStartupMessages({library(AnnotationDbi);library(org.Mm.eg.db);library(GO.db);library(fgsea);library(ggplot2)})
options(timeout=300);set.seed(42)
for(d in c('tables','figures','database'))dir.create(file.path(out,d),showWarnings=FALSE)
wr<-function(x,n)write.table(x,file.path(out,'tables',n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
rd<-function(n)read.delim(file.path(src,'tables',n),check.names=FALSE,stringsAsFactors=FALSE)
sets<-c('Core','Naive','ActCD4','Treg');arms<-c('Full','NoBG')
allids<-unique(unlist(lapply(sets,function(s)rd(paste0(s,'_Full_DEG.tsv'))$gene_id)))
valid<-intersect(allids,keys(org.Mm.eg.db,keytype='ENSEMBL'))
mp<-unique(AnnotationDbi::select(org.Mm.eg.db,keys=valid,keytype='ENSEMBL',columns='ENTREZID'))
mp<-mp[!is.na(mp$ENTREZID),];mp$ENTREZID<-as.character(mp$ENTREZID)
amb<-names(which(table(mp$ENSEMBL)>1));wr(mp,'ID_mapping_all.tsv')
mp<-mp[!mp$ENSEMBL %in% amb,];wr(mp,'ID_mapping_unambiguous.tsv')
go<-AnnotationDbi::select(org.Mm.eg.db,keys=unique(mp$ENTREZID),keytype='ENTREZID',columns=c('GOALL','ONTOLOGYALL'))
go<-unique(go[!is.na(go$GOALL)&go$ONTOLOGYALL=='BP',c('ENTREZID','GOALL')])
gosets<-lapply(split(as.character(go$ENTREZID),go$GOALL),unique)
gonames<-AnnotationDbi::select(GO.db,keys=names(gosets),keytype='GOID',columns='TERM')
# Two public KEGG REST requests; retain the exact responses used in this run.
fetch<-function(url,name){
 dest<-file.path(out,'database',name)
 for(attempt in 1:3){
  ok<-tryCatch({code<-download.file(url,paste0(dest,'.part'),method='libcurl',quiet=TRUE,mode='wb');isTRUE(code==0)&&file.info(paste0(dest,'.part'))$size>100},error=function(e)FALSE)
  if(isTRUE(ok)){if(!file.rename(paste0(dest,'.part'),dest))stop('Cannot commit download');return(dest)}
 }
 stop('Download failed: ',url)
}
kf<-fetch('https://rest.kegg.jp/link/pathway/mmu','KEGG_links.tsv')
kn<-fetch('https://rest.kegg.jp/list/pathway/mmu','KEGG_names.tsv')
k<-read.delim(kf,header=FALSE,stringsAsFactors=FALSE);stopifnot(ncol(k)==2,all(grepl('^mmu:',k[[1]])))
k[[1]]<-sub('^mmu:','',k[[1]]);k[[2]]<-sub('^path:','',k[[2]])
knm<-read.delim(kn,header=FALSE,stringsAsFactors=FALSE);knm[[1]]<-sub('^path:','',knm[[1]])
ksets<-lapply(split(k[[1]],k[[2]]),unique)
pws<-list(GO_BP=gosets,KEGG=ksets)
desc<-c(setNames(gonames$TERM,gonames$GOID),setNames(knm[[2]],knm[[1]]))
saveRDS(pws,file.path(out,'database/pathways.rds'))
wr(data.frame(pathway=names(desc),description=unname(desc)),'pathway_names.tsv')
nonrp<-function(x)!grepl('ribosom|translation|translational',x,ignore.case=TRUE)
# Genuine weighted running enrichment, hit positions and ranked metric.
curves<-function(res,stats,pw,stem){
 z<-res[is.finite(res$padj)&res$padj<.05 & nonrp(res$description),,drop=FALSE]
 z<-head(z[order(z$padj,z$pval,z$pathway),],5)
 if(!nrow(z)){writeLines('No FDR<0.05 NonRP pathway.',file.path(out,'figures',paste0(stem,'_no_curve.txt')));return(invisible(NULL))}
 z$label<-sprintf('P%02d',seq_len(nrow(z)))
 wr(z[,setdiff(names(z),'leadingEdge'),drop=FALSE],paste0(stem,'_labels.tsv'))
 n<-length(stats);cols<-c('#0072B2','#D55E00','#009E73','#CC79A7','#E69F00')[seq_len(nrow(z))]
 hit<-lapply(z$pathway,function(id)which(names(stats)%in%pw[[id]]))
 es<-lapply(hit,function(h){
   inc<-rep(-1/(n-length(h)),n);w<-abs(stats[h]);stopifnot(sum(w)>0)
   inc[h]<-w/sum(w);cumsum(inc)
 })
 observed<-vapply(es,function(v){if(max(v)>-min(v))max(v) else min(v)},numeric(1))
 stopifnot(all(abs(observed-z$ES)<1e-5))
 draw<-function(){
  layout(matrix(1:3,3,1),heights=c(3,1.1,1.4))
  par(mar=c(.5,4,1,1))
  plot(1:n,es[[1]],type='n',ylim=range(c(0,unlist(es))),xlab='',ylab='Running ES',xaxt='n')
  abline(h=0,col='grey75');for(i in seq_along(es))lines(1:n,es[[i]],col=cols[i],lwd=1.7)
  legend('topright',legend=z$label,col=cols,lty=1,lwd=2,bty='n',cex=.8)
  par(mar=c(.5,4,0,1));plot(c(1,n),c(.5,nrow(z)+.5),type='n',axes=FALSE,xlab='',ylab='Hits')
  for(i in seq_along(hit))segments(hit[[i]],i-.35,hit[[i]],i+.35,col=cols[i])
  axis(2,at=seq_len(nrow(z)),labels=z$label,las=1,cex.axis=.7)
  par(mar=c(4,4,0,1));plot(seq_len(n),stats,type='l',col='grey35',xlab='Gene rank',ylab='Rank score');abline(h=0,col='grey75')
 }
 for(ext in c('png','pdf')){
  f<-file.path(out,'figures',paste0(stem,'.',ext))
  if(ext=='png')png(f,width=8,height=7,units='in',res=350) else pdf(f,width=8,height=7)
  tryCatch(draw(),finally=dev.off())
 }
 writeLines(c('Top five FDR<0.05 NonRP pathways, selected by adjusted P value then P value. NonRP excludes pathway names containing ribosome/translation terms for display only; all pathways remain in the result table.',
 'Upper: weighted running ES (weight exponent 1). Middle: separate hit tracks. Lower: complete ranked statistic. Positive rank indicates higher expression in Silica. Curves are checked against fgsea ES.',
 paste(z$label,z$pathway,z$description,paste0('NES=',signif(z$NES,4)),paste0('FDR=',signif(z$padj,4)),sep=' | ')),file.path(out,'figures',paste0(stem,'_caption.txt')))
}
summ<-list();aud<-list();counter<-0L
for(ss in sets)for(aa in arms){
 tag<-paste(ss,aa,sep='_');tt<-rd(paste0(tag,'_DEG.tsv'));lc<-rd(paste0(tag,'_logCPM.tsv'))
 tt$avglog<-rowMeans(as.matrix(lc[match(tt$gene_id,lc$gene_id),-1,drop=FALSE]))
 tt$ENTREZID<-mp$ENTREZID[match(tt$gene_id,mp$ENSEMBL)]
 tt<-tt[!is.na(tt$ENTREZID),];tt<-tt[order(-tt$avglog,tt$gene_id),];tt<-tt[!duplicated(tt$ENTREZID),]
 tt$rank<-sign(tt$logFC)*sqrt(pmax(tt$F,0))
 tt<-tt[is.finite(tt$rank),];tt<-tt[order(-tt$rank,tt$ENTREZID),]
 wr(tt[,c('gene_id','symbol','ENTREZID','rank','avglog','FDR')],paste0(tag,'_mapped.tsv'))
 aud[[tag]]<-data.frame(model=tag,mapped_genes=nrow(tt),ties=sum(duplicated(tt$rank)))
 # NoIg is rank-set sensitivity, not a newly fitted DEG model.
 variants<-if(aa=='NoBG')c('Standard','NoIg') else 'Standard'
 for(vv in variants){
  x<-tt
  if(vv=='NoIg')x<-x[!grepl('^(Igh|Igl|Igk)',x$symbol)&x$symbol!='Jchain',]
  stats<-setNames(x$rank,x$ENTREZID);stopifnot(length(stats)>100,!anyDuplicated(names(stats)))
  for(db in names(pws)){
   stem<-paste(tag,vv,db,sep='_');message('START=',stem)
   pw<-lapply(pws[[db]],intersect,y=names(stats));pw<-pw[lengths(pw)>=15 & lengths(pw)<=500]
   set.seed(42)
   f<-as.data.frame(fgseaMultilevel(pw,stats,minSize=15,maxSize=500,eps=0,nproc=2))
   f$description<-unname(desc[f$pathway]);f$description[is.na(f$description)]<-f$pathway[is.na(f$description)]
   f<-f[order(f$padj,f$pval,na.last=TRUE),]
   curves(f,stats,pw,paste0(stem,'_GSEA'))
   f$leadingEdge<-vapply(f$leadingEdge,paste,character(1),collapse=';')
   wr(f,paste0(stem,'_GSEA.tsv'))
   counter<-counter+1L;summ[[counter]]<-data.frame(model=tag,variant=vv,db=db,tested=nrow(f),FDR05=sum(f$padj<.05,na.rm=TRUE),status=if(ss=='Treg')'exploratory_low_cells' else 'primary')
   if(vv!='Standard')next
   # ORA: separate up/down FDR<0.05 lists and database-annotated tested universe.
   universe<-intersect(names(stats),unique(unlist(pws[[db]])))
   opw<-lapply(pws[[db]],intersect,y=universe);opw<-opw[lengths(opw)>=15 & lengths(opw)<=500]
   for(direction in c('Up','Down')){
    selected<-intersect(x$ENTREZID[x$FDR<.05 & if(direction=='Up')x$logFC>0 else x$logFC<0],universe)
    o<-data.frame(pathway=names(opw),size=lengths(opw),hits=vapply(opw,function(g)sum(g%in%selected),integer(1)))
    o$selected_n<-length(selected);o$universe_n<-length(universe)
    o$PValue<-phyper(o$hits-1,o$size,length(universe)-o$size,length(selected),lower.tail=FALSE)
    o$FDR<-p.adjust(o$PValue,'BH');o$description<-unname(desc[o$pathway]);o<-o[order(o$FDR,o$PValue),]
    oname<-paste0(stem,'_ORA_',direction);wr(o,paste0(oname,'.tsv'))
    z<-head(o[o$FDR<.05 & o$hits>0 & nonrp(o$description),],10)
    if(!nrow(z)){writeLines('No FDR<0.05 NonRP ORA term.',file.path(out,'figures',paste0(oname,'_no_plot.txt')));next}
    z$label<-sprintf('P%02d',seq_len(nrow(z)));z$ratio<-z$hits/z$selected_n;wr(z,paste0(oname,'_labels.tsv'))
    p<-ggplot(z,aes(ratio,factor(label,levels=rev(label)),size=hits,color=-log10(pmax(FDR,1e-300))))+geom_point()+
      scale_color_gradient(low='#0072B2',high='#D55E00')+labs(x='Gene ratio',y=NULL,size='Hits',color='-log10 FDR')+theme_classic(base_size=12)
    for(ext in c('png','pdf'))ggsave(file.path(out,'figures',paste0(oname,'.',ext)),p,width=6,height=5,dpi=350,bg='white')
    writeLines(c(paste('ORA',direction,stem,'; BH FDR<0.05 NonRP terms. Gene ratio=overlap/selected genes. Each model uses its own mapped, tested, database-annotated universe.'),paste(z$label,z$description,sep=' = ')),file.path(out,'figures',paste0(oname,'_caption.txt')))
   }
  }
 }
}
wr(do.call(rbind,summ),'GSEA_summary.tsv');wr(do.call(rbind,aud),'mapping_summary.tsv')
writeLines(c(paste('Input:',src),paste('Generated:',Sys.time()),
 'GO BP uses GOALL ancestor annotations from installed org.Mm.eg.db and GO.db; KEGG uses mouse REST pathways downloaded and retained in database/.',
 'Ambiguous one-Ensembl-to-many-Entrez mappings are dropped. Where multiple Ensembl genes map to one Entrez, keep the greatest mean logCPM across all six samples; ties break by Ensembl ID, without selecting on DE significance.',
 'GSEA uses all mapped tested genes and sign(logFC)*sqrt(QLF). Tied ranks use deterministic Entrez order, no artificial jitter. FDR is BH within each model/variant/database.',
 'NoIg drops Igh/Igl/Igk-prefixed genes and Jchain from NoBG ranks and pathway overlaps. NoIg is not decontamination and does not refit edgeR or provide new gene-level P values.',
 'ORA uses FDR<0.05 DEG separately by direction. Universe is mapped tested genes annotated to that database. All size-eligible pathways, including zero-hit terms, enter BH.',
 'Gene-set size is 15-500 after universe intersection. NonRP only selects display terms; full tables retain ribosome/translation terms.',
 'GSEA tests ranked gene-set enrichment; its significance does not increase the biological sample size. Treg remains exploratory due to 13 cells in one control sample.',
 'The Mre11a function hypothesis will be evaluated separately; pathway enrichment alone does not identify its causal role.'),file.path(out,'tables/methods.txt'))
capture.output(sessionInfo(),file=file.path(out,'tables/sessionInfo.txt'))
message('GSE314943_STAGE7_ANALYSIS_PASS')
RSCRIPT
conda run -n R_4.5.2 Rscript -e 'invisible(parse(file=commandArgs(TRUE)[1]));cat("R_PARSE_PASS\n")' "$OUT/analysis.R" > "$OUT/syntax.log" 2>&1
echo 'SYNTAX_CHECK=PASS'
echo "OUTPUT=$OUT"
conda run --no-capture-output -n R_4.5.2 Rscript "$OUT/analysis.R" "$SRC" "$OUT"
tar -C "$OUT" -czf "$ROOT/tmp/GSE314943_stage7_return.tar.gz.part" tables figures database analysis.R syntax.log
mv "$ROOT/tmp/GSE314943_stage7_return.tar.gz.part" "$ROOT/tmp/GSE314943_stage7_return.tar.gz"
echo 'GSE314943_STAGE7_PASS'
