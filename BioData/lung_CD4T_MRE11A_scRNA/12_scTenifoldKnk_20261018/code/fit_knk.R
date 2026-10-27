suppressPackageStartupMessages({library(scTenifoldKnk);library(Matrix);library(jsonlite)})
O <- '/mnt/BioData/bio/wangmingyuan/lung_CD4T_MRE11A_scRNA/12_scTenifoldKnk_20261018'
jobs <- fromJSON(file.path(O,'jobs.json'),simplifyVector=FALSE)
a <- commandArgs(TRUE); j <- jobs[[as.integer(a[1])]]
d <- file.path(O,'runs',j$id);dir.create(d,recursive=TRUE,showWarnings=FALSE)
writeLines(capture.output(sessionInfo()),file.path(d,'sessionInfo.txt'),useBytes=TRUE)
X <- as(readMM(file.path(j$path,'counts.mtx')),'CsparseMatrix');rownames(X)<-readLines(file.path(j$path,'genes.txt'));colnames(X)<-read.delim(file.path(j$path,'cells.tsv'))$cell
stopifnot(all(is.finite(X@x)),ncol(X)==j$cells,nrow(X)==j$genes,all(Matrix::rowSums(X)>0),all(Matrix::colSums(X)>0))
targets <- unlist(j$targets); warns<-character(); t0<-Sys.time()
result <- tryCatch(withCallingHandlers({
 first<-targets[1]
 out<-scTenifoldKnk(X,gKO=first,qc=FALSE,nc_nNet=10,nc_nCells=min(500L,as.integer(floor(.8*ncol(X)))),nc_nComp=3,td_K=3,ma_method='manifold',dr_direction=FALSE,nCores=4,seed=j$seed)
 WT<-as.matrix(out$tensorNetworks$WT)
 saveRDS(out$tensorNetworks$WT,file.path(d,'WT.rds'))
 dirs<-knockoutDirection(X,gKO=as.list(targets),genes=rownames(WT),regressLibSize=TRUE)
 audits<-list()
 for(g in targets){
  edges<-sum(WT[g,]!=0)
  if(edges==0){audits[[g]]<-list(gene=g,status='NO_OUTGOING_EDGES',outgoing_edges=0);next}
  if(g==first)DR<-out$diffRegulation else {KO<-WT;KO[g,]<-0;set.seed(j$seed);ma<-scTenifoldNet::manifoldAlignment(WT,KO,d=2,nCores=4);DR<-scTenifoldKnk::dRegulation(ma,gKO=g)}
  stopifnot(!anyDuplicated(DR$gene),all(is.finite(DR$distance)),all(is.finite(DR$p.value)),all(DR$p.value>=0&DR$p.value<=1))
  DR$directionScore<-as.numeric(dirs[g,DR$gene]);DR$direction<-ifelse(DR$directionScore>0,'up',ifelse(DR$directionScore<0,'down','zero'))
  DR$is_knockout_target<-DR$gene==g
  write.table(DR,file.path(d,paste0(g,'.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
  audits[[g]]<-list(gene=g,status='success',outgoing_edges=edges,genes=nrow(DR),nominal_nonself=sum(DR$p.value<.05&DR$gene!=g),fdr_nonself=sum(DR$p.adj<.05&DR$gene!=g))
 }
 list(status='success',targets=audits)
},warning=function(w){warns<<-c(warns,conditionMessage(w))}),error=function(e){list(status='error',message=conditionMessage(e))})
result$warnings<-unique(warns);result$elapsed_seconds<-as.numeric(difftime(Sys.time(),t0,units='secs'));result$job<-j
write_json(result,file.path(d,'complete.json'),pretty=TRUE,auto_unbox=TRUE)
print(result)
if(result$status!='success')quit(status=1)