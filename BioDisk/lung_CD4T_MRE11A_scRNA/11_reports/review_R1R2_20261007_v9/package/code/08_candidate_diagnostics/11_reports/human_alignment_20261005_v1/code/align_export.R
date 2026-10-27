P <- '/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA'
O <- file.path(P,'11_reports/human_alignment_20261005_v1'); H <- file.path(P,'11_reports/human_first_20261004_v1/results')
suppressPackageStartupMessages(library(Matrix))
wr <- function(d,n) {con<-gzfile(file.path(O,'results',n),'wt');write.table(d,con,sep='\t',quote=FALSE,row.names=FALSE,na='NA');close(con)}
z<-readRDS(file.path(H,'stage2_aggregates.rds')); genes<-c('IKZF2','STAT1','IL15','IL15RA','IL2RB','IL2RG','TNF','CLSTN1','ADAM15','PLAU','DSC2','CXCL2','ICAM2','VTN','VEGFB','GSTP1','RTN4')
out<-list();for(key in names(c(z$receiver,z$sender))){a<-c(z$receiver,z$sender)[[key]];for(g in genes){present<-g%in%rownames(a$pb);d<-a$sm;d$key<-key;d$gene<-g;d$present<-present;d$pct<-if(present)a$pct[g,] else NA;d$mean_CP10k<-if(present)a$avg[g,] else NA;d$mean_log1pCP10k<-if(present)a$logavg[g,] else NA;out[[length(out)+1]]<-d}}
wr(do.call(rbind,out),'alignment_gene_donors.tsv.gz')
pre<-file.path(P,'02_analysis/scenic3_20261001_123452_Kvp9h5/input/human_mac_');m<-read.delim(paste0(pre,'cells.tsv'));g<-read.delim(paste0(pre,'genes.tsv'))$gene;x<-as(readMM(paste0(pre,'counts.mtx')),'CsparseMatrix');rownames(x)<-g;x<-x[c('IL15','IL15RA'),,drop=FALSE]
out<-list();for(key in names(z$sender)){a<-z$sender[[key]];for(s in a$sm$sample){i<-which(a$take&m$sample==s);marg<-if(length(i))rowMeans(x[c('IL15','IL15RA'),i,drop=FALSE]>0) else c(NA,NA);both<-if(length(i))mean(colSums(x[c('IL15','IL15RA'),i,drop=FALSE]>0)==2) else NA;out[[length(out)+1]]<-data.frame(sender=key,sample=s,cells=length(i),IL15=marg[1],IL15RA=marg[2],joint=both)}}
wr(do.call(rbind,out),'IL15_sender_joint_coverage.tsv.gz');rm(x)
# Program scores excluding either target TF, identical frozen genes, no new threshold search.
mem<-read.delim(gzfile(file.path(H,'stage2_targetset_members.tsv.gz')));kk<-c('receiver','filter','model','test','direction');out<-list()
for(ids in split(seq_len(nrow(mem)),interaction(mem[kk],drop=TRUE))){v<-mem[ids,];a<-z$receiver[[paste(v$receiver[1],v$filter[1],sep='_')]];for(tf in c('IKZF2','STAT1')){gs<-setdiff(v$gene,tf);gs<-intersect(gs,rownames(a$logavg));d<-a$sm;for(k in kk)d[[k]]<-v[[k]][1];d$excluded_tf<-tf;d$genes_n<-length(gs);d$target_score<-if(length(gs))colMeans(a$logavg[gs,,drop=FALSE]) else NA;out[[length(out)+1]]<-d}}
wr(do.call(rbind,out),'TF_excluded_program_donors.tsv.gz')
sink(file.path(O,'provenance/scTenifoldKnk_installed_API.txt'));cat('VERSION',as.character(packageVersion('scTenifoldKnk')),'\n');ns<-asNamespace('scTenifoldKnk');for(n in getNamespaceExports('scTenifoldKnk')){cat('\nFUNCTION',n,'\n');print(get(n,ns))};sink()
cat('EXPORT_COMPLETE\n')
