suppressPackageStartupMessages(library(Matrix))
P <- '/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA';O <- file.path(P,'11_reports/phenotype_candidate_discovery_20261005_v1');H<-file.path(P,'11_reports/human_first_20261004_v1/results');R<-file.path(O,'results')
rd<-function(p)read.delim(p,check.names=FALSE,stringsAsFactors=FALSE)
wr<-function(d,n)write.table(d,file.path(R,n),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
L<-c('ALCAM','HBEGF','HMGB1');res<-'/mnt/BioDisk/wangmingyuan/IL17_Candida_integrated_v1/07_synthesis/supplementary_bioinformatics_v1/a6_consensus_regulatory_chain_v1/01_nichenet_resources_v1'
lr<-readRDS(file.path(res,'lr_network.rds'));wr(as.data.frame(lr[lr$from%in%L|lr$from=='NUCB2',]),'selected_LR_prior_original_provenance.tsv')
lt<-readRDS(file.path(res,'ligand_target_matrix.rds'));members<-rd(file.path(H,'stage1_complete_program_members.tsv'));gs<-unique(members$gene[members$program%in%c('Mitophagy','TReg')]);q<-expand.grid(gene=gs,ligand=L,stringsAsFactors=FALSE);q$weight<-mapply(function(g,l)if(g%in%rownames(lt)&&l%in%colnames(lt))lt[g,l]else NA_real_,q$gene,q$ligand);wr(q,'selected_full_program_prior_weights.tsv');rm(lt);gc()
pre<-file.path(P,'02_analysis/scenic3_20261001_123452_Kvp9h5/input/human_mac_');m<-rd(paste0(pre,'cells.tsv'));g<-rd(paste0(pre,'genes.tsv'))$gene;x<-as(readMM(paste0(pre,'counts.mtx')),'CsparseMatrix');stopifnot(nrow(x)==length(g),ncol(x)==nrow(m));rownames(x)<-g
rev<-rd(file.path(H,'stage1_sender_cell_dictionary.tsv.gz'));rev<-rev[match(m$cell,rev$cell),];stopifnot(all(rev$cell==m$cell));des<-rd(file.path(H,'stage1_donor_design.tsv'));m$label<-rev$label;m$review_scope<-rev$review_scope;m$cluster_review<-rev$r4
keep<-m$sample%in%des$sample;gs<-intersect(c(L,'C1QA','C1QB','C1QC','FABP4','MARCO','PPARG','MRC1','TREM2','GPNMB','FCN1','S100A8','S100A9','EPCAM','CD3D'),g);xx<-x[gs,,drop=FALSE];mature<-intersect(c('C1QA','C1QB','C1QC','MRC1','TREM2','GPNMB'),g);mono<-intersect(c('FCN1','S100A8','S100A9'),g);m$mac_marker2<-colSums(x[mature,,drop=FALSE]>0)>=2;m$mono_marker2<-colSums(x[mono,,drop=FALSE]>0)>=2
out<-list();contexts<-list()
for(ss in des$sample)for(cl in sort(unique(m$cluster_review))){ii<-which(keep&m$sample==ss&m$cluster_review==cl);if(!length(ii))next;n<-length(ii);z<-xx[,ii,drop=FALSE];lab<-unique(m$label[ii]);scope<-unique(m$review_scope[ii]);stopifnot(length(lab)==1,length(scope)==1)
 out[[length(out)+1]]<-data.frame(sample=ss,group=des$group[match(ss,des$sample)],series=des$series[match(ss,des$sample)],cluster=cl,label=lab,scope=scope,cells=n,gene=gs,pct=as.numeric(rowMeans(z>0)),mean_CP10k=as.numeric(rowMeans(z%*%Diagonal(x=1e4/m$library[ii]))))
 for(l in L){positive<-ii[x[l,ii]>0];contexts[[length(contexts)+1]]<-data.frame(sample=ss,cluster=cl,label=lab,scope=scope,gene=l,cells=n,positive=length(positive),mature_marker2_positive=sum(m$mac_marker2[positive]),mono_marker2_positive=sum(m$mono_marker2[positive]))}
}
wr(do.call(rbind,out),'selected_cluster_gene_donors_SERVER_ONLY.tsv');wr(do.call(rbind,contexts),'selected_positive_cell_identity_SERVER_ONLY.tsv')
shared<-file.path(P,'02_analysis/mechanism_step2_20260927_154919_EsN4gk/shared.py');file.copy(shared,file.path(O,'provenance/original_program_definition_shared.py'),overwrite=TRUE)
cat('RAW_CONTEXT_COMPLETE\n');writeLines(capture.output(sessionInfo()),file.path(O,'provenance/raw_context_R_sessionInfo.txt'))
