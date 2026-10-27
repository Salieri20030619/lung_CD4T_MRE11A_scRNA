P<-'/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA';O<-file.path(P,'11_reports/human_first_20261004_v1');A<-file.path(P,'02_analysis');p1<-file.path(A,'mechanism_step1_20260927_152218_s4V5eF');p3<-file.path(A,'mechanism_step3_20260928_103057_lIRL38');p5<-file.path(A,'mechanism_step5_20260928_153825_Q4sB1x')
suppressPackageStartupMessages(library(Matrix));rd<-function(p)read.delim(p,check.names=FALSE)
g<-rd(file.path(p1,'comm_input/CD4_genes.tsv'));gm<-rd(file.path(p3,'tables/CD4_gene_mapping.tsv'));m<-rd(file.path(p1,'comm_input/CD4_cells.tsv'))
con<-gzfile(file.path(p1,'comm_input/CD4_counts.mtx.gz'),'rt');x<-as(readMM(con),'CsparseMatrix');close(con)
bgp<-strsplit('SFTPA1 SFTPA2 SFTPB SFTPC SFTPD SCGB1A1 SCGB3A1 SCGB3A2 ABCA3 SLC34A2 NAPSA LYZ LYZ1 LYZ2 S100A8 S100A9 CSF1R LST1 TYROBP FCER1G',' ')[[1]]
bgi<-g$canonical%in%bgp|grepl('^(IGH[VDJMAEG]|IGK[VJC]|IGL[VJC])',g$canonical);m$bg<-colSums(x[bgi,,drop=FALSE])/m$library
u<-sort(unique(gm$mapped));x<-sparseMatrix(i=match(gm$mapped,u),j=seq_len(nrow(gm)),x=1,dims=c(length(u),nrow(gm)))%*%x;rownames(x)<-u
lt<-readRDS(file.path(P,'00_raw/nichenet_v2_mouse/ligand_target_matrix_nsga2r_final_mouse.rds'));old<-rd(file.path(p5,'tables/geneset_status.tsv'));rows<-list()
for(rec in c('Core','ActCD4'))for(f in c('All','BG1')){
 take<-if(rec=='Core')rep(TRUE,nrow(m))else m$subtype==rec;if(f=='BG1')take<-take&m$bg<=.01;mm<-m[take,];xx<-x[,take,drop=FALSE];ss<-sort(unique(m$sample));n<-as.numeric(table(factor(mm$sample,levels=ss)));ag<-sparseMatrix(i=seq_len(nrow(mm)),j=match(mm$sample,ss),x=1,dims=c(nrow(mm),length(ss)));pct<-sweep(as.matrix((xx>0)%*%ag),2,n,'/')
 de<-rd(file.path(p3,'tables',paste0('DE_CD4_',rec,'_',f,'.tsv')));ec<-rownames(pct)[rowSums(pct[,grepl('^Saline_',ss)]>=.1)>=2];ed<-rownames(pct)[rowSums(pct[,grepl('^Silica_',ss)]>=.1)>=2];bg<-Reduce(intersect,list(de$gene,union(ec,ed),rownames(lt)))
 o<-old[old$receiver==rec&old$filter==f&old$gate=='Pct10'&old$omit=='Full',];stopifnot(all(o$background==length(bg)));rows[[length(rows)+1]]<-data.frame(receiver=rec,filter=f,gate='Pct10',omit='Full',gene=bg)
}
write.table(do.call(rbind,rows),file.path(O,'results/stage3_mouse_background_members.tsv'),sep='\t',row.names=FALSE,quote=FALSE);cat('MOUSE_BACKGROUND_RECONSTRUCTION_EXACT\n')
