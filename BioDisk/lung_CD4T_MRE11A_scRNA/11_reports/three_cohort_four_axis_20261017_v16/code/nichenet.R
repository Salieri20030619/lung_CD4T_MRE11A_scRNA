suppressPackageStartupMessages(library(nichenetr))
p <- '/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA';o<-file.path(p,'11_reports/three_cohort_four_axis_20261017_v16')
rd<-function(n)read.delim(file.path(o,n),stringsAsFactors=FALSE,check.names=FALSE)
w<-function(x,n)write.table(as.data.frame(x),file=file.path(o,'tables',n),sep='\t',row.names=FALSE,quote=FALSE,na='NA',fileEncoding='UTF-8')
canon<-function(x){x<-toupper(x);a<-c(MRE11A='MRE11',ATP5A1='ATP5F1A',ATP5B='ATP5F1B',ATP5C1='ATP5F1C',ATP5G1='ATP5MC1',ATP5O='ATP5PO',PARK2='PRKN');i<-x%in%names(a);x[i]<-a[x[i]];x}
lt<-readRDS(file.path(p,'00_raw/nichenet_v2_mouse/ligand_target_matrix_nsga2r_final_mouse.rds'));rownames(lt)<-canon(rownames(lt));colnames(lt)<-canon(colnames(lt));stopifnot(!anyDuplicated(rownames(lt)),!anyDuplicated(colnames(lt)))
lr<-rd('tables/expressed_curated_LR.tsv');bg<-rd('inputs/receiver_background.tsv');de<-rd('tables/DE_all_models.tsv');de<-de[de$model%in%c('314_disease','joint_stratified','274_time') & de$contrast!='time_omnibus',]
acts<-list();links<-list();status<-list()
for(fam in unique(de$family)){
 d<-de[de$family==fam,];ds<-if(d$model[1]=='274_time')'GSE274372' else 'GSE314943';sc<-d$scope[1];b<-intersect(bg$gene[bg$dataset==ds & bg$scope==sc],rownames(lt));l<-lr[lr$dataset==ds & lr$scope==sc & lr$sender=='AllMac',];potential<-intersect(unique(l$ligand_parts[!grepl(';',l$ligand_parts)]),colnames(lt))
 for(dir in c('up','down')){
  gs<-intersect(d$gene[d$PValue<.05 & (if(dir=='up') d$logFC>0 else d$logFC<0)],b);base<-data.frame(family=fam,model=d$model[1],dataset=ds,scope=sc,direction=dir,n_geneset=length(gs),n_background=length(b),n_potential=length(potential))
  if(length(gs)<5 || length(gs)>=length(b) || length(potential)==0){base$status<-'insufficient_gene_set_or_ligands';status[[length(status)+1]]<-base;next}
  a<-as.data.frame(predict_ligand_activities(geneset=gs,background_expressed_genes=b,ligand_target_matrix=lt,potential_ligands=potential));a<-a[order(-a$aupr_corrected,a$test_ligand),];a$activity_rank<-seq_len(nrow(a));a$selected<-a$aupr_corrected>0 & a$activity_rank<=30
  for(k in names(base))a[[k]]<-base[[k]];acts[[length(acts)+1]]<-a
  for(lig in a$test_ligand[a$selected]){
   weights<-lt[b,lig];ord<-order(-weights,names(weights));ord<-head(ord[weights[ord]>0],250);t<-data.frame(gene=b[ord],ligand=lig,weight=as.numeric(weights[ord]),target_rank=seq_along(ord),aupr_corrected=a$aupr_corrected[match(lig,a$test_ligand)],activity_rank=a$activity_rank[match(lig,a$test_ligand)]);for(k in names(base))t[[k]]<-base[[k]];links[[length(links)+1]]<-t
  }
  base$status<-'success';status[[length(status)+1]]<-base;cat(fam,dir,length(gs),length(potential),sum(a$selected),'\n')
 }
}
w(do.call(rbind,acts),'NicheNet_activity_all.tsv');w(do.call(rbind,links),'NicheNet_top250_targets.tsv');w(do.call(rbind,status),'NicheNet_run_inventory.tsv');writeLines(capture.output(sessionInfo()),file.path(o,'logs','R_session_NicheNet.txt'));writeLines('success',file.path(o,'logs','NicheNet_complete.txt'))