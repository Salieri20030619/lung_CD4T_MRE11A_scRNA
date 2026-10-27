suppressPackageStartupMessages(library(nichenetr))
P<-'/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA';O<-file.path(P,'11_reports/death_priority_20261018_v2');V<-file.path(P,'11_reports/death_priority_20261018_v31')
rd<-function(n)read.delim(file.path(O,'tables',n),stringsAsFactors=FALSE,check.names=FALSE)
wr<-function(x,n)write.table(as.data.frame(x),file.path(O,'tables',n),sep='\t',quote=FALSE,row.names=FALSE,na='NA',fileEncoding='UTF-8')
canon<-function(x){x<-toupper(x);a<-c(MRE11A='MRE11',ATP5A1='ATP5F1A',ATP5B='ATP5F1B',ATP5C1='ATP5F1C',ATP5G1='ATP5MC1',ATP5O='ATP5PO',PARK2='PRKN');i<-x%in%names(a);x[i]<-a[x[i]];x}
lt<-readRDS(file.path(P,'00_raw/nichenet_v2_mouse/ligand_target_matrix_nsga2r_final_mouse.rds'));rownames(lt)<-canon(rownames(lt));colnames(lt)<-canon(colnames(lt))
stopifnot(!anyDuplicated(rownames(lt)),!anyDuplicated(colnames(lt)))
lr<-rd('joint_LR_eligibility.tsv');lr$eligible<-tolower(as.character(lr$eligible))=='true';bg<-rd('joint_receiver_background.tsv');cand<-rd('frozen_candidate_pairs.tsv')
traj<-read.delim(gzfile(file.path(V,'tables/gene_trajectory_all.tsv.gz')));traj<-traj[traj$variant=='K15_R1',]
activities<-list();links<-list();runs<-list()
for(source in c('AllMac','MacStrict'))for(sc in c('Core','Naive','ActCD4','Treg')){
 b<-intersect(bg$gene[bg$scope==sc],rownames(lt));l<-lr[lr$source==source&lr$scope==sc&lr$eligible,];pot<-intersect(unique(l$ligand_parts[!grepl(';',l$ligand_parts)]),colnames(lt))
 de<-read.delim(file.path(O,'inputs',paste0('joint_DE_',sc,'.tsv')))
 for(response in c('joint_disease_DE','joint_state_trajectory')){
  if(response=='joint_state_trajectory'&sc=='Core')next
  d<-if(response=='joint_disease_DE')data.frame(gene=de$gene,P=de$PValue,effect=de$logFC) else data.frame(gene=traj$gene[traj$scope==sc],P=traj$P[traj$scope==sc],effect=traj$mean_slope[traj$scope==sc])
  for(direction in c('up','down')){
   gs<-intersect(d$gene[d$P<.05&if(direction=='up')d$effect>0 else d$effect<0],b)
   base<-data.frame(source=source,scope=sc,response=response,direction=direction,n_geneset=length(gs),n_background=length(b),n_ligands=length(pot),n_donors=10)
   if(length(gs)<5||length(gs)>=length(b)||!length(pot)){base$status<-'NOT_EVALUABLE';runs[[length(runs)+1]]<-base;next}
   a<-as.data.frame(predict_ligand_activities(gs,b,lt,pot));a<-a[order(-a$aupr_corrected,a$test_ligand),];a$activity_rank<-seq_len(nrow(a));a$selected<-a$aupr_corrected>0&a$activity_rank<=30
   for(k in names(base))a[[k]]<-base[[k]];activities[[length(activities)+1]]<-a
   targets<-unique(cand$gene[if(sc=='Core')rep(TRUE,nrow(cand)) else cand$scope==sc])
   for(lig in a$test_ligand[a$selected]){
    weights<-lt[b,lig];ord<-order(-weights,names(weights));ranks<-setNames(seq_along(ord),b[ord])
    for(g in targets){
     t<-data.frame(gene=g,ligand=lig,weight=if(g%in%b)as.numeric(lt[g,lig]) else NA,target_rank=if(g%in%b)unname(ranks[g]) else NA,in_response=g%in%gs,in_background=g%in%b,activity_rank=a$activity_rank[match(lig,a$test_ligand)],aupr_corrected=a$aupr_corrected[match(lig,a$test_ligand)])
     t$prior_top250<-is.finite(t$weight)&t$weight>0&is.finite(t$target_rank)&t$target_rank<=250
     for(k in names(base))t[[k]]<-base[[k]];links[[length(links)+1]]<-t
    }
   }
   base$status<-'success';runs[[length(runs)+1]]<-base
   cat('NICHENET',source,sc,response,direction,length(gs),length(pot),'\n')
  }
 }
}
wr(do.call(rbind,activities),'joint_NicheNet_activity.tsv');wr(do.call(rbind,links),'joint_NicheNet_candidate_links.tsv');wr(do.call(rbind,runs),'joint_NicheNet_inventory.tsv');writeLines('success',file.path(O,'logs/NicheNet_complete.txt'))
