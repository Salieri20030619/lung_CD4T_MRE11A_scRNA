suppressPackageStartupMessages(library(nichenetr))
P<-'/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA'
V<-file.path(P,'11_reports/death_priority_20261018_v2')
O<-'/mnt/BioData/bio/wangmingyuan/lung_CD4T_MRE11A_scRNA/18_macrophage_supplement_20261019'
rd<-function(n)read.delim(file.path(V,'tables',n),stringsAsFactors=FALSE)
wr<-function(x,n)write.table(as.data.frame(x),file.path(O,'tables',n),sep='\t',quote=FALSE,row.names=FALSE,na='NA',fileEncoding='UTF-8')
canon<-function(x){x<-toupper(x);a<-c(MRE11A='MRE11',ATP5A1='ATP5F1A',ATP5B='ATP5F1B',ATP5C1='ATP5F1C',ATP5G1='ATP5MC1',ATP5O='ATP5PO',PARK2='PRKN');i<-x%in%names(a);x[i]<-a[x[i]];x}
lt<-readRDS(file.path(P,'00_raw/nichenet_v2_mouse/ligand_target_matrix_nsga2r_final_mouse.rds'));rownames(lt)<-canon(rownames(lt));colnames(lt)<-canon(colnames(lt))
lr<-rd('joint_LR_eligibility.tsv');lr$eligible<-tolower(as.character(lr$eligible))=='true'
bg<-rd('joint_receiver_background.tsv');mods<-rd('module_members.tsv');mods<-mods[tolower(as.character(mods$retained))=='true',]
traj<-read.delim(gzfile(file.path(P,'11_reports/death_priority_20261018_v31/tables/gene_trajectory_all.tsv.gz')));traj<-traj[traj$variant=='K15_R1',]
activities<-list();runs<-list();members<-list();links<-list()
eps<-c('TAct','TReg','Inhib','Naive','Cycle','MitoROS','MitoISR','ApoExec','PyroExec')
for(source in c('AllMac','MacStrict'))for(sc in c('Naive','ActCD4','Treg')){
 b<-intersect(bg$gene[bg$scope==sc],rownames(lt))
 l<-lr[lr$source==source&lr$scope==sc&lr$eligible,];pot<-intersect(unique(l$ligand_parts[!grepl(';',l$ligand_parts)]),colnames(lt))
 de<-read.delim(file.path(V,'inputs',paste0('joint_DE_',sc,'.tsv')))
 for(response in c('disease','trajectory')){
  d<-if(response=='disease')data.frame(gene=de$gene,P=de$PValue,effect=de$logFC) else data.frame(gene=traj$gene[traj$scope==sc],P=traj$P[traj$scope==sc],effect=traj$mean_slope[traj$scope==sc])
  for(ep in eps)for(direction in c('up','down')){
   gs<-intersect(intersect(d$gene[!is.na(d$P)&d$P<.05&if(direction=='up')d$effect>0 else d$effect<0],mods$gene[mods$module==ep]),b)
   base<-data.frame(source=source,scope=sc,response=response,program=ep,direction=direction,n_geneset=length(gs),n_background=length(b),n_ligands=length(pot),genes=paste(gs,collapse=';'))
   if(length(gs)<5||length(gs)>=length(b)||!length(pot)){base$status<-'INSUFFICIENT_RESPONSE_GENES';runs[[length(runs)+1]]<-base;next}
   a<-as.data.frame(predict_ligand_activities(gs,b,lt,pot));a<-a[order(-a$aupr_corrected,a$test_ligand),];a$activity_rank<-seq_len(nrow(a));a$selected<-a$aupr_corrected>0&a$activity_rank<=10
   for(k in names(base))a[[k]]<-base[[k]]
   activities[[length(activities)+1]]<-a
   for(lig in a$test_ligand[a$selected]){
    weights<-lt[b,lig];ord<-order(-weights,names(weights));ranks<-setNames(seq_along(ord),b[ord])
    for(g in gs){z<-data.frame(ligand=lig,gene=g,weight=as.numeric(lt[g,lig]),target_rank=unname(ranks[g]));for(k in names(base))z[[k]]<-base[[k]];links[[length(links)+1]]<-z}
   }
   base$status<-'success';runs[[length(runs)+1]]<-base
  }
 }
}
wr(do.call(rbind,runs),'targeted_nichenet_inventory.tsv')
if(length(activities))wr(do.call(rbind,activities),'targeted_nichenet_activity.tsv')
if(length(links))wr(do.call(rbind,links),'targeted_nichenet_targets.tsv')
writeLines('success',file.path(O,'logs/targeted_nichenet_complete.txt'))
cat('TARGETED_COMPLETE',length(runs),length(activities),'\n')
