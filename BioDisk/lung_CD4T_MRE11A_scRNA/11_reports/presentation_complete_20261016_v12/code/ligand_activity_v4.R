args<-commandArgs(TRUE);O<-args[1];P<-dirname(dirname(O));suppressPackageStartupMessages(library(nichenetr))
canon<-function(x){x<-toupper(trimws(x));a<-c(MRE11A='MRE11',ATP5A1='ATP5F1A',ATP5B='ATP5F1B',ATP5C1='ATP5F1C',ATP5G1='ATP5MC1',ATP5O='ATP5PO',PARK2='PRKN');j<-x%in%names(a);x[j]<-a[x[j]];x}
lt<-readRDS(file.path(P,'00_raw/nichenet_v2_mouse/ligand_target_matrix_nsga2r_final_mouse.rds'));rownames(lt)<-canon(rownames(lt));colnames(lt)<-canon(colnames(lt))
elig<-read.delim(file.path(O,'02_sender/ligand_scope_eligibility.tsv'),check.names=FALSE)
elig$expression_supported<-tolower(as.character(elig$expression_supported))=='true'
genes<-read.delim(file.path(P,'02_analysis/mechanism_step1_20260927_152218_s4V5eF/comm_input/CD4_genes.tsv'),check.names=FALSE)$canonical
all<-list();sets<-list();status<-list();counter<-0
for(scope in c('Core','Naive','ActCD4','Treg')){
 de<-read.delim(gzfile(file.path(O,'03_receiver',paste0(scope,'_DE_all_genes.tsv.gz'))),check.names=FALSE)
 bg<-intersect(intersect(de$canonical,rownames(lt)),genes)
 potential<-unique(elig$ligand[elig$scope==scope & elig$expression_supported]);potential<-intersect(potential,colnames(lt))
 for(direction in c('Up','Down')){
  gs<-intersect(de$canonical[de$PValue<0.05 & if(direction=='Up')de$logFC>0 else de$logFC<0],bg)
  sets[[paste(scope,direction)]]<-data.frame(scope=scope,direction=direction,gene=bg,in_target=bg%in%gs)
  evaluable<-length(gs)>=10 && length(bg)>=500 && length(setdiff(bg,gs))>=10 && length(potential)>=2
  status[[paste(scope,direction)]]<-data.frame(scope=scope,direction=direction,background=length(bg),targets=length(gs),potential_ligands=length(potential),evaluable=evaluable)
  if(!evaluable)next
  act<-as.data.frame(nichenetr::predict_ligand_activities(geneset=gs,background_expressed_genes=bg,ligand_target_matrix=lt,potential_ligands=potential))
  act<-act[order(-act$aupr_corrected,act$test_ligand),];act$rank<-seq_len(nrow(act));act$scope<-scope;act$direction<-direction;act$target_rule<-'rawP<0.05 and sign; no new fold-change gate';act$target_n<-length(gs);act$background_n<-length(bg)
  all[[paste(scope,direction)]]<-act
  cat('ACTIVITY_DONE',scope,direction,'targets',length(gs),'ligands',nrow(act),'\n');flush.console()
 }
}
write.table(do.call(rbind,all),file.path(O,'02_sender/ligand_activity_all_states.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(do.call(rbind,sets),file.path(O,'02_sender/receiver_target_background_members.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(do.call(rbind,status),file.path(O,'02_sender/ligand_activity_evaluability.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(c('Predictive ligand-target ranking only; corrected AUPR is not P and not causal probability.','Prior is unsigned. Down-target ranking does not establish an inhibitory ligand effect.','New all-nonempty case-expression coverage; legacy >=10% detection is a robustness flag, not an eligibility gate.',capture.output(sessionInfo())),file.path(O,'02_sender/methods_and_session.txt'))
cat('LIGAND_ACTIVITY_ALL_STATES_PASS\n')