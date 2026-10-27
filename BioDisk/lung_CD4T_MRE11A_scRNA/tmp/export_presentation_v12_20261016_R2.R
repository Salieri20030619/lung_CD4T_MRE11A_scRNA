args <- commandArgs(TRUE)
stopifnot(length(args)==1)
O <- args[1]
if(dir.exists(O)) stop('Output already exists; preserve previous export')
dir.create(O, recursive=TRUE)
P <- '/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA'
M <- file.path(P,'02_analysis/mechanism_step1_20260927_152218_s4V5eF')
V <- file.path(P,'11_reports/customer_sequence_20261015_v2')
for(f in c('tables/module_definition.tsv','tables/module_coverage.tsv','tables/score_controls.tsv.gz','methods.txt','analyze.py')) {
 stopifnot(file.copy(file.path(M,f),file.path(O,basename(f)),overwrite=FALSE))
}
stopifnot(file.copy(file.path(V,'code/compute_v2.py'),file.path(O,'compute_v2.py'),overwrite=FALSE))
canon <- function(x) { x<-toupper(trimws(x));a<-c(MRE11A='MRE11',ATP5A1='ATP5F1A',ATP5B='ATP5F1B',ATP5C1='ATP5F1C',ATP5G1='ATP5MC1',ATP5O='ATP5PO',PARK2='PRKN');j<-x%in%names(a);x[j]<-a[x[j]];x }
f <- file.path(P,'00_raw/nichenet_v2_mouse/ligand_target_matrix_nsga2r_final_mouse.rds')
lt<-readRDS(f);rownames(lt)<-canon(rownames(lt));colnames(lt)<-canon(colnames(lt))
stopifnot(!anyDuplicated(rownames(lt)),!anyDuplicated(colnames(lt)))
defs<-read.delim(file.path(O,'module_definition.tsv'),stringsAsFactors=FALSE)
targets<-unique(c(defs$gene[defs$module%in%c('PyroExec','ApoExec','Inflam')],c('IL1R2','PENK','TNFRSF9','LAG3','TTC39C','FOXP1','IKZF2','STAT1','AHNAK','MDFIC','ST6GALNAC3')))
out<-expand.grid(target=targets,ligand=colnames(lt),stringsAsFactors=FALSE)
out$in_prior<-out$target%in%rownames(lt)
out$weight<-NA_real_;j<-which(out$in_prior)
out$weight[j]<-lt[cbind(match(out$target[j],rownames(lt)),match(out$ligand[j],colnames(lt)))]
con<-gzfile(file.path(O,'nichenet_death_candidate_prior_weights.tsv.gz'),'wt')
write.table(out,con,sep='\t',quote=FALSE,row.names=FALSE,na='NA')
close(con)
writeLines(c(paste('matrix_source',f),paste('matrix_dimensions',paste(dim(lt),collapse=' x ')),paste('export_rows',nrow(out)),paste('targets_in_prior',sum(targets%in%rownames(lt))),paste('targets_requested',length(targets)),'No predictive model refit; unchanged prior weights extracted.',capture.output(sessionInfo())),file.path(O,'export_provenance.txt'))
writeLines('EXPORT_COMPLETE',file.path(O,'EXPORT_COMPLETE.txt'))
cat('EXPORT_COMPLETE rows=',nrow(out),' targets=',length(targets),'\n')
