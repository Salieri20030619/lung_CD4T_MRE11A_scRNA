O<-'/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/death_priority_20261018_v2'
p<-readRDS(file.path(O,'inputs/pathways_symbol.rds'));rows<-lapply(names(p),function(db)do.call(rbind,lapply(names(p[[db]]),function(k)data.frame(database=db,pathway=k,gene=p[[db]][[k]]))))
con<-gzfile(file.path(O,'tables/pathway_symbol_members.tsv.gz'),'wt',encoding='UTF-8');write.table(do.call(rbind,rows),con,sep='\t',row.names=FALSE,quote=FALSE);close(con)
