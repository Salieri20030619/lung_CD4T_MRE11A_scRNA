P<-'/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA';O<-file.path(P,'11_reports/manuscript_graphics_20261006_v5')
suppressPackageStartupMessages({library(DOSE);library(jsonlite)})
rd<-function(f)read.delim(f,check.names=FALSE,stringsAsFactors=FALSE)
ref<-rd(file.path(P,'10_tables/CD4_final/pathway_v1/GO_KEGG_GSEA.tsv'));ans<-list()
for(nm in c('Core_KEGG','Core_GO','Primary_KEGG','Primary_GO')){
 g<-readRDS(file.path(O,'tables',paste0('Human_',nm,'_frozen_gsea.rds')));r<-g@result;tag<-strsplit(nm,'_')[[1]];db<-if(tag[2]=='GO')'GO_BP'else'KEGG'
 s<-ref[ref$analysis==tag[1]&ref$database==db,];stopifnot(nrow(r)==nrow(s),setequal(r$ID,s$ID));s<-s[match(r$ID,s$ID),]
 delta<-sapply(c('NES','enrichmentScore','pvalue','p.adjust'),function(k){a<-r[[k]];b<-s[[k]];stopifnot(identical(is.na(a),is.na(b)));max(abs(a-b)/pmax(abs(a),abs(b),1e-300),na.rm=TRUE)})
 stopifnot(all(delta<1e-10));ans[[nm]]<-list(rows=nrow(r),maximum_relative_difference=as.list(delta))
}
es<-rd(file.path(O,'tables/GSEA_running_ES_QA.tsv'));stopifnot(nrow(es)==20,max(es$absolute_delta)<1e-8,all(es$hits==es$size))
write_json(list(original_result_comparisons=ans,curve_checks=nrow(es),max_ES_absolute_delta=max(es$absolute_delta),all_hit_counts_match=TRUE),file.path(O,'provenance/GSEA_final_numeric_QA.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
cat('GSEA_FULL_RESULT_AND_20_CURVE_QA_PASSED\n')
