suppressPackageStartupMessages({library(enrichplot);library(DOSE)})
O<-'/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/manuscript_graphics_20261006_v5'
g<-readRDS(file.path(O,'tables/Human_Core_KEGG_frozen_gsea.rds'));ix<-order(-abs(g@result$NES))[1:5]
p<-gseaplot2(g,geneSetID=g@result$ID[ix],subplots=1:3,rel_heights=c(1.5,.4,1),base_size=11,pvalue_table=FALSE,ES_geom='line')
print(class(p));str(p,max.level=1);cat('P_NAMES',names(p),'\n');cat('ATTRS',names(attributes(p)),'\n');str(attributes(p),max.level=2)
