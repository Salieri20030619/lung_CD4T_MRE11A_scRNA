suppressPackageStartupMessages({library(org.Mm.eg.db);library(GO.db);library(AnnotationDbi)})
o <- '/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/three_cohort_four_axis_20261017_v16'
w <- function(x,n) write.table(as.data.frame(x),file=file.path(o,'resources',n),sep='\t',quote=FALSE,row.names=FALSE,na='NA',fileEncoding='UTF-8')
t <- AnnotationDbi::select(GO.db,keys=keys(GO.db),columns=c('TERM','ONTOLOGY'),keytype='GOID')
t <- t[grepl('intrinsic apoptotic|pyropt|mitochondrial.*apopt|apopt.*mitochondri',t$TERM,ignore.case=TRUE) & t$ONTOLOGY=='BP',]
w(t,'GO_selected_terms.tsv')
a <- AnnotationDbi::select(org.Mm.eg.db,keys=intersect(t$GOID,keys(org.Mm.eg.db,keytype='GOALL')),keytype='GOALL',columns=c('SYMBOL','ENTREZID','EVIDENCEALL'))
a <- merge(a,t,by.x='GOALL',by.y='GOID',all.x=TRUE);a$canonical <- toupper(a$SYMBOL);w(a,'GO_mouse_membership.tsv')
w(AnnotationDbi::metadata(org.Mm.eg.db),'orgMm_metadata.tsv')
writeLines(capture.output(sessionInfo()),file.path(o,'resources','R_session_resources.txt'))
load(file.path(o,'resources','CellChatDB.mouse.rda'))
w(CellChatDB.mouse$interaction,'CellChat_interaction.tsv')
z <- CellChatDB.mouse$complex;z$complex <- rownames(z);w(z,'CellChat_complex.tsv')
cat('GO terms',nrow(t),'genes',length(unique(a$canonical)),'CellChat',nrow(CellChatDB.mouse$interaction),'\n');print(t)