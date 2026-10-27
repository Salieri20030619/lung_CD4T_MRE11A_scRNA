P <- '/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA';O<-file.path(P,'11_reports/manuscript_graphics_20261006_v5')
suppressPackageStartupMessages({library(enrichplot);library(DOSE);library(ggplot2);library(cowplot);library(jsonlite)})
cols<-c('#0072B2','#B33A17','#007D5C','#8B3B88','#8F6500')
spec<-list(selection='Top five finite NES among all originally tested pathways, ordered by descending absolute NES then pathway ID; no P/FDR filtering or pathway exclusion',human_scopes=c('Core','Primary'),databases=c('KEGG','GO_BP'),input='Original saved GSEA objects; source-adjusted existing donor pseudobulk rank sign(logFC)*sqrt(QLF)',reanalysis=FALSE,pvalue_table=FALSE,ES_geom='line',base_size=11,relative_heights=c(1.5,.4,1),subplots=1:3,width_inches=10,height_inches=8,png_dpi=300)
write_json(spec,file.path(O,'provenance/GSEA_display_spec.json'),pretty=TRUE,auto_unbox=TRUE)
manifest<-list();checks<-list();allsel<-list()
for(nm in c('Core_KEGG','Core_GO','Primary_KEGG','Primary_GO')){
 g<-readRDS(file.path(O,'tables',paste0('Human_',nm,'_frozen_gsea.rds')));r<-g@result
 ix<-order(-abs(r$NES),r$ID,na.last=TRUE)[1:5];z<-r[ix,];code<-paste0(if(grepl('KEGG',nm))'K'else'G',1:5)
 short<-if(grepl('KEGG',nm))c('Ribosome','OXPHOS',"Parkinson disease",'Cardiac muscle contraction',"Huntington disease")[match(z$ID,c('KEGG_RIBOSOME','KEGG_OXIDATIVE_PHOSPHORYLATION','KEGG_PARKINSONS_DISEASE','KEGG_CARDIAC_MUSCLE_CONTRACTION','KEGG_HUNTINGTONS_DISEASE'))]else rep(NA_character_,5)
 go_short<-c('GOBP_NEGATIVE_REGULATION_OF_SYNCYTIUM_FORMATION_BY_PLASMA_MEMBRANE_FUSION'='Syncytium formation (-reg.)','GOBP_REGULATION_OF_SYNCYTIUM_FORMATION_BY_PLASMA_MEMBRANE_FUSION'='Syncytium formation (reg.)','GOBP_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT'='ATP-coupled electron transport','GOBP_OXIDATIVE_PHOSPHORYLATION'='OXPHOS','GOBP_PROTON_MOTIVE_FORCE_DRIVEN_ATP_SYNTHESIS'='Proton-driven ATP synthesis')
 go_short<-c(go_short,'GOBP_CYTOPLASMIC_TRANSLATION'='Cytoplasmic translation')
 if(!grepl('KEGG',nm))short<-unname(go_short[z$ID])
 short[is.na(short)]<-z$Description[is.na(short)]
 short<-gsub('_',' ',short,fixed=TRUE)
 labs<-paste(code,short)
 for(i in seq_len(nrow(z))){
  gs<-enrichplot:::gsInfo(g,z$ID[i]);es<-if(max(gs$runningScore)>abs(min(gs$runningScore)))max(gs$runningScore)else min(gs$runningScore)
  delta<-abs(es-z$enrichmentScore[i]);stopifnot(delta<1e-8)
  hits<-sum(gs$position==1);stopifnot(hits==z$setSize[i])
  checks[[length(checks)+1]]<-data.frame(specification=nm,ID=z$ID[i],ES_saved=z$enrichmentScore[i],ES_curve=es,absolute_delta=delta,hits,size=z$setSize[i])
 }
 plotg<-g;plotg@result$Description[match(z$ID,plotg@result$ID)]<-labs
 title<-paste('Human',sub('_',' | ',nm,fixed=TRUE),'| Top 5 by |NES|')
 pp<-gseaplot2(plotg,geneSetID=z$ID,title=title,color=cols,base_size=11,rel_heights=c(1.5,.4,1),subplots=1:3,pvalue_table=FALSE,ES_geom='line')
 # Preserve the requested three-layer plot; separate statistics footer uses original values.
 ft<-paste(sprintf('%s   NES %+.3f   P %.3g   BH %.3g',code,z$NES,z$pvalue,z$p.adjust),collapse='    |    ')
 footer<-ggdraw()+draw_label(paste(sprintf('%s   NES %+.3f   P %.3g   BH %.3g',code,z$NES,z$pvalue,z$p.adjust),collapse='\n'),x=.015,y=.5,hjust=0,size=10.5,fontfamily='sans')+
  draw_label(sprintf('%s tested pathways\n%s ranked genes\nPositive rank: CWP\nBH: full original family',nrow(r),length(g@geneList)),x=.985,y=.5,hjust=1,size=10.5,fontfamily='sans')
 panels<-if(inherits(pp,'gglist'))unclass(pp)else list(pp)
 panels[[1]]<-panels[[1]]+theme(legend.text=element_text(size=12),axis.text=element_text(size=11))
 panels[[3]]<-panels[[3]]+theme(axis.text=element_text(size=11))
 coreplot<-plot_grid(plotlist=panels,ncol=1,rel_heights=c(1.5,.4,1),align='v',axis='lr')
 combo<-plot_grid(coreplot,footer,ncol=1,rel_heights=c(1,.17))
 stem<-paste0('Human_',nm,'_Top5_GSEA')
 ggsave(file.path(O,'figures',paste0(stem,'.pdf')),combo,width=10,height=8,device=cairo_pdf)
 ggsave(file.path(O,'figures',paste0(stem,'.png')),combo,width=10,height=8,dpi=300,bg='white')
 z$display_code<-code;z$display_label<-short;z$specification<-nm;z$family_size<-nrow(r);z$ranked_genes<-length(g@geneList);allsel[[nm]]<-z
 manifest[[nm]]<-list(stem=stem,title=title,tested=nrow(r),ranked_genes=length(g@geneList),selected=z[,c('display_code','ID','Description','display_label','NES','pvalue','p.adjust','setSize')])
 cat('PLOTTED',nm,'\n');flush.console()
}
write.table(do.call(rbind,allsel),file.path(O,'tables/GSEA_Top5_original_statistics.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(do.call(rbind,checks),file.path(O,'tables/GSEA_running_ES_QA.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write_json(manifest,file.path(O,'provenance/GSEA_figure_manifest.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
writeLines(capture.output(sessionInfo()),file.path(O,'provenance/GSEA_plot_sessionInfo.txt'))
writeLines('complete',file.path(O,'provenance/gsea.done'))
