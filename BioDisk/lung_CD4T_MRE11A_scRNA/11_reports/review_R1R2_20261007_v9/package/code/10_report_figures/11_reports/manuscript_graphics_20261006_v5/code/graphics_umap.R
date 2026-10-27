P<-'/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA';O<-file.path(P,'11_reports/manuscript_graphics_20261006_v5')
suppressPackageStartupMessages({library(ggplot2);library(cowplot);library(jsonlite)})
rd<-function(f)read.delim(f,check.names=FALSE,stringsAsFactors=FALSE)
b<-readRDS(file.path(O,'tables/Human_lung_embedding_frozen.rds'))$metadata
c<-readRDS(file.path(O,'tables/Human_CD4_embedding_frozen.rds'))$metadata
r<-rd(file.path(P,'02_analysis/mechanism_step2_20260927_154919_EsN4gk/input/human_cells.tsv'))
stopifnot(nrow(b)==192575,nrow(c)==6523,length(unique(b$subject_id))==20,length(unique(c$subject_id))==19,setequal(r$cell,c$cell))
m<-match(c$cell,r$cell);stopifnot(all(c$cd4_label==r$subtype[m]),all(c$subject_id==r$sample[m]),all(c$group==r$group[m]))
smap<-c('reused_GSE135893'='G135','reused_GSE227136'='G227','GSE288887_new'='G288');c$source<-unname(smap[c$source_type]);b$source<-unname(smap[b$source_type]);stopifnot(all(c$source==r$series[m]))
expect<-rd(file.path(P,'10_tables/broad_clustering_v1/GSE288887_preliminary_lineage_summary.tsv'));ct<-table(b$lineage_prelim);stopifnot(all(as.integer(ct[expect$lineage_prelim])==expect$cells))
stopifnot(all(is.finite(c$UMAP_1)),all(is.finite(b$UMAP_1)),!anyDuplicated(b$cell),!anyDuplicated(c$cell),all(c$cell%in%b$cell))
immune<-b[b$lineage_prelim%in%c('T','NK','B','My','Mast'),];stopifnot(nrow(immune)==57706)
lc<-c('T'='#B83238','NK'='#68459E','B'='#2166AC','My'='#A65D00','Mast'='#C03B82','Epi'='#007F72','Endo'='#364D80','Str'='#566B2F')
cc<-c('CD4-CCR7'='#2166AC','CD4-GZMA'='#B83238','CD4-IL7R'='#007F72','Treg'='#70429B','CD4state'='#956018')
sc<-c('G135'='#275A87','G227'='#A54B24','G288'='#006C53')
set.seed(20261006)
um<-function(d,col,palette,title,subtitle,small=.1,labels=TRUE){
 d<-d[sample.int(nrow(d)),];d$display<-d[[col]]
 pp<-ggplot(d,aes(UMAP_1,UMAP_2,color=display))+geom_point(size=small,alpha=.75,shape=16)+scale_color_manual(values=palette,drop=FALSE)+coord_equal(clip='off')+
 labs(title=title,subtitle=subtitle,color=NULL,x='UMAP 1',y='UMAP 2')+theme_classic(base_size=13)+theme(plot.title=element_text(face='bold',size=14),plot.subtitle=element_text(size=11),legend.position='bottom',legend.text=element_text(size=11),legend.key.width=unit(.45,'cm'),legend.margin=margin(0,0,0,0),axis.text=element_blank(),axis.ticks=element_blank(),plot.margin=margin(8,8,8,8))+
 guides(color=guide_legend(nrow=if(length(palette)>5)2 else 1,override.aes=list(size=3,alpha=1)))
 if(labels){
  med<-aggregate(cbind(UMAP_1,UMAP_2)~display,d,median)
  med$label<-sub('CD4-','',med$display,fixed=TRUE)
  pp<-pp+geom_label(data=med,aes(label=label),fill='white',label.padding=unit(.1,'lines'),size=4.5,fontface='bold',show.legend=FALSE,linewidth=0)
 }
 pp+theme(plot.subtitle=element_text(size=13),legend.text=element_text(size=13))
}
pa<-um(b,'lineage_prelim',lc,'A  Human lung atlas','192,575 cells | 20 donors | frozen whole-lung UMAP',.08)
pb<-um(immune,'lineage_prelim',lc[c('T','NK','B','My','Mast')],'B  Major immune lineages','57,706 cells | same embedding as A',.12)+xlim(range(b$UMAP_1))+ylim(range(b$UMAP_2))
pc<-um(c,'cd4_label',cc,'C  Human CD4 subgroups','6,523 cells | 19 donors | frozen CD4 UMAP',.38)
pd<-um(c,'source',sc,'D  CD4 cohort source','Same embedding as C | G135 / G227 / G288',.38,FALSE)
combo<-plot_grid(pa,pb,pc,pd,ncol=2,align='hv')
savefig<-function(plot,stem,w,h){ggsave(file.path(O,'figures',paste0(stem,'.pdf')),plot,width=w,height=h,device=cairo_pdf);ggsave(file.path(O,'figures',paste0(stem,'.png')),plot,width=w,height=h,dpi=300,bg='white')}
savefig(combo,'Human_UMAP_atlas_CD4',12,10)
savefig(pa,'Human_lung_UMAP',8,7);savefig(pc,'Human_CD4_UMAP',8,7)
qa<-list(species='Homo sapiens',dataset='GSE288887 integrated collection, including reused GSE135893 and GSE227136',whole_lung_cells=nrow(b),whole_lung_donors=length(unique(b$subject_id)),whole_lung_samples=length(unique(b$sample_id)),immune_cells=nrow(immune),CD4_cells=nrow(c),CD4_donors=length(unique(c$subject_id)),CD4_samples=length(unique(c$sample_id)),whole_lung_source='04_preprocessing/GSE288887_broad_clustered_scored.rds',CD4_source='06_CD4T_analysis/CD4_final/GSE288887_CD4_final.rds',lineage_annotation='lineage_prelim inherited from original broad marker annotation',CD4_annotation='cd4_label, exact match to frozen disease-analysis roster',coordinates='Original umap@cell.embeddings, no recomputation or coordinate transformation; panels A/B share one embedding and C/D share the separate CD4 embedding',CD4_cell_set_match=TRUE,CD4_donor_group_source_subtype_match=TRUE,lineage_counts_match=TRUE,lineage_counts=as.list(ct),CD4_subgroup_counts=as.list(table(c$cd4_label)),new_UMAP_runs=0,new_clustering_runs=0)
write_json(qa,file.path(O,'provenance/UMAP_QA.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(O,'provenance/UMAP_sessionInfo.txt'))
writeLines('complete',file.path(O,'provenance/umap.done'));cat('UMAP_PLOTTED\n')
