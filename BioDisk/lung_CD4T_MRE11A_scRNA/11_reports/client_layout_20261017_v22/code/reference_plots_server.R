suppressPackageStartupMessages({library(DOSE);library(enrichplot);library(ggplot2);library(cowplot);library(circlize)})
O<-'/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/client_layout_20261017_v22';B<-'/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/client_rebuild_20261017_v20';W<-'/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/feedback_revision_20261017_v19';dir.create(file.path(O,'figures'),showWarnings=FALSE);dir.create(file.path(O,'plot_data'),showWarnings=FALSE)
rd<-function(p)read.delim(p,check.names=FALSE,stringsAsFactors=FALSE)
wr<-function(x,n)write.table(as.data.frame(x),file.path(O,'plot_data',n),sep='\t',quote=FALSE,row.names=FALSE,na='NA',fileEncoding='UTF-8')
all<-rd(file.path(B,'tables/GSEA_all.tsv'));mem<-rd(file.path(B,'tables/enrichment_members.tsv'));qa<-list();selected<-list()
for(sc in c('Treg','ActCD4')){
 z<-all[all$scope==sc&all$database=='KEGG',];r<-rd(file.path(B,paste0('tables/GSEA_ranks_',sc,'.tsv')));stats<-setNames(r$rank_stat,r$gene);sets<-split(mem$gene[mem$database=='KEGG'],mem$pathway[mem$database=='KEGG']);sets<-lapply(sets,function(x)intersect(unique(x),names(stats)))
 tab<-data.frame(ID=z$pathway,Description=sub(' - Mus musculus (house mouse)','',z$description,fixed=TRUE),setSize=z$size,enrichmentScore=z$ES,NES=z$NES,pvalue=z$pval,p.adjust=z$padj,qvalue=z$padj,rank=0,leading_edge='',core_enrichment=gsub(';','/',z$leadingEdge,fixed=TRUE),stringsAsFactors=FALSE);rownames(tab)<-tab$ID
 obj<-methods::new('gseaResult',result=tab,organism='Mus musculus',setType='KEGG',geneSets=sets,geneList=stats,keytype='SYMBOL',params=list(exponent=1),readable=FALSE)
 eligible<-tab[is.finite(tab$pvalue)&tab$pvalue<.05,];eligible<-eligible[order(-abs(eligible$NES),eligible$ID),];top<-head(eligible,5);ids<-top$ID;id<-if(sc=='Treg')'F09'else'F10'
 # Same three-layer gseaplot2 geometry and proportions as user reference.
 cols<-c('#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F');pal<-setNames(cols,top$Description)
 panels<-lapply(1:3,function(k)enrichplot::gseaplot2(obj,geneSetID=ids,title=if(k==1)paste(sc,'| KEGG GSEA')else'',base_size=11,rel_heights=c(1.5,.4,1),subplots=k,pvalue_table=FALSE,ES_geom='line',color=pal))
 legend<-cowplot::get_legend(panels[[1]]+theme(legend.position='top',legend.justification='left',legend.text=element_text(size=10),legend.key.height=grid::unit(.42,'cm'))+guides(color=guide_legend(ncol=1))+scale_color_manual(values=pal,labels=function(x)vapply(x,function(s)paste(strwrap(s,width=62),collapse='\n'),character(1))))
 panels[[1]]<-panels[[1]]+theme(legend.position='none',plot.title=element_text(size=13,face='bold'))
 panels[[2]]<-panels[[2]]+theme(legend.position='none');panels[[3]]<-panels[[3]]+theme(legend.position='none')
 curves<-cowplot::plot_grid(plotlist=panels,ncol=1,rel_heights=c(1.5,.4,1),align='v',axis='lr')
 pv<-function(x)ifelse(x<=1e-8,'<=1e-8',formatC(x,format='g',digits=3));labels<-data.frame(pathway=top$Description,NES=formatC(top$NES,format='f',digits=2),P=pv(top$pvalue),FDR=formatC(top$p.adjust,format='g',digits=3));names(labels)<-c('Pathway','NES','P','FDR')
 labplot<-ggplot()+theme_void()+xlim(0,1)+ylim(0,6)
 for(j in 1:4){xx<-c(.01,.68,.80,.94)[j];labplot<-labplot+annotate('text',x=xx,y=5.7,label=names(labels)[j],hjust=if(j==1)0 else .5,size=3.1,fontface='bold');for(i in 1:5)labplot<-labplot+annotate('text',x=xx,y=5.6-i,label=if(j==1)paste0(i,'. ',labels[i,j])else labels[i,j],hjust=if(j==1)0 else .5,size=2.8)}
 final<-cowplot::plot_grid(legend,curves,labplot,ncol=1,rel_heights=c(.95,3.7,1.15))
 ggsave(file.path(O,'figures',paste0(id,'.pdf')),final,width=7,height=9.0,limitsize=FALSE,bg='white');ggsave(file.path(O,'figures',paste0(id,'.png')),final,width=7,height=9.0,dpi=300,limitsize=FALSE,bg='white')
 saveRDS(obj,file.path(O,paste0(sc,'_gseaResult_from_saved_statistics.rds')));top$scope<-sc;top$database<-'KEGG';top$figure<-id;wr(top,paste0(id,'.tsv'));selected[[length(selected)+1]]<-top
 for(pid in ids){hit<-names(stats)%in%sets[[pid]];inc<-ifelse(hit,abs(stats)/sum(abs(stats[hit])),-1/sum(!hit));es<-cumsum(inc);mx<-max(es);mn<-min(es);val<-if(mx> -mn)mx else mn;stopifnot(abs(val-tab[pid,'enrichmentScore'])<1e-8);qa[[length(qa)+1]]<-data.frame(scope=sc,pathway=pid,ES_original=tab[pid,'enrichmentScore'],ES_reconstructed=val,error=abs(val-tab[pid,'enrichmentScore']))}
}
wr(do.call(rbind,qa),'GSEA_curve_QA.tsv');wr(do.call(rbind,selected),'GSEA_top5_selected.tsv')
# Chord overview: separate overlapping sender definitions, no sum across AllMac/AMlike.
a<-rd(file.path(B,'plot_data/F11.tsv'));a<-a[a$scope%in%c('Naive','ActCD4','Treg'),];wr(a,'F11.tsv')
chord<-function(edges,cols,labels,cex=.85,gap=8){
 circos.clear();circos.par(start.degree=90,gap.degree=gap,cell.padding=c(0,0,0,0),track.margin=c(.01,.01),points.overflow.warning=FALSE)
 chordDiagram(edges,grid.col=cols,transparency=.3,directional=1,direction.type=c('diffHeight','arrows'),link.arr.type='big.arrow',annotationTrack='grid',preAllocateTracks=list(track.height=.16),link.sort=TRUE,link.decreasing=TRUE)
 circos.trackPlotRegion(track.index=1,bg.border=NA,panel.fun=function(x,y){sec<-get.cell.meta.data('sector.index');xlim<-get.cell.meta.data('xlim');circos.text(mean(xlim),.3,labels[[sec]],facing='bending.inside',niceFacing=TRUE,cex=cex)})
 circos.clear()
}
draw_overview<-function(){
 par(mfrow=c(2,2),mar=c(2.5,1.3,3,1.3),oma=c(0,0,1,0),xpd=NA)
 for(sender in c('AllMac','AMlike'))for(gr in c('Saline','Silica')){
  z<-a[a$sender==sender&a$group==gr,];edges<-data.frame(from=sender,to=z$scope,value=z$eligible_interactions);nodes<-c(sender,'Naive','ActCD4','Treg');cols<-setNames(c('#3278a8','#4c9e9c','#dc8451','#8b74aa'),nodes);labels<-setNames(nodes,nodes);chord(edges,cols,labels,.9,9);title(main=paste(sender,'|',gr),cex.main=1.05,line=1);mtext(paste(z$scope,z$eligible_interactions,sep='=',collapse='  '),side=1,line=.4,cex=.78)
 }
}
pdf(file.path(O,'figures/F11.pdf'),width=7,height=7.4);draw_overview();dev.off();png(file.path(O,'figures/F11.png'),width=7,height=7.4,units='in',res=300);draw_overview();dev.off()
# Detailed receptor chord: Treg routes from the already assessed 18-route set.
r<-rd(file.path(W,'tables/R02_all18_routes_evidence.tsv'));r<-r[r$receiver_scope=='Treg'&r$ligand!='PROS1',];r<-unique(r[,c('ligand','receptor_required_subunits','interaction')]);wr(r,'F12.tsv');lig<-unique(r$ligand);rec<-unique(r$receptor_required_subunits);nodes<-c(paste0('L:',lig),paste0('R:',rec));cols<-setNames(c(rep('#d78942',length(lig)),rep('#468da6',length(rec))),nodes);labs<-setNames(c(lig,rec),nodes);edges<-data.frame(from=paste0('L:',r$ligand),to=paste0('R:',r$receptor_required_subunits),value=1)
draw_detail<-function(){par(mar=c(3,2,3,2),xpd=NA);chord(edges,cols,labs,1.05,5);title(main='Macrophage ligand - Treg receptor candidates',cex.main=1.05,line=.8);legend('bottom',legend=c('Sender ligand','Treg receptor'),fill=c('#d78942','#468da6'),horiz=TRUE,bty='n',inset=c(0,-.09),cex=.95)}
pdf(file.path(O,'figures/F12.pdf'),width=7,height=7);draw_detail();dev.off();png(file.path(O,'figures/F12.png'),width=7,height=7,units='in',res=300);draw_detail();dev.off()
writeLines(capture.output(sessionInfo()),file.path(O,'R_plot_session.txt'));writeLines('success',file.path(O,'R_plots_complete.txt'))
quit(save='no',status=0,runLast=FALSE)
