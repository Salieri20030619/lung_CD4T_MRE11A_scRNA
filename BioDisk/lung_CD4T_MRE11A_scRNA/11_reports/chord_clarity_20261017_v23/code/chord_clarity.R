suppressPackageStartupMessages(library(circlize))
P<-'/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports';B<-file.path(P,'client_layout_20261017_v22');O<-file.path(P,'chord_clarity_20261017_v23');dir.create(file.path(O,'figures'),recursive=TRUE,showWarnings=FALSE);dir.create(file.path(O,'plot_data'),showWarnings=FALSE)
rd<-function(p)read.delim(p,stringsAsFactors=FALSE,check.names=FALSE)
chord<-function(edges,cols,labs,order,linkcols,cex=1,gap=13){
 circos.clear();circos.par(start.degree=90,gap.degree=gap,cell.padding=c(0,0,0,0),track.margin=c(.01,.01),points.overflow.warning=FALSE)
 chordDiagram(edges,order=order,grid.col=cols,col=linkcols,transparency=.35,directional=1,direction.type=c('diffHeight','arrows'),diffHeight=mm_h(3),link.arr.type='big.arrow',link.arr.length=.14,link.arr.width=.14,link.border='white',link.lwd=1.6,annotationTrack='grid',preAllocateTracks=list(track.height=.17),link.sort=FALSE,link.decreasing=TRUE)
 circos.trackPlotRegion(track.index=1,bg.border=NA,panel.fun=function(x,y){sec<-get.cell.meta.data('sector.index');circos.text(mean(get.cell.meta.data('xlim')),.3,labs[[sec]],facing='bending.inside',niceFacing=TRUE,cex=cex)})
 circos.clear()
}
a<-rd(file.path(B,'plot_data/F11.tsv'));statecols<-c(Naive='#1696A3',ActCD4='#E28036',Treg='#8E5BA6')
draw_overview<-function(){
 par(mfrow=c(2,2),mar=c(3.1,1.2,3,1.2),oma=c(.8,0,0,0),xpd=NA)
 for(sender in c('AllMac','AMlike'))for(gr in c('Saline','Silica')){
  z<-a[a$sender==sender&a$group==gr,];z<-z[match(c('Naive','ActCD4','Treg'),z$scope),];e<-data.frame(from=sender,to=z$scope,value=z$eligible_interactions);nodes<-c(sender,rev(z$scope));co<-c(setNames('#526675',sender),statecols);lab<-setNames(nodes,nodes)
  chord(e,co,lab,nodes,statecols[z$scope],.9,17);title(main=paste(sender,'-> CD4 |',gr),line=.7,cex.main=.96)
  legend('bottom',legend=paste(z$scope,z$eligible_interactions,sep=': '),fill=statecols[z$scope],bty='n',horiz=TRUE,inset=c(0,-.17),cex=.67,x.intersp=.45,y.intersp=.8)
 }
 invisible(NULL)
}
pdf(file.path(O,'figures/F11.pdf'),width=7,height=7.8);draw_overview();dev.off();png(file.path(O,'figures/F11.png'),width=7,height=7.8,units='in',res=300);draw_overview();dev.off()
r<-rd(file.path(B,'plot_data/F12.tsv'));lig<-unique(r$ligand);rec<-unique(r$receptor_required_subunits);pal<-setNames(c('#0072B2','#E69F00','#009E73','#CC79A7','#D55E00','#56B4E9'),lig);nodes<-c(paste0('L:',lig),paste0('R:',rev(rec)));cols<-c(setNames(pal,paste0('L:',lig)),setNames(rep('#617789',length(rec)),paste0('R:',rec)));labs<-setNames(nodes,nodes);e<-data.frame(from=paste0('L:',r$ligand),to=paste0('R:',r$receptor_required_subunits),value=1)
draw_detail<-function(){
 par(mar=c(8,2,3,2),xpd=NA);chord(e,cols,labs,nodes,pal[r$ligand],.98,8)
 title(main='Macrophage ligands -> Treg receptors',line=.8,cex.main=1.08)
 legend('bottom',legend=paste(r$ligand,'->',r$receptor_required_subunits),fill=pal[r$ligand],ncol=2,bty='n',inset=c(0,-.12),cex=.82,x.intersp=.6,y.intersp=1.1)
 mtext('L = sender ligand; R = receiver receptor; arrowhead points to R',side=1,line=7.0,cex=.8)
}
pdf(file.path(O,'figures/F12.pdf'),width=7,height=8.6);draw_detail();dev.off();png(file.path(O,'figures/F12.png'),width=7,height=8.6,units='in',res=300);draw_detail();dev.off()
file.copy(file.path(B,'plot_data/F11.tsv'),file.path(O,'plot_data/F11.tsv'),overwrite=TRUE);file.copy(file.path(B,'plot_data/F12.tsv'),file.path(O,'plot_data/F12.tsv'),overwrite=TRUE)
writeLines('success',file.path(O,'R_complete.txt'));quit(save='no',status=0,runLast=FALSE)