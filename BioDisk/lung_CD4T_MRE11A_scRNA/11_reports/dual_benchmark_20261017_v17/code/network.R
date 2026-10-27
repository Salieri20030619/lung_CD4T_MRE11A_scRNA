suppressPackageStartupMessages(library(jsonlite))
p<-'/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA';b<-file.path(p,'11_reports/three_cohort_four_axis_20261017_v16');o<-file.path(p,'11_reports/dual_benchmark_20261017_v17')
genes<-fromJSON(file.path(o,'genes.json'));a<-read.delim(file.path(b,'tables/nominal_DE_axis_membership.tsv'),stringsAsFactors=FALSE);a<-a[a$gene%in%genes & a$contrast!='time_omnibus',];acts<-read.delim(file.path(b,'tables/NicheNet_activity_all.tsv'),stringsAsFactors=FALSE);bg<-read.delim(file.path(b,'inputs/receiver_background.tsv'),stringsAsFactors=FALSE)
lt<-readRDS(file.path(p,'00_raw/nichenet_v2_mouse/ligand_target_matrix_nsga2r_final_mouse.rds'));canon<-function(x){x<-toupper(x);al<-c(MRE11A='MRE11',ATP5A1='ATP5F1A',ATP5B='ATP5F1B',ATP5C1='ATP5F1C',ATP5G1='ATP5MC1',ATP5O='ATP5PO',PARK2='PRKN');i<-x%in%names(al);x[i]<-al[x[i]];x};rownames(lt)<-canon(rownames(lt));colnames(lt)<-canon(colnames(lt));out<-list()
for(fam in unique(a$family))for(dir in c('up','down')){
 z<-a[a$family==fam & a$direction==dir,];if(!nrow(z))next
 ac<-acts[acts$family==fam & acts$direction==dir & acts$aupr_corrected>0,];if(!nrow(ac))next
 ds<-ac$dataset[1];sc<-ac$scope[1];back<-intersect(bg$gene[bg$dataset==ds & bg$scope==sc],rownames(lt));gg<-intersect(z$gene,back)
 for(i in seq_len(nrow(ac))){lig<-ac$test_ligand[i];w<-lt[back,lig];r<-rank(-w,ties.method='min');names(r)<-back;ii<-match(gg,back);yes<-w[ii]>0;if(!any(yes))next;out[[length(out)+1]]<-data.frame(gene=gg[yes],family=fam,direction=dir,scope=sc,dataset=ds,ligand=lig,weight=as.numeric(w[ii][yes]),target_rank=as.numeric(r[ii][yes]),activity_rank=ac$activity_rank[i],aupr_corrected=ac$aupr_corrected[i],v16_reference_rank=ac$activity_rank[i]<=30 & r[ii][yes]<=250)}
}
z<-do.call(rbind,out);write.table(z,gzfile(file.path(o,'tables/network_links_no_rank_exclusion.tsv.gz')),sep='\t',quote=FALSE,row.names=FALSE,na='NA');cat('LINKS',nrow(z),'GENES',length(unique(z$gene)),'\n')