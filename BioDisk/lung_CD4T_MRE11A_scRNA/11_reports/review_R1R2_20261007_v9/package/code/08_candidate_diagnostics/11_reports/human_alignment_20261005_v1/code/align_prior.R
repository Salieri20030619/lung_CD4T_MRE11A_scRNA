P<-'/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA';O<-file.path(P,'11_reports/human_alignment_20261005_v1')
lt<-readRDS('/mnt/BioDisk/wangmingyuan/IL17_Candida_integrated_v1/07_synthesis/supplementary_bioinformatics_v1/a6_consensus_regulatory_chain_v1/01_nichenet_resources_v1/ligand_target_matrix.rds')
lig<-strsplit('IL15 TNF CLSTN1 ADAM15 PLAU DSC2 CXCL2 ICAM2 VTN VEGFB GSTP1 RTN4',' ')[[1]]
out<-do.call(rbind,lapply(lig,function(l)data.frame(ligand=l,target=rownames(lt),weight=if(l%in%colnames(lt))lt[,l] else NA)))
con<-gzfile(file.path(O,'results/fixed_ligand_complete_prior_targets.tsv.gz'),'wt');write.table(out,con,sep='\t',quote=FALSE,row.names=FALSE);close(con)
write.table(out[out$target%in%c('IKZF2','STAT1'),],file.path(O,'results/ligand_to_TF_prior_weights.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
cat('PRIOR_COMPLETE\n')
