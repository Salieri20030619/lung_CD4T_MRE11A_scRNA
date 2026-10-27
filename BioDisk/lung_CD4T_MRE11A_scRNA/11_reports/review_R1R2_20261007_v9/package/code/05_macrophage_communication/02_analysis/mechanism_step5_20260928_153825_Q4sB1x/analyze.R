args <- commandArgs(TRUE);p1 <- args[1];p3 <- args[2];out <- args[3];cache <- args[4]
suppressPackageStartupMessages({library(Matrix);library(nichenetr)})
set.seed(314943)
rd <- function(p)read.delim(p,check.names=FALSE,stringsAsFactors=FALSE)
wr <- function(x,n)write.table(x,file.path(out,'tables',n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
canon <- function(x){x <- toupper(x);x[x=='MRE11A'] <- 'MRE11';x}
targets <- c('IL1R2','PENK','TTC39C','FOXP1','TNFRSF9','LAG3','MRE11')
g <- rd(file.path(p1,'comm_input/CD4_genes.tsv'));m <- rd(file.path(p1,'comm_input/CD4_cells.tsv'))
g3 <- rd(file.path(p3,'tables/CD4_gene_mapping.tsv'))
stopifnot(identical(g$canonical,g3$canonical),!anyNA(g3$mapped),!anyDuplicated(m$cell),nrow(m)==2311)
con <- gzfile(file.path(p1,'comm_input/CD4_counts.mtx.gz'),'rt');x <- as(readMM(con),'CsparseMatrix');close(con)
stopifnot(nrow(x)==nrow(g),ncol(x)==nrow(m),all(x@x>=0),all(is.finite(x@x)),all(abs(x@x-round(x@x))<1e-6),all(m$library>0))
bgp <- strsplit('SFTPA1 SFTPA2 SFTPB SFTPC SFTPD SCGB1A1 SCGB3A1 SCGB3A2 ABCA3 SLC34A2 NAPSA LYZ LYZ1 LYZ2 S100A8 S100A9 CSF1R LST1 TYROBP FCER1G',' ')[[1]]
bgi <- g$canonical%in%bgp|grepl('^(IGH[VDJMAEG]|IGK[VJC]|IGL[VJC])',g$canonical)
m$bg <- as.numeric(Matrix::colSums(x[bgi,,drop=FALSE]))/m$library
u <- sort(unique(g3$mapped));pr <- sparseMatrix(i=match(g3$mapped,u),j=seq_len(nrow(g3)),x=1,dims=c(length(u),nrow(g3)))
x <- pr%*%x;rownames(x) <- u;colnames(x) <- m$cell
prior <- file.path(cache,'ligand_target_matrix_nsga2r_final_mouse.rds');stopifnot(file.exists(prior))
lt <- readRDS(prior);stopifnot(!anyDuplicated(rownames(lt)),!anyDuplicated(colnames(lt)))
lr <- rd(file.path(p3,'tables/supported_LR_edges.tsv'))
old <- rd(file.path(p3,'tables/ligand_activities.tsv'));old_status <- rd(file.path(p3,'tables/geneset_status.tsv'))
activities <- list();statuses <- list();coverage <- list();sets <- list();reproduction <- list()
make_sets <- function(de,background,test,direction,exclude=character()) {
 bg <- setdiff(background,exclude)
 p <- if(test=='FDR05')de$FDR else de$PValue
 keep <- is.finite(p)&p<.05&(if(direction=='Up')de$logFC>=.5 else de$logFC<= -.5)
 list(bg=bg,gs=intersect(de$gene[keep],bg))
}
for(receiver in c('Core','ActCD4'))for(filt in c('All','BG1')) {
 take <- if(receiver=='Core')rep(TRUE,nrow(m)) else m$subtype==receiver
 if(filt=='BG1')take <- take&m$bg<=.01
 mm <- m[take,,drop=FALSE];xx <- x[,take,drop=FALSE]
 ss <- sort(unique(m$sample));n <- as.numeric(table(factor(mm$sample,levels=ss)))
 agg <- sparseMatrix(i=seq_len(nrow(mm)),j=match(mm$sample,ss),x=1,dims=c(nrow(mm),length(ss)))
 pct <- sweep(as.matrix((xx>0)%*%agg),2,pmax(n,1),'/')
 stopifnot(all(n>=20),sum(grepl('^Saline_',ss))==3,sum(grepl('^Silica_',ss))==3)
 de <- rd(file.path(p3,'tables',paste0('DE_CD4_',receiver,'_',filt,'.tsv')))
 potential <- sort(intersect(unique(lr$from[lr$receiver==receiver & lr$filter==filt & lr$threshold=='Strict10']),colnames(lt)))
 for(cut in c(.1,.05)) {
   gate <- if(cut==.1)'Pct10' else 'Pct5'
   ec <- rownames(pct)[rowSums(pct[,grepl('^Saline_',ss),drop=FALSE]>=cut)>=2]
   ed <- rownames(pct)[rowSums(pct[,grepl('^Silica_',ss),drop=FALSE]>=cut)>=2]
   background <- Reduce(intersect,list(de$gene,union(ec,ed),rownames(lt)))
   coverage[[length(coverage)+1]] <- data.frame(receiver,filter=filt,gate,target=targets,in_background=targets%in%canon(background),DE_tested=targets%in%canon(de$gene))
   for(test in c('FDR05','Nominal05'))for(direction in c('Up','Down'))for(omit in c('Full','Drop7')) {
     exclude <- if(omit=='Drop7')rownames(lt)[canon(rownames(lt))%in%targets] else character()
     z <- make_sets(de,background,test,direction,exclude);gs <- z$gs;bg <- z$bg
     st <- data.frame(receiver,filter=filt,gate,test,direction,omit,background=length(bg),geneset=length(gs),ligands=length(potential),status='PASS')
     if(length(bg)<500)st$status <- 'SKIP_SMALL_BACKGROUND'
     if(length(gs)<10)st$status <- 'SKIP_GENESET_LT10'
     if(length(setdiff(bg,gs))<10)st$status <- 'SKIP_NEGATIVE_GENES_LT10'
     if(length(potential)<2)st$status <- 'SKIP_LIGANDS_LT2'
     statuses[[length(statuses)+1]] <- st
     if(length(gs))sets[[length(sets)+1]] <- data.frame(receiver,filter=filt,gate,test,direction,omit,gene=gs)
     if(gate=='Pct10'&&omit=='Full') {
       os <- old_status[old_status$receiver==receiver & old_status$filter==filt & old_status$test==test & old_status$direction==direction,,drop=FALSE]
       stopifnot(nrow(os)==1,st$background==os$background,st$geneset==os$geneset,st$ligands==os$ligands,st$status==os$status)
     }
     if(st$status!='PASS')next
     # Reuse verified original 10% scores rather than recomputing them.
     if(gate=='Pct10'&&omit=='Full') {
       act <- old[old$receiver==receiver & old$filter==filt & old$test==test & old$direction==direction,,drop=FALSE]
       stopifnot(nrow(act)==length(potential),setequal(act$test_ligand,potential))
       act <- act[,c('test_ligand','auroc','aupr','aupr_corrected','pearson','rank'),drop=FALSE]
       reproduction[[length(reproduction)+1]] <- data.frame(receiver,filter=filt,test,direction,status='GENESET_BACKGROUND_LIGAND_EXACT_MATCH_ORIGINAL_SCORES_REUSED')
     } else {
       act <- as.data.frame(nichenetr::predict_ligand_activities(geneset=gs,background_expressed_genes=bg,ligand_target_matrix=lt,potential_ligands=potential))
       stopifnot(all(c('test_ligand','aupr_corrected')%in%names(act)),any(is.finite(act$aupr_corrected)))
       act$rank <- rank(-act$aupr_corrected,ties.method='min')
       act <- act[,c('test_ligand','auroc','aupr','aupr_corrected','pearson','rank'),drop=FALSE]
     }
     act$receiver <- receiver;act$filter <- filt;act$gate <- gate;act$test <- test;act$direction <- direction;act$omit <- omit;act$geneset_n <- length(gs);act$background_n <- length(bg)
     activities[[length(activities)+1]] <- act
     cat('SENSITIVITY_DONE=',receiver,' ',filt,' ',gate,' ',test,' ',direction,' ',omit,' genes=',length(gs),'\n',sep='')
   }
 }
}
wr(do.call(rbind,activities),'ligand_activities.tsv');wr(do.call(rbind,statuses),'geneset_status.tsv');wr(do.call(rbind,coverage),'candidate_gate_coverage.tsv');wr(do.call(rbind,sets),'geneset_members.tsv');wr(do.call(rbind,reproduction),'baseline_reproduction.tsv')
writeLines(c('Step5: predefined 5% and10% target-expression gates, keeping >=2 libraries in either group, independent-library DE and sender/receptor eligibility unchanged. All and BG1 are separate.',
'Full uses all genes passing the frozen rule. Drop7 removes IL1R2/PENK/TTC39C/FOXP1/TNFRSF9/LAG3/MRE11 from BOTH the gene set and background; evaluates whether ligand ranking depends on the selected candidates. It is within-dataset sensitivity, not independent validation.',
'Only >=10 genes and>=500 background genes are scored. FDR and nominal P sets stay separately labeled. Failed FDR sets are not replaced by nominal sets. This is a prespecified follow-up prompted by observed low detection, not a new primary test.',
'Original 10% Full gene sets, backgrounds, ligand counts and statuses are checked against step3; their scores are reused unchanged. All other viable combinations are recomputed with the same NicheNet mouse prior.',
'Ranks are relative prioritization, not interaction P values. No new DE or causal sign is inferred. The original six-library design limits precision of sample correlations.',
'Output flags IL15/TNF plus all other ligands; no favorable threshold selection. No new database download.'),file.path(out,'methods.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
cat('ANALYSIS_LOGIC_PASS FROZEN_DE_AND_PRIOR_BASELINE_MATCHED\n')

