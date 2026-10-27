args <- commandArgs(TRUE);p1 <- args[1];p3 <- args[2];out <- args[3]
suppressPackageStartupMessages(library(edgeR))
rd <- function(p)read.delim(p,check.names=FALSE,stringsAsFactors=FALSE)
wr <- function(x,n)write.table(x,file.path(out,'tables',paste0(n,'.tsv')),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
canon <- function(x){x <- toupper(trimws(x));a <- c(MRE11A='MRE11',ATP5A1='ATP5F1A',ATP5B='ATP5F1B',ATP5C1='ATP5F1C',ATP5G1='ATP5MC1',ATP5O='ATP5PO',PARK2='PRKN');i <- x%in%names(a);x[i] <- unname(a[x[i]]);x}
targets <- c('IL1R2','PENK','TTC39C','FOXP1','TNFRSF9','LAG3','MRE11')
# ULM with intercept, using all expressed genes; unlinked weights are zero.
# Matches t value of the weight coefficient in lm(response ~ weight).
ulm <- function(expr,network,exclude=character(),minsize=10) {
 expr <- expr[!rownames(expr)%in%exclude,,drop=FALSE]
 if(nrow(expr)<500)return(NULL)
 net <- network[network$target%in%rownames(expr),,drop=FALSE]
 sz <- table(net$tf);tfs <- names(sz)[sz>=minsize]
 if(!length(tfs))return(NULL)
 ans <- matrix(NA_real_,length(tfs),ncol(expr),dimnames=list(tfs,colnames(expr)))
 Y <- sweep(expr,2,colMeans(expr),'-');sy <- colSums(Y^2);n <- nrow(expr)
 for(tf in tfs){q <- net[net$tf==tf,,drop=FALSE];w <- numeric(n);w[match(q$target,rownames(expr))] <- q$mor;w <- w-mean(w);sw <- sum(w*w)
  if(sw<1e-12)next
  r <- as.numeric(crossprod(w,Y))/sqrt(sw*sy);r[sy<1e-12] <- NA_real_;r <- pmax(-1+1e-12,pmin(1-1e-12,r));ans[tf,] <- r*sqrt((n-2)/(1-r*r))
 }
 list(scores=ans,sizes=as.integer(sz[tfs]),universe=n)
}
partial <- function(x,y,group) {
 ok <- is.finite(x)&is.finite(y);x <- x[ok];y <- y[ok];g <- group[ok]
 if(length(x)<5||length(unique(g))!=2||min(table(g))<2)return(NA_real_)
 Z <- model.matrix(~factor(g));u <- lm.fit(Z,rank(x))$residuals;v <- lm.fit(Z,rank(y))$residuals
 if(sum(u*u)<1e-10||sum(v*v)<1e-10)return(NA_real_)
 max(-1,min(1,cor(u,v)))
}
assoc <- function(x,y,g){r <- partial(x,y,g);loo <- vapply(seq_along(x),function(i)partial(x[-i],y[-i],g[-i]),numeric(1));ok <- is.finite(loo);c(rho=r,LOO_min=if(any(ok))min(loo[ok]) else NA_real_,LOO_max=if(any(ok))max(loo[ok]) else NA_real_,LOO_n=sum(ok))}
effect <- function(x,g){a <- x[g=='Saline'];b <- x[g=='Silica'];delta <- mean(b)-mean(a);loo <- vapply(seq_along(x),function(i)mean(x[-i][g[-i]=='Silica'])-mean(x[-i][g[-i]=='Saline']),numeric(1));c(control=mean(a),case=mean(b),delta=delta,LOO_min=min(loo),LOO_max=max(loo))}
net <- rd(file.path(out,'network.tsv'));stopifnot(all(c('tf','target','mor','confidence')%in%names(net)),all(net$confidence%in%c('A','B','C')))
net$tf_original <- net$tf;net$target_original <- net$target;net$tf <- canon(net$tf);net$target <- canon(net$target)
stopifnot(all(is.finite(net$mor)),all(net$mor!=0))
# Conflicting signs are excluded, identical duplicate edges retain best confidence.
key <- paste(net$tf,net$target,sep=':');bad <- names(which(vapply(split(sign(net$mor),key),function(x)length(unique(x))>1,logical(1))))
wr(net[key%in%bad,,drop=FALSE],'network_conflicting_edges')
net <- net[!key%in%bad,];net <- net[order(net$confidence,net$tf,net$target),];net <- net[!duplicated(net[,c('tf','target')]),]
wr(net,'network_canonical')
fits <- readRDS(file.path(p3,'objects/pseudobulk_models.rds'));stopifnot(all(c('CD4_Core_All','CD4_ActCD4_All','CD4_Core_BG1','CD4_ActCD4_BG1')%in%names(fits)))
mods <- rd(file.path(p1,'tables/module_definition.tsv'));mods$gene <- canon(mods$gene)
ss <- rd(file.path(p1,'tables/sample_scores.tsv'));ss <- ss[ss$dataset=='GSE314943',]
mac <- rd(file.path(p3,'tables/ligand_receptor_sample_expression.tsv'))
stopifnot(all(c('lineage','scope','filter','canonical','sample','group','cells','mean_CP10k')%in%names(mac)))
mac <- mac[mac$lineage=='Mac'&mac$scope=='AMlike'&mac$filter=='All'&mac$canonical%in%c('IL15','TNF'),]
mac$gene <- as.character(mac$canonical);mac$CP10k <- mac$mean_CP10k
expected_samples <- fits[['CD4_ActCD4_All']]$metadata$sample
stopifnot(length(expected_samples)==6,!anyDuplicated(mac[,c('sample','gene')]))
for(lig in c('IL15','TNF')){
 mm <- mac[mac$gene==lig,]
 missing <- setdiff(expected_samples,mm$sample)
 if(length(missing))stop('Missing AMlike ',lig,' samples in step3 expression table: ',paste(missing,collapse=','))
 stopifnot(nrow(mm)==6,setequal(mm$sample,expected_samples),all(mm$cells>=20),all(is.finite(mm$CP10k)),all(mm$CP10k>=0),all(mm$group==sub('_.*','',mm$sample)))
}
wr(mac,'AMlike_ligand_input_audit')
cat('ANALYSIS_LOGIC_PASS MAC_LIGAND_INPUT_2_GENES_6_SAMPLES_EACH\n')
res <- list();activity <- list();coverage <- list();effects <- list();associations <- list();edges <- list();scope_audit <- list();skips <- list();universe <- list()
for(key in names(fits)[grepl('^CD4_',names(fits))]){
 f <- fits[[key]];md <- f$metadata;scope <- unique(md$scope);filt <- unique(md$filter)
 stopifnot(length(scope)==1,length(filt)==1,nrow(md)==6,sum(md$group=='Saline')==3,sum(md$group=='Silica')==3,all(md$cells>=20),identical(colnames(f$y$counts),md$sample))
 # Use prior edgeR normalization, retaining original full-library offsets.
 E <- cpm(f$y,log=TRUE,prior.count=2);cg <- canon(rownames(E));stopifnot(!anyDuplicated(cg));rownames(E) <- cg
 # Gene-centered logCPM measures relative regulatory activity across six samples.
 E <- E-rowMeans(E)
 scope_audit[[key]] <- md;universe[[key]] <- data.frame(key,gene=rownames(E))
 q <- ss[ss$scope==if(scope=='Core')'All' else scope,];q <- q[q$filter==filt,];q <- q[match(md$sample,q$sample),]
 stopifnot(nrow(q)==6,!anyNA(q$sample),identical(q$sample,md$sample),all(q$cells==md$cells))
 de <- rd(file.path(p3,'tables',paste0('DE_',key,'.tsv')));de$canonical <- canon(de$gene)
 for(level in c('ABC','AB')){
  nn <- net[net$confidence%in%strsplit(level,'')[[1]],]
  for(mode in c('Full','Drop7')){
   rm <- if(mode=='Drop7')targets else character();z <- ulm(E,nn,rm)
   if(is.null(z)){skips[[length(skips)+1]] <- data.frame(key,level,mode,reason='NO_TF_WITH_10_TARGETS');next}
   A <- z$scores;stopifnot(all(is.finite(A)))
   for(tf in rownames(A)){
    x <- A[tf,];sz <- z$sizes[match(tf,rownames(A))]
    activity[[length(activity)+1]] <- data.frame(key,scope,filter=filt,level,mode,tf,sample=md$sample,group=md$group,cells=md$cells,activity=as.numeric(x),n_targets=sz)
    ef <- effect(x,md$group);effects[[length(effects)+1]] <- data.frame(key,scope,filter=filt,level,mode,tf,n_targets=sz,control=ef['control'],case=ef['case'],delta=ef['delta'],LOO_min=ef['LOO_min'],LOO_max=ef['LOO_max'],row.names=NULL)
   }
   coverage[[length(coverage)+1]] <- data.frame(key,level,mode,tf=rownames(A),n_targets=z$sizes,universe=z$universe)
   if(mode!='Drop7')next
   for(tf in rownames(A)){
    x <- A[tf,]
    for(gene in targets){v <- q[[paste0('expr_',gene)]];stopifnot(length(v)==6);r <- assoc(x,v,md$group);associations[[length(associations)+1]] <- data.frame(key,level,tf,kind='Candidate',feature=gene,exclusion='Drop7',rho=r['rho'],LOO_min=r['LOO_min'],LOO_max=r['LOO_max'],LOO_n=r['LOO_n'],row.names=NULL)}
    for(lig in c('IL15','TNF')){mm <- mac[mac$gene==lig,];mm <- mm[match(md$sample,mm$sample),];stopifnot(!anyNA(mm$sample),all(mm$cells>=20));r <- assoc(x,mm$CP10k,md$group);associations[[length(associations)+1]] <- data.frame(key,level,tf,kind='AMlike',feature=lig,exclusion='Drop7',rho=r['rho'],LOO_min=r['LOO_min'],LOO_max=r['LOO_max'],LOO_n=r['LOO_n'],row.names=NULL)}
   }
   for(mod in unique(mods$module)){
    mg <- unique(mods$gene[mods$module==mod]);zz <- ulm(E,nn,union(targets,mg));if(is.null(zz))next
    for(tf in rownames(zz$scores)){
     val <- q[[paste0('score_',mod)]];stopifnot(length(val)==6);r <- assoc(zz$scores[tf,],val,md$group)
     associations[[length(associations)+1]] <- data.frame(key,level,tf,kind='Module',feature=mod,exclusion='Drop7PlusModule',rho=r['rho'],LOO_min=r['LOO_min'],LOO_max=r['LOO_max'],LOO_n=r['LOO_n'],row.names=NULL)
     coverage[[length(coverage)+1]] <- data.frame(key,level,mode=paste0('Drop7Plus_',mod),tf,n_targets=zz$sizes[match(tf,rownames(zz$scores))],universe=zz$universe)
    }
   }
   # Prior TF->candidate and TF->module edges. These are neither causal nor candidate->TF edges.
   ee <- nn[nn$tf%in%rownames(A)&nn$target%in%unique(c(targets,mods$gene)),]
   if(nrow(ee)){
    di <- match(ee$target,de$canonical)
    ee$key <- key;ee$level <- level;ee$target_tested <- !is.na(di);ee$target_logFC <- de$logFC[di];ee$target_P <- de$PValue[di];ee$target_FDR <- de$FDR[di];ee$candidate <- ee$target%in%targets
    ee$modules <- vapply(ee$target,function(g)paste(unique(mods$module[mods$gene==g]),collapse=';'),character(1))
    edges[[length(edges)+1]] <- ee
   }
   res[[paste(key,level)]] <- list(scores=A,metadata=md,expression=E)
  }
 }
 cat('TF_SCOPE_DONE=',key,'\n',sep='')
}
bind <- function(x)if(length(x))do.call(rbind,x) else data.frame(note=character())
wr(bind(activity),'TF_sample_activity');wr(bind(effects),'TF_group_effects');wr(bind(coverage),'TF_coverage');wr(bind(associations),'TF_associations');wr(bind(edges),'TF_target_edges');wr(bind(scope_audit),'sample_audit');wr(bind(skips),'skipped');wr(bind(universe),'expression_universe')
# Candidate coverage audit explicitly distinguishes no prior edge from unexpressed target.
cc <- expand.grid(key=names(fits)[grepl('^CD4_',names(fits))],gene=targets,stringsAsFactors=FALSE)
cc$prior_TFs_ABC <- vapply(cc$gene,function(g)paste(sort(unique(net$tf[net$target==g])),collapse=';'),character(1))
cc$tested_expression <- vapply(seq_len(nrow(cc)),function(i)cc$gene[i]%in%canon(rownames(fits[[cc$key[i]]]$y$counts)),logical(1));wr(cc,'candidate_prior_coverage')
saveRDS(res,file.path(out,'objects/TF_scores_and_expression.rds'))
writeLines(c('Mouse GSE314943, independent-library pseudobulk, original step3 edgeR TMM/full-library offsets; 3 Saline vs 3 Silica, min20 cells. Each eligible CD4 scope analyzed separately; Core remains composition-sensitive.',
'DoRothEA mouse ABC primary, AB sensitivity. Uppercase canonical matching is within mouse, not a human orthology map. Opposing-sign duplicate edges excluded; identical duplicates keep best confidence. Exact prior snapshot and SHA256 supplied.',
'ULM implemented as t statistic of lm(centered logCPM ~ signed TF weights), intercept included, unlinked genes weight0, >=10 measured regulon targets. logCPM prior.count2; each gene centered across six samples within scope/filter. Scores describe relative regulatory activity, not TF RNA abundance or proof of binding.',
'Full scores provide baseline; Drop7 excludes all seven candidate genes from expression universe and regulons. Module correlations recompute activity after additionally excluding every member of the corresponding module from universe/regulon. Matched-control genes used in module scoring can still overlap regulons; correlations are exploratory.',
'Group differences and leave-one-sample-out directions are descriptive. Six samples: no TF disease-significance declarations or BH selection. ULM gene-level t statistics are activity scores, not biological-replicate P values. Rank correlations residualize condition; LOO refits correlation on remaining samples.',
'Candidate is not assumed to be a TF. Prior arrows are TF->target; candidate-TF and macrophage ligand-TF links are unsigned observational associations. IL15/TNF links are not receptor-to-TF signaling paths. Existing target DEG P/FDR retained for provenance.',
'No SCENIC, human macrophage validation, human TF analysis or final mechanism novelty audit is claimed completed by this stage.',
'https://saezlab.github.io/decoupleR/reference/run_ulm.html',
'https://saezlab.github.io/dorothea/articles/dorothea.html'),file.path(out,'methods.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
cat('TF_ANALYSIS_PASS\n')

