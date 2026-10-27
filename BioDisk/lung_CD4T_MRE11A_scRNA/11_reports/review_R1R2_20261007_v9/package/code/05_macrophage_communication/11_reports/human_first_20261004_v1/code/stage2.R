P <- '/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA'
O <- file.path(P,'11_reports/human_first_20261004_v1'); A <- file.path(P,'02_analysis')
suppressPackageStartupMessages({library(Matrix);library(edgeR);library(nichenetr)})
set.seed(20261004)
stopifnot(file.exists(file.path(O,'provenance/STAGE2_SPEC.txt')), !file.exists(file.path(O,'provenance/stage2_complete.txt')))
rd <- function(p) read.delim(p,check.names=FALSE,stringsAsFactors=FALSE)
wr <- function(d,n) {p<-file.path(O,'results',n);if(grepl('gz$',n)){con<-gzfile(p,'wt');on.exit(close(con));write.table(d,con,sep='\t',quote=FALSE,row.names=FALSE,na='NA')}else write.table(d,p,sep='\t',quote=FALSE,row.names=FALSE,na='NA')}
bind <- function(x) if(length(x)) do.call(rbind,x) else data.frame()
load <- function(pre,gcol,lineage){
 m<-rd(paste0(pre,'cells.tsv'));g<-rd(paste0(pre,'genes.tsv'))[[gcol]];x<-as(readMM(paste0(pre,'counts.mtx')),'CsparseMatrix')
 stopifnot(ncol(x)==nrow(m),nrow(x)==length(g),!anyDuplicated(m$cell),all(x@x>=0),all(x@x==round(x@x)))
 sym<-as.character(convert_alias_to_symbols(g,'human',verbose=FALSE));bad<-is.na(sym)|!nzchar(sym);sym[bad]<-g[bad]
 u<-sort(unique(sym));pr<-sparseMatrix(i=match(sym,u),j=seq_along(sym),x=1,dims=c(length(u),length(sym)));x<-pr%*%x;rownames(x)<-u;colnames(x)<-m$cell
 wr(data.frame(lineage,original=g,mapped=sym,duplicate_after_mapping=duplicated(sym)|duplicated(sym,fromLast=TRUE)) ,paste0('stage2_',lineage,'_mapping.tsv'))
 stopifnot(all(abs(colSums(x)-m$library)<1e-6));list(x=x,m=m)
}
cd<-load(file.path(A,'mechanism_step2_20260927_154919_EsN4gk/input/human_'),'symbol','CD4')
mac<-load(file.path(A,'scenic3_20261001_123452_Kvp9h5/input/human_mac_'),'gene','Mac')
rev<-rd(file.path(P,'11_reports/technical_delivery_20261002_ikzf2_stat1_v1/revision_v2/final_release/prior_review_repaired/scenic3_review/results/reviewed_cells.tsv.gz'))
stopifnot(setequal(rev$cell,mac$m$cell));mac$m$review_scope<-rev$review_scope[match(mac$m$cell,rev$cell)]
bgpanel<-strsplit('SFTPA1 SFTPA2 SFTPB SFTPC SFTPD SCGB1A1 SCGB3A1 SCGB3A2 ABCA3 SLC34A2 NAPSA LYZ LYZ1 LYZ2 S100A8 S100A9 CSF1R LST1 TYROBP FCER1G',' ')[[1]]
bm<-rownames(cd$x)%in%bgpanel|grepl('^(IGH[VDJMAEG]|IGK[VJC]|IGL[VJC])',rownames(cd$x));cd$m$bg<-colSums(cd$x[bm,,drop=FALSE])/cd$m$library
design<-rd(file.path(O,'results/stage1_donor_design.tsv'));samples<-sort(design$sample);design<-design[match(samples,design$sample),]
scopes<-c('Primary','Core',sort(unique(cd$m$subtype)));senders<-c('MacAM','MonoMac','MacWide')
aggregate_scope<-function(obj,scope,filt,lineage){
 m<-obj$m;x<-obj$x
 if(lineage=='CD4') take<-if(scope=='Primary')rep(TRUE,nrow(m)) else if(scope=='Core')m$core else m$subtype==scope else take<-if(scope=='MacWide')m$review_scope%in%c('MacAM','MonoMac') else m$review_scope==scope
 take<-take&m$sample%in%samples;if(filt=='BG1')take<-take&m$bg<=.01
 x<-x[,take,drop=FALSE];m<-m[take,,drop=FALSE]
 ag<-sparseMatrix(i=seq_len(nrow(m)),j=match(m$sample,samples),x=1,dims=c(nrow(m),length(samples)))
 n<-as.integer(table(factor(m$sample,levels=samples)));lib<-as.numeric(tapply(m$library,factor(m$sample,levels=samples),sum));lib[is.na(lib)]<-0
 pb<-as.matrix(x%*%ag);pct<-sweep(as.matrix((x>0)%*%ag),2,pmax(n,1),'/');avg<-sweep(as.matrix((x%*%Diagonal(x=1e4/m$library))%*%ag),2,pmax(n,1),'/')
 logx<-x%*%Diagonal(x=1e4/m$library);logx@x<-log1p(logx@x);logavg<-sweep(as.matrix(logx%*%ag),2,pmax(n,1),'/')
 pct[,n==0]<-avg[,n==0]<-logavg[,n==0]<-NA_real_
 colnames(pb)<-colnames(pct)<-colnames(avg)<-colnames(logavg)<-samples
 sm<-data.frame(sample=samples,group=design$group,series=design$series,cells=n,library=lib,mean_depth=lib/pmax(n,1),eligible20=n>=20,scope,filter=filt,lineage)
 list(pb=pb,pct=pct,avg=avg,logavg=logavg,sm=sm,take=take,scope=scope,filter=filt)
}
sc<-list();for(s in scopes)for(f in c('All','BG1'))sc[[paste(s,f,sep='_')]]<-aggregate_scope(cd,s,f,'CD4')
ms<-setNames(lapply(senders,function(s)aggregate_scope(mac,s,'All','Mac')),senders)
wr(bind(lapply(c(sc,ms),`[[`,'sm')),'stage2_scope_sample_counts.tsv')
saveRDS(list(receiver=sc,sender=ms),file.path(O,'results/stage2_aggregates.rds'))
cat('AGGREGATES_READY\n');flush.console()
lr0<-as.data.frame(readRDS('/mnt/BioDisk/wangmingyuan/IL17_Candida_integrated_v1/07_synthesis/supplementary_bioinformatics_v1/a6_consensus_regulatory_chain_v1/01_nichenet_resources_v1/lr_network.rds'));lr<-unique(lr0[,c('from','to')]);lr$from<-as.character(lr$from);lr$to<-as.character(lr$to)
lt<-readRDS('/mnt/BioDisk/wangmingyuan/IL17_Candida_integrated_v1/07_synthesis/supplementary_bioinformatics_v1/a6_consensus_regulatory_chain_v1/01_nichenet_resources_v1/ligand_target_matrix.rds')
wr(data.frame(resource=c('human_LT_targets','human_LT_ligands','human_LR_pairs'),n=c(nrow(lt),ncol(lt),nrow(lr))),'stage2_prior_dimensions.tsv')
expr<-function(z,group,cut=.1){j<-which(z$sm$eligible20&z$sm$group==group);if(length(j)<2)return(character());rownames(z$pct)[rowSums(z$pct[,j,drop=FALSE]>=cut,na.rm=TRUE)>=2]}
countsupport<-function(z,genes,group,cut=.1){j<-which(z$sm$eligible20&z$sm$group==group);v<-rep(0,length(genes));i<-genes%in%rownames(z$pct);if(length(j))v[i]<-rowSums(z$pct[genes[i],j,drop=FALSE]>=cut,na.rm=TRUE);v}
ex<-list();for(key in names(c(sc,ms))){z<-c(sc,ms)[[key]];gs<-intersect(if(key%in%senders)lr$from else lr$to,rownames(z$pb));for(j in seq_along(samples))ex[[length(ex)+1]]<-data.frame(key,gene=gs,sample=samples[j],group=z$sm$group[j],cells=z$sm$cells[j],pct=z$pct[gs,j],mean_CP10k=z$avg[gs,j])}
wr(bind(ex),'stage2_LR_expression_all_donors.tsv.gz');rm(ex)
edges<-list();potential<-list()
for(key in names(sc))for(sen in senders)for(cut in c(.1,.05)){
 z<-sc[[key]];s<-ms[[sen]];v<-lr
 v$receiver<-z$scope;v$filter<-z$filter;v$sender<-sen;v$gate<-if(cut==.1)'Pct10' else 'Pct5'
 v$ligand_present<-v$from%in%rownames(s$pb);v$receptor_present<-v$to%in%rownames(z$pb);v$ligand_prior<-v$from%in%colnames(lt)
 v$ligand_case_donors<-countsupport(s,v$from,'CWP',cut);v$receptor_case_donors<-countsupport(z,v$to,'CWP',cut)
 v$pair_supported<-v$ligand_case_donors>=2&v$receptor_case_donors>=2&v$ligand_prior
 v$status<-ifelse(!v$ligand_present|!v$receptor_present,'FEATURE_NOT_MAPPED',ifelse(!v$ligand_prior,'LIGAND_NOT_IN_LT_PRIOR',ifelse(v$pair_supported,'PAIR_EXPRESSION_COMPATIBLE','BELOW_EXPRESSION_GATE_NOT_BIOLOGICAL_ABSENCE')))
 v$complex_status<-'PAIR_ONLY_COMPLEX_NOT_EXHAUSTIVELY_CURATED';edges[[length(edges)+1]]<-v
 if(cut==.1)potential[[paste(key,sen,sep='|')]]<-sort(unique(v$from[v$pair_supported]))
}
edges<-bind(edges);wr(edges,'stage2_all_LR_supported_and_unsupported.tsv.gz')
# Fixed complex audits: component expression and same-cell co-detection are separate from pair support.
complexes<-list(IL18=c('IL18R1','IL18RAP'),IL1A=c('IL1R1','IL1RAP'),IL1B=c('IL1R1','IL1RAP'),IL2=c('IL2RB','IL2RG'),IL4=c('IL4R','IL2RG'),IL7=c('IL7R','IL2RG'),IL9=c('IL9R','IL2RG'),IL15=c('IL2RB','IL2RG'),IL21=c('IL21R','IL2RG'),IL10=c('IL10RA','IL10RB'),TGFB1=c('TGFBR1','TGFBR2'),TGFB2=c('TGFBR1','TGFBR2'),TGFB3=c('TGFBR1','TGFBR2'),IL6=c('IL6R','IL6ST'),IL12A=c('IL12RB1','IL12RB2'),IL12B=c('IL12RB1','IL12RB2'),IL23A=c('IL23R','IL12RB1'),IFNG=c('IFNGR1','IFNGR2'),CSF2=c('CSF2RA','CSF2RB'))
cr<-list();for(key in names(sc))for(lig in names(complexes)){
 z<-sc[[key]];gs<-complexes[[lig]];present<-all(gs%in%rownames(cd$x));marg<-if(present)all(vapply(gs,function(g)g%in%expr(z,'CWP'),logical(1)))else FALSE
 for(ss in samples){ii<-which(z$take&cd$m$sample==ss);n<-length(ii);co<-if(present&&n)sum(colSums(cd$x[gs,ii,drop=FALSE]>0)==length(gs))else NA_integer_;cr[[length(cr)+1]]<-data.frame(receiver=z$scope,filter=z$filter,ligand=lig,sample=ss,cells=n,chains=paste(gs,collapse=';'),features_present=present,all_chains_marginal_gate=marg,codetected=co,pct=co/max(n,1),note=if(lig=='IL15')'IL15RA trans-presentation not established' else 'RNA component audit; no protein/function conclusion')}
}
wr(bind(cr),'stage2_receptor_complex_audit.tsv.gz')
de_status<-list();sets<-list();bgs<-list();set_status<-list();acts<-list();candidates<-list();programs<-list();links<-list()
for(key in names(sc)){
 z<-sc[[key]];sel<-z$sm$eligible20;md<-z$sm[sel,,drop=FALSE];md$group<-factor(md$group,levels=c('Control','CWP'));md$series<-factor(md$series)
 for(model in c('Base','Source')){
  st<-data.frame(receiver=z$scope,filter=z$filter,model,n_control=sum(md$group=='Control'),n_case=sum(md$group=='CWP'),samples=paste(md$sample,collapse=';'),df=NA_integer_,tested=0,status='INSUFFICIENT_GROUP_REPLICATION')
  if(min(st$n_control,st$n_case)<2){de_status[[length(de_status)+1]]<-st;next}
  X<-if(model=='Source'&&nlevels(md$series)>1)model.matrix(~group+series,md)else model.matrix(~group,md);st$df<-nrow(X)-qr(X)$rank
  if(qr(X)$rank<ncol(X)||st$df<2){st$status<-'NON_IDENTIFIABLE';de_status[[length(de_status)+1]]<-st;next}
  y<-DGEList(z$pb[,sel,drop=FALSE],lib.size=md$library);keep<-filterByExpr(y,design=X)
  if(sum(keep)<100){st$status<-'TOO_FEW_EXPRESSED_GENES';de_status[[length(de_status)+1]]<-st;next}
  y<-y[keep,,keep.lib.sizes=TRUE];y<-calcNormFactors(y);y<-estimateDisp(y,X,robust=TRUE);fit<-glmQLFit(y,X,robust=TRUE)
  tt<-topTags(glmQLFTest(fit,coef=which(colnames(X)=='groupCWP')),n=Inf,sort.by='none')$table
  d<-data.frame(gene=rownames(tt),tt,receiver=z$scope,filter=z$filter,model,check.names=FALSE);wr(d,paste0('stage2_DE_',key,'_',model,'.tsv.gz'))
  st$status<-'FIT';st$tested<-nrow(d);de_status[[length(de_status)+1]]<-st
  bg<-Reduce(intersect,list(d$gene,union(expr(z,'Control'),expr(z,'CWP')),rownames(lt)));bgs[[length(bgs)+1]]<-data.frame(receiver=z$scope,filter=z$filter,model,gene=bg)
  for(test in c('FDR05','Nominal05'))for(direction in c('Up','Down')){
   pv<-if(test=='FDR05')d$FDR else d$PValue;gs<-intersect(d$gene[pv<.05 & if(direction=='Up')d$logFC>=.5 else d$logFC<= -.5],bg)
   if(length(gs))sets[[length(sets)+1]]<-data.frame(receiver=z$scope,filter=z$filter,model,test,direction,gene=gs)
   pp<-sort(unique(unlist(potential[paste(key,senders,sep='|')])));status<-if(length(bg)<500)'BACKGROUND_LT500' else if(length(gs)<10)'TARGETSET_LT10' else if(length(setdiff(bg,gs))<10)'NEGATIVESET_LT10' else if(length(pp)<2)'LIGANDS_LT2' else 'COMPUTED'
   set_status[[length(set_status)+1]]<-data.frame(receiver=z$scope,filter=z$filter,model,test,direction,background_n=length(bg),geneset_n=length(gs),candidate_n=length(pp),status)
   if(length(gs))programs[[length(programs)+1]]<-data.frame(receiver=z$scope,filter=z$filter,model,test,direction,sample=samples,target_score=colMeans(z$logavg[gs,,drop=FALSE]),cells=z$sm$cells)
   ac<-NULL
   if(status=='COMPUTED'){
    ac<-as.data.frame(predict_ligand_activities(geneset=gs,background_expressed_genes=bg,ligand_target_matrix=lt,potential_ligands=pp));ac$rank<-rank(-ac$aupr_corrected,ties.method='min');ac$receiver<-z$scope;ac$filter<-z$filter;ac$model<-model;ac$test<-test;ac$direction<-direction;ac$geneset_n<-length(gs);ac$background_n<-length(bg);acts[[length(acts)+1]]<-ac
    con<-gzfile(file.path(O,'results',paste0('stage2_LT_',key,'_',model,'_',test,'_',direction,'.tsv.gz')),'wt')
    ll<-expand.grid(target=gs,ligand=pp,stringsAsFactors=FALSE);ll$weight<-as.numeric(lt[gs,pp,drop=FALSE]);ll$receiver<-z$scope;ll$filter<-z$filter;ll$model<-model;ll$test<-test;ll$direction<-direction;write.table(ll,con,sep='\t',quote=FALSE,row.names=FALSE);close(con)
   }
   for(sen in senders){v<-data.frame(ligand=colnames(lt),receiver=z$scope,filter=z$filter,model,test,direction,sender=sen);v$LR_expression_supported<-v$ligand%in%potential[[paste(key,sen,sep='|')]];v$geneset_n<-length(gs);v$background_n<-length(bg);v$aupr_corrected<-if(is.null(ac))NA_real_ else ac$aupr_corrected[match(v$ligand,ac$test_ligand)];v$rank<-if(is.null(ac))NA_real_ else ac$rank[match(v$ligand,ac$test_ligand)];v$status<-ifelse(!v$LR_expression_supported,'NO_SOURCE_LR_GATE_SUPPORT',if(status!='COMPUTED')rep(status,nrow(v)) else ifelse(is.finite(v$aupr_corrected)&v$aupr_corrected>0,'POSITIVE_ACTIVITY_PREDICTION','COMPUTED_NONPOSITIVE_ACTIVITY'));candidates[[length(candidates)+1]]<-v}
  }
  cat('DE_NICHENET_DONE',key,model,'genes',nrow(d),'\n');flush.console()
 }
}
wr(bind(de_status),'stage2_DE_status.tsv');wr(bind(bgs),'stage2_background_members.tsv.gz');wr(bind(sets),'stage2_targetset_members.tsv.gz');wr(bind(set_status),'stage2_targetset_status.tsv');wr(bind(acts),'stage2_ligand_activities.tsv.gz');wr(bind(candidates),'stage2_all_candidates_with_negative_results.tsv.gz');wr(bind(programs),'stage2_disease_program_donor_scores.tsv.gz')
writeLines(capture.output(sessionInfo()),file.path(O,'provenance/stage2_R_sessionInfo.txt'));writeLines('COMPLETE',file.path(O,'provenance/stage2_complete.txt'));cat('STAGE2_COMPLETE\n')
