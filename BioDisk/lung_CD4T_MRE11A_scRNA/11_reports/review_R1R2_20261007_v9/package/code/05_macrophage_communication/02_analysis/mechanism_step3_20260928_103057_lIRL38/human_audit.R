args <- commandArgs(TRUE);prev2 <- args[1];out <- args[2]
rd <- function(f)read.delim(f,check.names=FALSE,stringsAsFactors=FALSE)
wr <- function(x,n)write.table(x,file.path(out,'tables',n),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
t <- rd(file.path(prev2,'tables/human_sample_scores.tsv'))
# Correct the previous reporting gate: BH is mathematically valid with fewer than10 samples/group.
m <- rd(file.path(prev2,'tables/human_source_models.tsv'));m$FDR <- NA_real_
for(ix in split(seq_len(nrow(m)),interaction(m$scope,m$filter,m$model,m$family,drop=TRUE))) {
 j <- ix[is.finite(m$p_nominal[ix])];if(length(j))m$FDR[j] <- p.adjust(m$p_nominal[j],'BH')
}
m$FDR_status <- 'BH_ALL_FINITE_WITHIN_SCOPE_FILTER_MODEL_FAMILY'
wr(m,'human_source_models_BH.tsv')
calc <- function(v,g,mod,adjust) {
 base <- c('group',if(length(unique(v$series))>1)'series','depth_rank')
 if(!is.null(adjust))base <- c(base,adjust)
 v$depth_rank <- rank(v$library);v$group <- factor(v$group);v$series <- factor(v$series)
 X <- model.matrix(reformulate(base),v);df <- nrow(v)-qr(X)$rank-1
 if(qr(X)$rank<ncol(X)||df<3)return(c(rho=NA,p=NA,df=df))
 u <- qr.resid(qr(X),rank(v[[paste0('expr_',g)]]));w <- qr.resid(qr(X),rank(v[[paste0('score_',mod)]]))
 if(sd(u)<1e-10||sd(w)<1e-10)return(c(rho=NA,p=NA,df=df))
 r <- max(-1,min(1,cor(u,w)));p <- if(abs(r)>=1)0 else 2*pt(-abs(r)*sqrt(df/(1-r*r)),df)
 c(rho=r,p=p,df=df)
}
res <- list();loo <- list()
for(filt in c('All','BG1')) {
 v <- t[t$scope=='Core' & t$filter==filt & t$cells>=20,,drop=FALSE]
 # Use all observed cells in each Core subtype, including <20-cell subtypes, to compute proportions.
 sub <- t[grepl('^Subtype:',t$scope)&t$filter==filt & t$scope!='Subtype:CD4state',,drop=FALSE]
 for(label in unique(sub$scope)) {
   q <- sub[sub$scope==label,,drop=FALSE];cv <- paste0('frac_',make.names(sub('^Subtype:','',label)))
   v[[cv]] <- q$cells[match(v$sample,q$sample)];v[[cv]][is.na(v[[cv]])] <- 0;v[[cv]] <- v[[cv]]/v$cells
 }
 vars <- grep('^frac_',names(v),value=TRUE)
 stopifnot(max(abs(rowSums(v[,vars,drop=FALSE])-1))<1e-8)
 for(g in c('IL1R2','PENK','TTC39C','FOXP1','TNFRSF9','LAG3','MRE11'))for(mod in sub('^score_','',grep('^score_',names(t),value=TRUE)))for(adj in c('Base','JointComp',vars)) {
  w <- v[is.finite(v[[paste0('expr_',g)]]) & is.finite(v[[paste0('score_',mod)]]),,drop=FALSE]
  rr <- calc(w,g,mod,if(adj=='Base')NULL else if(adj=='JointComp')head(vars,-1) else adj)
  res[[length(res)+1]] <- data.frame(filter=filt,gene=g,module=mod,adjust=adj,n=nrow(w),rho=rr['rho'],p_nominal=rr['p'],df=rr['df'])
  if(adj=='Base')for(i in seq_len(nrow(w))) {
    r <- calc(w[-i,,drop=FALSE],g,mod,NULL)
    loo[[length(loo)+1]] <- data.frame(filter=filt,gene=g,module=mod,removed_sample=w$sample[i],rho=r['rho'],full_rho=rr['rho'])
  }
 }
}
r <- do.call(rbind,res);r$FDR <- NA_real_
for(ix in split(seq_len(nrow(r)),interaction(r$filter,r$adjust,drop=TRUE))) {
 j <- ix[is.finite(r$p_nominal[ix])];if(length(j))r$FDR[j] <- p.adjust(r$p_nominal[j],'BH')
}
wr(r,'human_composition_correlations.tsv');wr(do.call(rbind,loo),'human_correlation_LOO.tsv')
cat('HUMAN_CORRELATION_AUDIT_PASS\n')

