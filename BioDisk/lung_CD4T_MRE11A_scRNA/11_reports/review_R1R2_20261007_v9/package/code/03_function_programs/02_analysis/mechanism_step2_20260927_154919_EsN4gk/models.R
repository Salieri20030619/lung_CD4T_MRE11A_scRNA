out <- commandArgs(TRUE)[1]
t <- read.delim(file.path(out,'tables/human_sample_scores.tsv'),check.names=FALSE)
wr <- function(x,name) write.table(x,file.path(out,'tables',name),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
metrics <- grep('^(score_|expr_|pct_)',colnames(t),value=TRUE)
targets <- c('IL1R2','PENK','TTC39C','FOXP1','TNFRSF9','LAG3','MRE11')
mods <- sub('^score_','',grep('^score_',colnames(t),value=TRUE))
results <- list();assocs <- list();idx <- 0L;jdx <- 0L
for(scope in unique(t$scope)) for(filt in c('All','BG1')) {
  z <- t[t$scope==scope & t$filter==filt & t$cells>=20,,drop=FALSE]
  for(arm in c('Source','SourceDep','WithinG227')) {
    v <- if(arm=='WithinG227') z[z$series=='G227',,drop=FALSE] else z
    v$group <- factor(v$group,levels=c('Control','CWP'));v$series <- factor(v$series)
    for(metric in metrics) {
      w <- v[is.finite(v[[metric]]) & is.finite(v$library),,drop=FALSE]
      w <- droplevels(w);nc <- sum(w$group=='Control');nd <- sum(w$group=='CWP')
      row <- data.frame(scope=scope,filter=filt,model=arm,metric=metric,n_control=nc,n_case=nd,
        coefficient=NA_real_,se=NA_real_,CI_low=NA_real_,CI_high=NA_real_,p_nominal=NA_real_,residual_df=NA_real_,status='INSUFFICIENT_SAMPLES',stringsAsFactors=FALSE)
      if(nc>=2 && nd>=2) {
        terms <- c('group',if(arm!='WithinG227' && nlevels(w$series)>1) 'series',if(arm=='SourceDep') 'log1p(library)')
        form <- reformulate(terms,response=metric);X <- model.matrix(form,w)
        if(qr(X)$rank!=ncol(X) || nrow(X)-ncol(X)<2) row$status <- 'NON_IDENTIFIABLE' else {
          fit <- lm(form,data=w);co <- summary(fit)$coefficients
          row$status <- 'FIT'
          if('groupCWP'%in%rownames(co)) {
            b <- co['groupCWP',];df <- df.residual(fit);delta <- qt(.975,df)*b[2]
            row$coefficient <- b[1];row$se <- b[2];row$CI_low <- b[1]-delta;row$CI_high <- b[1]+delta;row$p_nominal <- b[4];row$residual_df <- df
          } else row$status <- 'GROUP_COEFFICIENT_MISSING'
        }
      }
      idx <- idx+1L;results[[idx]] <- row
    }
  }
  # Sample-level partial rank correlations adjust condition, source and mean library depth.
  for(g in targets) for(mod in mods) {
    ex <- paste0('expr_',g);sc <- paste0('score_',mod)
    v <- z[is.finite(z[[ex]]) & is.finite(z[[sc]]) & is.finite(z$library),,drop=FALSE]
    row <- data.frame(scope=scope,filter=filt,gene=g,module=mod,n=nrow(v),positive_samples=sum(v[[ex]]>0),rho=NA_real_,p_nominal=NA_real_,residual_df=NA_real_,status='INSUFFICIENT',stringsAsFactors=FALSE)
    if(nrow(v)>=10 && sum(v[[ex]]>0)>=4 && length(unique(v$group))==2 && sd(v[[ex]])>0 && sd(v[[sc]])>0) {
      v$group <- factor(v$group);v$series <- factor(v$series);v$depth_rank <- rank(v$library)
      terms <- c('group',if(nlevels(v$series)>1) 'series','depth_rank')
      X <- model.matrix(reformulate(terms),v);rankX <- qr(X)$rank;df <- nrow(v)-rankX-1L
      if(rankX==ncol(X) && df>=3) {
        a <- qr.resid(qr(X),rank(v[[ex]]));b <- qr.resid(qr(X),rank(v[[sc]]))
        if(sd(a)>1e-10 && sd(b)>1e-10) {
          r <- max(-1,min(1,cor(a,b)));row$rho <- r;row$residual_df <- df
          row$p_nominal <- if(abs(r)>=1) 0 else 2*pt(-abs(r)*sqrt(df/(1-r*r)),df)
          row$status <- 'EXPLORATORY_PARTIAL_RANK'
        }
      } else row$status <- 'NON_IDENTIFIABLE'
    }
    jdx <- jdx+1L;assocs[[jdx]] <- row
  }
}
r <- do.call(rbind,results);r$family <- sub('_.*$','',r$metric);r$FDR <- NA_real_
r$FDR_status <- 'NOT_COMPUTED_GROUP_N_LT10'
for(ix in split(seq_len(nrow(r)),interaction(r$scope,r$filter,r$model,r$family,drop=TRUE))) {
  keep <- ix[r$n_control[ix]>=10 & r$n_case[ix]>=10 & is.finite(r$p_nominal[ix])]
  if(length(keep)) {r$FDR[keep] <- p.adjust(r$p_nominal[keep],'BH');r$FDR_status[keep] <- 'BH_WITHIN_SCOPE_FILTER_MODEL_FAMILY'}
}
wr(r,'human_source_models.tsv');a <- do.call(rbind,assocs);a$FDR <- NA_real_
for(ix in split(seq_len(nrow(a)),interaction(a$scope,a$filter,drop=TRUE))) {
  keep <- ix[is.finite(a$p_nominal[ix])]
  if(length(keep)) a$FDR[keep] <- p.adjust(a$p_nominal[keep],'BH')
}
wr(a,'human_partial_correlations.tsv')
# Fixed-composition source-adjusted effects, explicitly retaining low-count subtype flags.
s <- read.delim(file.path(out,'tables/human_fixed_composition.tsv'),check.names=FALSE)
fres <- list();j <- 0L
for(ix in split(seq_len(nrow(s)),interaction(s$scope,s$filter,s$metric,drop=TRUE))) {
  v <- s[ix,,drop=FALSE];v <- v[is.finite(v$value),,drop=FALSE]
  if(!nrow(v))next
  row <- data.frame(scope=v$scope[1],filter=v$filter[1],metric=v$metric[1],n_control=sum(v$group=='Control'),n_case=sum(v$group=='CWP'),all_samples_subtypes_ge20=all(v$all_subtypes_eligible20),coefficient=NA_real_,status='INSUFFICIENT',stringsAsFactors=FALSE)
  if(row$n_control>=2 && row$n_case>=2) {
    v$group <- factor(v$group,levels=c('Control','CWP'));v$series <- factor(v$series)
    f <- reformulate(c('group',if(nlevels(v$series)>1)'series'),response='value');X <- model.matrix(f,v)
    if(qr(X)$rank==ncol(X) && nrow(X)>ncol(X)) {row$coefficient <- unname(coef(lm(f,v))['groupCWP']);row$status <- 'DESCRIPTIVE_FIXED_COMPOSITION'}
  }
  j <- j+1L;fres[[j]] <- row
}
if(length(fres))wr(do.call(rbind,fres),'human_fixed_composition_effects.tsv')
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
cat('HUMAN_SOURCE_MODELS_PASS\n')

