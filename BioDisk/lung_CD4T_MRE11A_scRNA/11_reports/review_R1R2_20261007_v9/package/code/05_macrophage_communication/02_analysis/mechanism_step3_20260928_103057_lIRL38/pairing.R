# Sourced after analyze.R: uses sc/tl generated above, independent-library descriptive association.
partial_group <- function(x,y,grp) {
 ok <- is.finite(x)&is.finite(y);x <- x[ok];y <- y[ok];grp <- grp[ok]
 if(length(x)<5||length(unique(grp))!=2||min(table(grp))<2||sd(x)==0||sd(y)==0)return(NA_real_)
 X <- model.matrix(~factor(grp));u <- qr.resid(qr(X),rank(x));v <- qr.resid(qr(X),rank(y))
 if(sd(u)<1e-10||sd(v)<1e-10)return(NA_real_)
 cor(u,v)
}
pr <- list()
if(nrow(tl)) {
 q <- unique(tl[tl$test=='FDR05' & tl$ligand_rank<=20,c('receiver','filter','ligand','target','canonical')])
 for(i in seq_len(nrow(q))) {
   rz <- sc[[paste('CD4',q$receiver[i],q$filter[i],sep='_')]]
   if(!q$target[i]%in%rownames(rz$avg))next
   for(sk in names(sc)[vapply(sc,function(z)z$lineage=='Mac',logical(1))]) {
     sz <- sc[[sk]];lg <- q$ligand[i]
     if(!lg%in%rownames(sz$avg))next
     keep <- rz$sm$eligible20&sz$sm$eligible20
     x <- sz$avg[lg,keep];y <- rz$avg[q$target[i],keep];gr <- rz$sm$group[keep]
     r <- if(sum(keep)>=6)partial_group(x,y,gr) else NA_real_
     lr <- if(sum(keep)>=6)vapply(seq_along(x),function(j)partial_group(x[-j],y[-j],gr[-j]),numeric(1)) else NA_real_
     pr[[length(pr)+1]] <- data.frame(receiver=q$receiver[i],filter=q$filter[i],sender=sz$scope,ligand=lg,target=q$target[i],canonical=q$canonical[i],n=sum(keep),n_control=sum(gr=='Saline'),n_case=sum(gr=='Silica'),rho_group_adjusted=r,loo_min=if(any(is.finite(lr)))min(lr,na.rm=TRUE) else NA_real_,loo_max=if(any(is.finite(lr)))max(lr,na.rm=TRUE) else NA_real_,loo_same_fraction=if(is.finite(r)&&r!=0&&any(is.finite(lr)))mean(sign(lr[is.finite(lr)])==sign(r)) else NA_real_,evidence='DESCRIPTIVE_NOT_CAUSAL')
   }
 }
}
save_rows(pr,'paired_ligand_target_correlations.tsv',c('receiver','filter','sender','ligand','target','n','rho_group_adjusted'))
cat('MATCHED_SAMPLE_CHECK_PASS\n')

