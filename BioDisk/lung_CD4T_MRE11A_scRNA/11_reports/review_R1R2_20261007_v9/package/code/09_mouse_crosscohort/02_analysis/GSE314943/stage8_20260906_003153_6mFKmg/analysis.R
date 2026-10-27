args<-commandArgs(TRUE);src<-args[1];out<-args[2]
req<-c('Matrix','UCell','AnnotationDbi','org.Mm.eg.db','GO.db','ggplot2','BiocParallel')
miss<-req[!vapply(req,requireNamespace,logical(1),quietly=TRUE)]
if(length(miss)){
 options(repos=c(CRAN='https://cloud.r-project.org'),timeout=600)
 if(!requireNamespace('BiocManager',quietly=TRUE))install.packages('BiocManager')
 BiocManager::install(miss,ask=FALSE,update=FALSE)
}
stopifnot(all(vapply(req,requireNamespace,logical(1),quietly=TRUE)))
suppressPackageStartupMessages({library(Matrix);library(UCell);library(AnnotationDbi);library(org.Mm.eg.db);library(GO.db);library(ggplot2)})
for(d in c('tables','figures','objects'))dir.create(file.path(out,d),showWarnings=FALSE)
wr<-function(x,n)write.table(x,file.path(out,'tables',n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
genes<-read.delim(file.path(out,'input/genes.tsv'),stringsAsFactors=FALSE)
meta<-read.delim(file.path(out,'input/cells.tsv'),stringsAsFactors=FALSE,check.names=FALSE)
x<-as(readMM(file.path(out,'input/counts.mtx')),'CsparseMatrix')
rownames(x)<-genes$gene_id;colnames(x)<-meta$cell
stopifnot(ncol(x)==2311,nrow(x)==nrow(genes),!anyDuplicated(rownames(x)),!anyDuplicated(colnames(x)))
samples<-c(paste0('Saline_',1:3),paste0('Silica_',1:3))
stopifnot(identical(as.integer(table(factor(meta$sample,levels=samples))),c(304L,424L,164L,553L,390L,476L)))
defs<-data.frame(module=c('OxPhos','MitoOrg','OxStress','Inflam','PyroResp','InflamRsp','DDR'),
 GO=c('GO:0006119','GO:0007005','GO:0006979','GO:0141084','GO:0070269','GO:0006954','GO:0006974'))
terms<-AnnotationDbi::select(GO.db,keys=defs$GO,keytype='GOID',columns='TERM')
defs$description<-terms$TERM[match(defs$GO,terms$GOID)]
# GOALL includes descendant annotations; export the actual annotated Ensembl members.
ga<-AnnotationDbi::select(org.Mm.eg.db,keys=defs$GO,keytype='GOALL',columns='ENSEMBL')
ga<-unique(ga[!is.na(ga$ENSEMBL),c('GOALL','ENSEMBL')])
sigs<-setNames(lapply(defs$GO,function(g)unique(ga$ENSEMBL[ga$GOALL==g])),defs$module)
sigs$PyroExec<-genes$gene_id[genes$symbol %in% c('Casp1','Casp4','Gsdmd','Gsdme')]
defs<-rbind(defs,data.frame(module='PyroExec',GO='Manual',description='Candidate pyroptosis execution machinery: Casp1, Casp4, Gsdmd, Gsdme; transcript abundance only'))
mre<-genes$gene_id[genes$gene_id=='ENSMUSG00000031928'|genes$symbol %in% c('Mre11a','Mre11')]
stopifnot('ENSMUSG00000031928'%in%mre)
coverage<-list();members<-list()
for(n in names(sigs)){
 full<-setdiff(unique(sigs[[n]]),mre);present<-intersect(full,rownames(x));sigs[[n]]<-present
 detected<-present[Matrix::rowSums(x[present,,drop=FALSE])>0]
 coverage[[n]]<-data.frame(module=n,annotated_n=length(full),present_n=length(present),detected_n=length(detected),coverage=length(present)/max(1,length(full)))
 members[[n]]<-data.frame(module=n,gene_id=full,present=full%in%present,detected=full%in%detected)
}
wr(defs,'module_definitions.tsv');wr(do.call(rbind,coverage),'module_coverage.tsv');wr(do.call(rbind,members),'module_members.tsv')
# Small/empty modules remain explicitly unavailable; do not manufacture scores.
usable<-names(sigs)[lengths(sigs)>=3];stopifnot(length(usable)>0)
set.seed(42)
score<-as.data.frame(UCell::ScoreSignatures_UCell(x,features=sigs[usable],maxRank=1500,BPPARAM=BiocParallel::SerialParam()))
stopifnot(all(meta$cell %in% rownames(score)))
score<-score[meta$cell,,drop=FALSE];colnames(score)<-sub('_UCell$','',colnames(score))
stopifnot(all(usable%in%names(score)),all(is.finite(as.matrix(score))))
wr(data.frame(cell=meta$cell,score),'cell_scores.tsv')
saveRDS(list(scores=score,definitions=defs,members=sigs,metadata=meta),file.path(out,'objects/CD4_scores.rds'))
# Score summaries and detection are aggregated by sample, not cell-level tests.
pb<-read.delim(file.path(src,'tables/PB_sample_metadata.tsv'),stringsAsFactors=FALSE)
rows<-list();q<-0L
for(ss in c('Core','Naive','ActCD4','Treg'))for(s in samples){
 take<-meta$sample==s & if(ss=='Core')rep(TRUE,nrow(meta)) else meta$Subtype==ss
 z<-pb[pb$set==ss & pb$sample==s,,drop=FALSE];stopifnot(nrow(z)==1,sum(take)==z$cells)
 for(mod in usable){q<-q+1L;rows[[q]]<-data.frame(set=ss,sample=s,group=z$group,cells=sum(take),module=mod,
 score_mean=mean(score[take,mod]),score_median=median(score[take,mod]),Mre_CPM=z$Mre_CPM,Mre_pct=z$Mre_pct,
 status=if(ss=='Treg')'exploratory_low_cells' else 'primary')}
}
a<-do.call(rbind,rows);wr(a,'sample_scores.tsv')
# Component-level detection in every sample, retained for execution-module interpretation.
component<-list();q<-0L
for(mod in usable)for(s in samples){
 z<-x[sigs[[mod]],meta$sample==s,drop=FALSE]
 q<-q+1L;component[[q]]<-data.frame(module=mod,sample=s,gene_id=rownames(z),symbol=genes$symbol[match(rownames(z),genes$gene_id)],
 pct=100*Matrix::rowMeans(z>0),mean_umi=Matrix::rowMeans(z))
}
wr(do.call(rbind,component),'component_detection.tsv')
# Exact label tests: 20 group allocations; 720 pooled correlation permutations;
# 36 within-group permutations for correlation after regressing group from ranks.
perms<-function(v){if(length(v)==1)return(matrix(v,nrow=1));do.call(rbind,lapply(seq_along(v),function(i)cbind(v[i],perms(v[-i]))))}
allperm<-perms(1:6);p3<-perms(1:3);within<-do.call(rbind,lapply(seq_len(nrow(p3)),function(i)cbind(matrix(p3[i,],nrow=6,ncol=3,byrow=TRUE),p3+3)))
alloc<-combn(1:6,3);grp<-factor(rep(c('Saline','Silica'),each=3))
safecor<-function(v,w){if(sd(v)<1e-12||sd(w)<1e-12)return(NA_real_);cor(v,w)}
cor_exact<-function(v,w,adjust=FALSE){
 rv<-rank(v);rw<-rank(w)
 if(adjust){rv<-residuals(lm(rv~grp));rw<-residuals(lm(rw~grp))}
 r<-safecor(rv,rw);if(!is.finite(r))return(c(rho=NA_real_,PValue=NA_real_))
 pp<-if(adjust)within else allperm
 null<-apply(pp,1,function(ix)safecor(rv,rw[ix]))
 c(rho=r,PValue=mean(abs(null)>=abs(r)-1e-12))
}
tests<-list();cors<-list();q<-0L;qq<-0L
for(ss in unique(a$set))for(mod in usable){
 d<-a[a$set==ss & a$module==mod,];d<-d[match(samples,d$sample),];v<-d$score_mean
 delta<-mean(v[4:6])-mean(v[1:3]);null<-apply(alloc,2,function(i)mean(v[i])-mean(v[-i]))
 q<-q+1L;tests[[q]]<-data.frame(set=ss,module=mod,Saline=mean(v[1:3]),Silica=mean(v[4:6]),difference=delta,PValue=mean(abs(null)>=abs(delta)-1e-12),status=d$status[1])
 for(metric in c('Mre_CPM','Mre_pct'))for(ad in c(FALSE,TRUE)){
 z<-cor_exact(d[[metric]],v,ad);qq<-qq+1L;cors[[qq]]<-data.frame(set=ss,module=mod,metric=metric,method=if(ad)'GroupAdj' else 'Pooled',rho=z['rho'],PValue=z['PValue'],status=d$status[1])
 }
}
t<-do.call(rbind,tests);t$FDR<-NA_real_
for(st in unique(t$status)){i<-t$status==st;t$FDR[i]<-p.adjust(t$PValue[i],'BH')}
r<-do.call(rbind,cors);r$FDR<-NA_real_
for(st in unique(r$status))for(method in unique(r$method)){i<-r$status==st & r$method==method;r$FDR[i]<-p.adjust(r$PValue[i],'BH')}
wr(t,'score_group_tests.tsv');wr(r,'Mre11a_correlations.tsv')
# Within-set/module scaling is only for the heatmap; tests use raw UCell scores.
figsave<-function(p,n,caption,w=9,h=7){for(ext in c('png','pdf'))ggsave(file.path(out,'figures',paste0(n,'.',ext)),p,width=w,height=h,dpi=350,bg='white');writeLines(c(caption,paste(defs$module,defs$description,sep=' = ')),file.path(out,'figures',paste0(n,'_caption.txt')))}
a$Z<-ave(a$score_mean,interaction(a$set,a$module),FUN=function(v){if(sd(v)==0)rep(0,length(v))else as.numeric(scale(v))})
a$label<-c('Sal1','Sal2','Sal3','Si1','Si2','Si3')[match(a$sample,samples)]
a$label<-factor(a$label,levels=c('Sal1','Sal2','Sal3','Si1','Si2','Si3'))
p<-ggplot(a,aes(label,module,fill=Z))+geom_tile()+facet_wrap(~set,ncol=2)+scale_fill_gradient2(low='#0072B2',mid='white',high='#D55E00')+labs(x=NULL,y=NULL,fill='Z score')+theme_minimal(base_size=11)
figsave(p,'scores_heatmap','Sample-mean UCell scores, standardized separately within each cell set and module across six samples. Sal1-3=Saline_1-3; Si1-3=Silica_1-3. Core=all selected CD4; Naive=naive-like CD4; ActCD4=activated CD4; Treg=regulatory T cells (exploratory). PyroExec measures transcript abundance, not enzyme cleavage or cell death.')
d<-r[r$metric=='Mre_CPM',];p<-ggplot(d,aes(set,module,fill=rho))+geom_tile()+facet_wrap(~method)+scale_fill_gradient2(low='#0072B2',mid='white',high='#D55E00',limits=c(-1,1),na.value='grey85')+labs(x=NULL,y=NULL,fill='Rho')+theme_minimal(base_size=11)
figsave(p,'Mre11a_corr','Sample-level associations between Mre11a pseudobulk CPM and mean UCell score. Pooled=Spearman correlation over six samples; GroupAdj=correlation of rank residuals after adjusting for treatment group. NA (grey) denotes no variation. Each correlation uses six samples, not individual cells. Exact P values and BH FDR are in Mre11a_correlations.tsv; color alone does not indicate significance.')
for(ss in unique(a$set)){
 z<-a[a$set==ss,];p<-ggplot(z,aes(group,score_mean,color=group))+geom_point(position=position_jitter(width=.08,seed=42),size=2)+facet_wrap(~module,scales='free_y',ncol=4)+scale_color_manual(values=c(Saline='#0072B2',Silica='#D55E00'))+labs(x=NULL,y='UCell')+theme_classic(base_size=11)+theme(legend.position='none')
 figsave(p,paste0(ss,'_scores'),'Each point is a sample mean; each group contains three samples. Module panels have independent y scales. Scores use raw RNA counts with UCell maxRank=1500; Mre11a/Mre11 is excluded from all module memberships to avoid circular correlation.',10,6)
}
writeLines(c(paste('Input:',src),
 'GO modules use installed GO.db names and org.Mm.eg.db GOALL descendant annotations. Full annotated membership, matrix coverage, detected membership and component detection are retained.',
 'PyroExec is a predefined candidate machinery panel (Casp1,Casp4,Gsdmd,Gsdme), not a validated activity assay. GO:0070269 uses the installed term name, which may be pyroptotic inflammatory response.',
 'UCell uses all 29809 matrix features as its ranking background, raw counts, maxRank=1500, serial execution. Module membership includes genes present in the matrix even if unexpressed. Modules with fewer than three available features are omitted and shown in module_coverage.tsv.',
 'Target gene Mre11a/Mre11 is excluded from every module. No cell selection by target detection, no imputation, no smoothing.',
 'Group score comparisons use exact 3+3 sample permutations (20 allocations). Correlations use exact 720 permutations pooled and 36 within-group permutations for rank residual correlations adjusted for group.',
 'BH: group tests pooled across modules and Core/Naive/ActCD4; Treg separately. Correlations adjusted separately for Pooled and GroupAdj, pooling both target metrics and primary sets/modules; Treg separately.',
 'Core associations can reflect subtype composition; subtype results are provided separately. Treg remains exploratory due to low cell counts. GroupAdj still has only six samples.',
 'Transcript scores are not direct measurements of mitochondrial damage, inflammasome assembly, gasdermin cleavage or Mre11a causal activity.'),file.path(out,'tables/methods.txt'))
capture.output(sessionInfo(),file=file.path(out,'tables/sessionInfo.txt'))
message('GSE314943_STAGE8_ANALYSIS_PASS')
