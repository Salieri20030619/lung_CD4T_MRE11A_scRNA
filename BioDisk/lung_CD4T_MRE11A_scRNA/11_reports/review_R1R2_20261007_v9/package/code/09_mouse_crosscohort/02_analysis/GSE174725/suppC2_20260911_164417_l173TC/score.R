a<-commandArgs(TRUE);out<-a[1]
req<-c('Matrix','UCell','BiocParallel','AnnotationDbi','org.Hs.eg.db')
missing<-req[!vapply(req,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))stop(paste('MISSING_PACKAGES',paste(missing,collapse=',')))
if(as.character(utils::packageVersion('UCell'))!='2.14.0')stop('UCell 2.14.0 required for main-cohort scoring consistency')
suppressPackageStartupMessages(library(Matrix))
x<-as(Matrix::readMM(file.path(out,'input/counts.mtx')),'CsparseMatrix')
g<-read.delim(file.path(out,'input/genes.tsv'),stringsAsFactors=FALSE)
m<-read.delim(gzfile(file.path(out,'input/cells.tsv.gz')),stringsAsFactors=FALSE,check.names=FALSE)
rownames(x)<-g$symbol;colnames(x)<-m$cell
stopifnot(ncol(x)==5481,!anyDuplicated(rownames(x)),!anyDuplicated(colnames(x)))
defs<-data.frame(module=c('OxPhos','OxStress','Inflam','PyroResp'),GO=c('GO:0006119','GO:0006979','GO:0141084','GO:0070269'),description=c('oxidative phosphorylation','response to oxidative stress','inflammasome-mediated signaling pathway','pyroptotic inflammatory response'))
ga<-AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,keys=defs$GO,keytype='GOALL',columns='SYMBOL')
ga<-unique(ga[!is.na(ga$SYMBOL),c('GOALL','SYMBOL')])
sigs<-setNames(lapply(defs$GO,function(go)setdiff(unique(ga$SYMBOL[ga$GOALL==go]),c('MRE11','MRE11A'))),defs$module)
wr<-function(z,n)write.table(z,file.path(out,'tables',n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
coverage<-list();members<-list()
for(n in names(sigs)){
 full<-sigs[[n]];present<-intersect(full,rownames(x));detected<-present[Matrix::rowSums(x[present,,drop=FALSE])>0]
 coverage[[n]]<-data.frame(module=n,annotated_n=length(full),present_n=length(present),detected_n=length(detected),coverage=length(present)/max(1,length(full)),usable=length(present)>=3)
 members[[n]]<-data.frame(module=n,symbol=full,present=full%in%present,detected=full%in%detected)
 sigs[[n]]<-present
}
wr(defs,'module_definitions.tsv');wr(do.call(rbind,coverage),'module_coverage.tsv');wr(do.call(rbind,members),'module_members.tsv')
usable<-names(sigs)[lengths(sigs)>=3];stopifnot(length(usable)>0)
set.seed(42)
s<-as.data.frame(UCell::ScoreSignatures_UCell(x,features=sigs[usable],maxRank=1500,BPPARAM=BiocParallel::SerialParam()))
stopifnot(all(m$cell%in%rownames(s)));s<-s[m$cell,,drop=FALSE];names(s)<-sub('_UCell$','',names(s))
stopifnot(all(usable%in%names(s)),all(is.finite(as.matrix(s))))
wr(data.frame(cell=m$cell,s),'cell_scores.tsv')
components<-list();k<-0L
for(scope in c('Core','Expanded'))for(ss in unique(m$sample)){
 take<-m$sample==ss & if(scope=="Core") {flag<-tolower(as.character(m$core));stopifnot(all(flag %in% c("true","false")));flag=="true"} else rep(TRUE,nrow(m))
 for(n in names(sigs)){
  z<-x[sigs[[n]],take,drop=FALSE];k<-k+1L
  components[[k]]<-data.frame(scope=scope,sample=ss,module=n,symbol=rownames(z),cells=ncol(z),positive=Matrix::rowSums(z>0),UMI=Matrix::rowSums(z))
 }
}
wr(do.call(rbind,components),'component_detection.tsv')
wr(AnnotationDbi::metadata(org.Hs.eg.db::org.Hs.eg.db),'annotation_metadata.tsv')
capture.output(sessionInfo(),file=file.path(out,'tables/sessionInfo.txt'))
cat('HUMAN_UCELL_PASS\n')
