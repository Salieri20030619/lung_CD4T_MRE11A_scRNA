#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA
OUT=$(mktemp -d "$ROOT/02_analysis/GSE314943/stage12_$(date +%Y%m%d_%H%M%S)_XXXXXX")
trap 'rc=$?; if ((rc!=0)); then echo "GSE314943_STAGE12_FAIL EXIT_CODE=$rc"; fi' EXIT
cat > "$OUT/export.py" <<'PYTHON'
import sys
from pathlib import Path
import numpy as np,pandas as pd,scanpy as sc
from scipy.io import mmwrite
from scipy.sparse import vstack,csr_matrix
root=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/02_analysis/GSE314943')
out=Path(sys.argv[1]);(out/'input').mkdir()
ma=sc.read_h5ad(root/'stage11_20260906_100314_pALwAB/objects/macrophages_core.h5ad')
t=sc.read_h5ad(root/'stage5_20260905_222457_a2dApT/objects/CD4_core.h5ad')
assert ma.n_obs==1796 and t.n_obs==2311
assert ma.var_names.equals(t.var_names) and not ma.obs_names.isin(t.obs_names).any()
assert ma.var.symbol.equals(t.var.symbol)
frames=[]
for a,lineage,field in [(ma,'Mac','MacType'),(t,'CD4','Subtype')]:
 frames.append(pd.DataFrame({'cell':a.obs_names,'sample':a.obs['sample'].values,'group':a.obs['group'].values,'lineage':lineage,'subtype':a.obs[field].astype(str).values}))
m=pd.concat(frames,ignore_index=True);assert m.cell.is_unique
x=vstack([ma.layers['counts'],t.layers['counts']],format='csr')
# Sum raw counts for duplicate gene symbols, then normalize in R.
symbols=ma.var['symbol'].astype('string').fillna('').to_numpy(dtype=str)
valid=np.array([bool(s.strip()) for s in symbols]);unique=pd.Index(sorted(set(symbols[valid])))
idx=unique.get_indexer(symbols[valid]);proj=csr_matrix((np.ones(valid.sum()),(np.flatnonzero(valid),idx)),shape=(len(symbols),len(unique)))
x=x@proj
mmwrite(out/'input/counts.mtx',x.T)
pd.DataFrame({'symbol':unique}).to_csv(out/'input/genes.tsv',sep='\t',index=False)
m.to_csv(out/'input/cells.tsv',sep='\t',index=False)
pd.DataFrame({'gene_id':ma.var_names,'symbol':symbols,'retained':valid}).to_csv(out/'input/gene_mapping.tsv',sep='\t',index=False)
print('EXPORT_PASS',flush=True)
PYTHON
cat > "$OUT/analysis.R" <<'RSCRIPT'
out<-commandArgs(TRUE)[1]
options(repos=c(CRAN='https://cloud.r-project.org'),timeout=900)
if(!requireNamespace('CellChat',quietly=TRUE)){
 if(!requireNamespace('BiocManager',quietly=TRUE))install.packages('BiocManager')
 bioc<-c('BiocNeighbors','ComplexHeatmap');bioc<-bioc[!vapply(bioc,requireNamespace,logical(1),quietly=TRUE)]
 if(length(bioc))BiocManager::install(bioc,ask=FALSE,update=FALSE)
 if(!requireNamespace('remotes',quietly=TRUE))install.packages('remotes')
 remotes::install_github('jinworks/CellChat',upgrade='never',dependencies=NA)
}
req<-c('CellChat','Matrix','ggplot2','future')
stopifnot(all(vapply(req,requireNamespace,logical(1),quietly=TRUE)))
suppressPackageStartupMessages({library(CellChat);library(Matrix);library(ggplot2)})
future::plan('sequential');set.seed(42)
for(d in c('tables','figures','objects'))dir.create(file.path(out,d),showWarnings=FALSE)
wr<-function(t,n)write.table(t,file.path(out,'tables',n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
x<-as(readMM(file.path(out,'input/counts.mtx')),'CsparseMatrix')
g<-read.delim(file.path(out,'input/genes.tsv'));m<-read.delim(file.path(out,'input/cells.tsv'),stringsAsFactors=FALSE)
rownames(x)<-g$symbol;colnames(x)<-m$cell;rownames(m)<-m$cell
stopifnot(ncol(x)==4107,!anyDuplicated(rownames(x)),all(colSums(x)>0))
y<-x%*%Diagonal(x=1e4/Matrix::colSums(x));y@x<-log1p(y@x)
colnames(y)<-m$cell;rownames(y)<-g$symbol
utils::data('CellChatDB.mouse',package='CellChat',envir=environment())
db<-CellChatDB.mouse
keep<-db$interaction$annotation %in% c('Secreted Signaling','ECM-Receptor','Cell-Cell Contact')
db$interaction<-db$interaction[keep,,drop=FALSE];stopifnot(nrow(db$interaction)>100)
saveRDS(db,file.path(out,'objects/CellChatDB_used.rds'))
wr(db$interaction,'database_interactions.tsv')
wr(as.data.frame(table(m$subtype,m$sample)),'subtype_sample_counts.tsv')
# Six broad sample models and two pooled broad/refined models per condition.
configs<-data.frame(name=c(paste0('Sample_',unique(m$sample)),'Saline_Broad','Silica_Broad','Saline_Fine','Silica_Fine'),
 level=c(rep('Sample',6),'Broad','Broad','Fine','Fine'),stringsAsFactors=FALSE)
allout<-list();statuses<-list()
for(i in seq_len(nrow(configs))){
 cfg<-configs[i,];message('CELLCHAT_START=',cfg$name)
 if(cfg$level=='Sample'){
  sn<-sub('^Sample_','',cfg$name);take<-m$sample==sn;condition<-unique(m$group[take])
 }else{sn<-'Pooled';condition<-strsplit(cfg$name,'_')[[1]][1];take<-m$group==condition}
 md<-m[take,,drop=FALSE];md$label<-if(cfg$level=='Fine')md$subtype else md$lineage
 cnt<-table(md$label);eligible<-names(cnt)[cnt>=20];md<-md[md$label%in%eligible,,drop=FALSE]
 md$label<-factor(md$label);meta<-data.frame(labels=md$label,samples=md$sample,row.names=md$cell)
 if(length(eligible)<2)stop('Insufficient cell groups: ',cfg$name)
 cc<-createCellChat(object=y[,md$cell,drop=FALSE],meta=meta,group.by='labels');cc@DB<-db
 cc<-subsetData(cc)
 extra<-if('do.fast'%in%names(formals(identifyOverExpressedGenes)))list(do.fast=FALSE)else list()
 cc<-do.call(identifyOverExpressedGenes,c(list(object=cc),extra))
 cc<-identifyOverExpressedInteractions(cc)
 cc<-computeCommunProb(cc,type='truncatedMean',trim=.1,raw.use=TRUE,population.size=FALSE,nboot=100,seed.use=42,distance.use=FALSE)
 cc<-filterCommunication(cc,min.cells=20)
 cc<-computeCommunProbPathway(cc);cc<-aggregateNet(cc)
 saveRDS(cc,file.path(out,'objects',paste0(cfg$name,'.rds')))
 ar<-cc@net$prob;pv<-cc@net$pval
 if(length(ar)==0)stop('Empty communication array: ',cfg$name)
 dn<-dimnames(ar);stopifnot(length(dn)==3)
 d<-expand.grid(source=dn[[1]],target=dn[[2]],interaction=dn[[3]],KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE)
 d$prob<-as.vector(ar);d$pvalue<-as.vector(pv)
 macros<-if(cfg$level=='Fine')c('AMlike','C1qMac','Folr2Mac','Spp1Mac')else 'Mac'
 ttypes<-if(cfg$level=='Fine')c('Naive','ActCD4','Treg')else 'CD4'
 d<-d[(d$source%in%macros & d$target%in%ttypes)|(d$source%in%ttypes & d$target%in%macros),,drop=FALSE]
 ann<-db$interaction[match(d$interaction,rownames(db$interaction)),,drop=FALSE]
 d$ligand<-ann$ligand;d$receptor<-ann$receptor;d$pathway<-ann$pathway_name
 d$FDR<-p.adjust(d$pvalue,'BH');d$model<-cfg$name;d$level<-cfg$level;d$sample<-sn;d$condition<-condition
 wr(d,paste0(cfg$name,'_LR.tsv'));allout[[cfg$name]]<-d
 statuses[[cfg$name]]<-data.frame(model=cfg$name,cells=nrow(md),groups=length(eligible),tested_cross=nrow(d),positive=sum(d$prob>0),nominal_P05=sum(d$prob>0 & d$pvalue<.05),FDR05=sum(d$prob>0 & d$FDR<.05))
 if(cfg$level=='Fine')wr(data.frame(group=names(cnt),cells=as.integer(cnt),included=names(cnt)%in%eligible),paste0(cfg$name,'_eligibility.tsv'))
 rm(cc);gc();message('CELLCHAT_PASS=',cfg$name)
}
wr(do.call(rbind,statuses),'CellChat_summary.tsv')
z<-do.call(rbind,allout);wr(z,'all_cross_LR.tsv')
# Replicate support is counted within each condition, without treating cells as replicates.
ss<-z[z$level=='Sample',];ids<-unique(ss[,c('source','target','interaction')]);rows<-list()
for(i in seq_len(nrow(ids))){
 q<-ids[i,];d<-ss[ss$source==q$source & ss$target==q$target & ss$interaction==q$interaction,]
 row<-q
 for(gr in c('Saline','Silica')){
  v<-d[d$condition==gr,];row[[paste0(gr,'_tested')]]<-nrow(v)
  row[[paste0(gr,'_positive')]]<-sum(v$prob>0)
  row[[paste0(gr,'_P05')]]<-sum(v$prob>0 & v$pvalue<.05)
  row[[paste0(gr,'_FDR05')]]<-sum(v$prob>0 & v$FDR<.05)
 }
 rows[[i]]<-row
}
if(length(rows))wr(do.call(rbind,rows),'replicate_support.tsv')
# Descriptive pooled-condition bubble plots, with complete labels in captions.
for(lev in c('Broad','Fine')){
 d<-z[z$level==lev & z$prob>0 & z$pvalue<.05,]
 if(!nrow(d)){writeLines('No positive nominal P<0.05 cross-lineage interaction.',file.path(out,'figures',paste0(lev,'_no_plot.txt')));next}
 strength<-tapply(d$prob,d$interaction,max);chosen<-names(sort(strength,decreasing=TRUE))[seq_len(min(20,length(strength)))]
 d<-d[d$interaction%in%chosen,];labs<-data.frame(interaction=chosen,label=sprintf('L%02d',seq_along(chosen)))
 d$label<-factor(labs$label[match(d$interaction,labs$interaction)],levels=rev(labs$label));d$edge<-paste(d$source,d$target,sep='>')
 labs$ligand<-db$interaction$ligand[match(labs$interaction,rownames(db$interaction))];labs$receptor<-db$interaction$receptor[match(labs$interaction,rownames(db$interaction))]
 wr(labs,paste0(lev,'_plot_labels.tsv'))
 # Edge axis codes keep all labels <=10 characters.
 edges<-unique(d$edge);ecode<-setNames(sprintf('E%02d',seq_along(edges)),edges);d$edgecode<-unname(ecode[d$edge])
 wr(data.frame(label=unname(ecode),edge=names(ecode)),paste0(lev,'_edge_labels.tsv'))
 p<-ggplot(d,aes(edgecode,label,size=prob,color=condition))+geom_point(position=position_dodge(width=.55))+
  scale_color_manual(values=c(Saline='#0072B2',Silica='#D55E00'))+labs(x='Cell pairs',y=NULL,size='CC prob')+theme_classic(base_size=11)
 for(ext in c('png','pdf'))ggsave(file.path(out,'figures',paste0(lev,'_bubble.',ext)),p,width=max(8,length(edges)*.35),height=8,dpi=350,bg='white')
 writeLines(c('CellChat pooled-condition descriptive screen; top20 interactions ranked by maximum positive probability among nominal bootstrap P<0.05 interactions. CC prob=CellChat communication probability, not a between-condition effect or P value. Different model scales and pooled sample composition must be considered. Full tables include BH FDR.',
 paste(labs$label,labs$interaction,labs$ligand,labs$receptor,sep=' | '),paste(unname(ecode),names(ecode),sep=' = '),
 'AMlike=alveolar-like macrophages; C1qMac=C1q-positive macrophages; Folr2Mac=Folr2-positive macrophages; Spp1Mac=Spp1-positive macrophages; Naive=naive-like CD4; ActCD4=activated CD4; Treg=regulatory T cells. Subtype comparisons with insufficient per-sample counts remain descriptive.'),file.path(out,'figures',paste0(lev,'_caption.txt')))
}
writeLines(c('Inputs: 1796 core macrophages and2311 core CD4 from Stage11/Stage5. Gene-symbol duplicates summed at raw count level; log1p CP10k normalization. No Harmony-corrected expression used.',
 'CellChat mouse database from installed package, saved exactly in objects/CellChatDB_used.rds. Secreted, ECM-receptor and contact interactions; multimeric receptor/ligand complexes and cofactors handled by CellChat.',
 'Six sample-level Broad models (Mac/CD4), two pooled-condition Broad models and two pooled-condition Fine models. Each group requires20 cells per fitted model. Fine pooled eligibility does not establish three biological replicates.',
 'truncatedMean trim0.1, raw.use TRUE, population.size FALSE,100 bootstrap permutations, seed42; no spatial distance. Overexpressed interactions selected separately per model; missing sample/model interaction means not tested, not measured absence.',
 'BH correction per model across extracted cross-lineage LR tests. Zero empirical bootstrap P values reflect finite100 permutations. No formal between-condition differential communication test is claimed.',
 'replicate_support distinguishes tested count, positive probability, nominal P and FDR support. Do not compare mean probabilities across separately scaled sample models or interpret sample bootstrap P values as group differences.',
 'Spp1Mac controls9,7,0; C1qMac control15 in one sample; Folr2Mac and Treg also have low per-sample counts. Fine results are pooled descriptive candidates. Primary replicate assessment is Broad Mac/CD4.',
 'Communication is inferred from RNA expression and curated interactions; it does not establish physical interaction or Mre11a-mediated signaling.'),file.path(out,'tables/methods.txt'))
capture.output(sessionInfo(),file=file.path(out,'tables/sessionInfo.txt'))
message('GSE314943_STAGE12_ANALYSIS_PASS')
RSCRIPT
conda run -n batch_bench python -m py_compile "$OUT/export.py"
conda run -n R_4.5.2 Rscript -e 'invisible(parse(file=commandArgs(TRUE)[1]));cat("R_PARSE_PASS\n")' "$OUT/analysis.R" > "$OUT/syntax.log" 2>&1
echo 'SYNTAX_CHECK=PASS'
echo "OUTPUT=$OUT"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
conda run --no-capture-output -n batch_bench python -u "$OUT/export.py" "$OUT"
conda run --no-capture-output -n R_4.5.2 Rscript "$OUT/analysis.R" "$OUT"
tar -C "$OUT" -czf "$ROOT/tmp/GSE314943_stage12_return.tar.gz.part" tables figures analysis.R export.py syntax.log
mv "$ROOT/tmp/GSE314943_stage12_return.tar.gz.part" "$ROOT/tmp/GSE314943_stage12_return.tar.gz"
echo 'GSE314943_STAGE12_PASS'
