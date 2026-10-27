#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA
CONDA=/mnt/BioDisk/wangmingyuan/miniconda3/bin/conda
mkdir -p "$ROOT/logs" "$ROOT/tmp" "$ROOT/02_analysis/GSE174725"
exec 9>"$ROOT/tmp/GSE174725_suppC2.lock"
if ! flock -n 9; then echo SUPPC2_ALREADY_RUNNING; exit 0; fi
RUN=$(mktemp -d "$ROOT/02_analysis/GSE174725/suppC2_$(date +%Y%m%d_%H%M%S)_XXXXXX")
export RUN
trap 'rc=$?; echo "GSE174725_SUPPC2_FAIL EXIT_CODE=$rc OUTPUT=$RUN"; exit "$rc"' ERR
mkdir -p "$RUN/input" "$RUN/tables" "$RUN/figures"
SRC="$ROOT/02_analysis/screen2_Trefine_nuzVQu/GSE174725/T_refined.h5ad"
echo "OUTPUT=$RUN"
test -s "$SRC"
cat > "$RUN/export.py" <<'C2_SOURCE'
from pathlib import Path
import sys
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse,io
out=Path(sys.argv[1]);src=Path(sys.argv[2]);a=ad.read_h5ad(src)
assert a.n_obs==10286 and a.var_names.is_unique and a.obs_names.is_unique
assert 'MRE11' in a.var_names
cl=a.obs.F08.astype(str);a=a[cl.isin(['0','1','12','2','3'])].copy()
x=sparse.csr_matrix(a.layers['counts']);x.eliminate_zeros()
assert np.isfinite(x.data).all() and (x.data>=0).all() and np.allclose(x.data,np.round(x.data))
m=a.obs.copy();m['cell']=a.obs_names.astype(str);m['F08']=m.F08.astype(str)
m['subtype']=m.F08.map({'0':'MemCD4','1':'NaiveCD4','12':'Treg','2':'Mix2','3':'Mix3'})
m['core']=m.F08.isin(['0','1','12']);m['MRE11_UMI']=x[:,list(a.var_names).index('MRE11')].toarray().ravel();m['library_UMI']=np.asarray(x.sum(axis=1)).ravel()
m['group']=m['sample'].map({'GSM5324866':'Exposed','GSM5324867':'Exposed','GSM5324868':'Silicosis','GSM5324869':'Silicosis','GSM5324870':'Silicosis'})
assert m.group.notna().all() and (m.library_UMI>0).all()
old=pd.read_csv(out/'input/expected_samples.tsv',sep='\t');old=old[old.dataset=='GSE174725']
for scope in ['Core','Expanded']:
 z=m[m.core] if scope=='Core' else m
 for ss,v in z.groupby('sample'):
  e=old[(old.scope==scope)&(old['sample']==ss)].iloc[0]
  assert len(v)==e.n_cells and v.MRE11_UMI.sum()==e.target_counts and np.isclose(v.library_UMI.sum(),e.total_counts)
assert m.core.sum()==2958 and len(m)==5481
m.to_csv(out/'input/cells.tsv.gz',sep='\t',index=False)
pd.DataFrame({'symbol':a.var_names}).to_csv(out/'input/genes.tsv',sep='\t',index=False)
io.mmwrite(out/'input/counts.mtx',x.T)
(out/'input/provenance.txt').write_text('Object='+str(src)+'\nCounts=layers[counts]; full measured feature background; pre-fixed F08 Core=0,1,12; Expanded adds 2,3.\n')
print('ANALYSIS_LOGIC_PASS INPUT_COUNTS_AND_SCOPE',flush=True)
C2_SOURCE
cat > "$RUN/score.R" <<'C2_SOURCE'
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
 take<-m$sample==ss & if(scope=='Core')m$core else rep(TRUE,nrow(m))
 for(n in names(sigs)){
  z<-x[sigs[[n]],take,drop=FALSE];k<-k+1L
  components[[k]]<-data.frame(scope=scope,sample=ss,module=n,symbol=rownames(z),cells=ncol(z),positive=Matrix::rowSums(z>0),UMI=Matrix::rowSums(z))
 }
}
wr(do.call(rbind,components),'component_detection.tsv')
wr(AnnotationDbi::metadata(org.Hs.eg.db::org.Hs.eg.db),'annotation_metadata.tsv')
capture.output(sessionInfo(),file=file.path(out,'tables/sessionInfo.txt'))
cat('HUMAN_UCELL_PASS\n')
C2_SOURCE
cat > "$RUN/summarize.py" <<'C2_SOURCE'
from pathlib import Path
import sys
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path(sys.argv[1]);m=pd.read_csv(out/'input/cells.tsv.gz',sep='\t');s=pd.read_csv(out/'tables/cell_scores.tsv',sep='\t');mods=list(s.columns[1:]);m=m.merge(s,on='cell',validate='one_to_one')
rows=[]
for scope in ['Core','Expanded']:
 z=m[m.core] if scope=='Core' else m
 for subtype in ['All']+sorted(z.subtype.unique().tolist()):
  zz=z if subtype=='All' else z[z.subtype==subtype]
  for ss,v in zz.groupby('sample'):
   for mod in mods:
    rows.append(dict(scope=scope,subtype=subtype,sample=ss,group=v.group.iloc[0],module=mod,cells=len(v),score_mean=v[mod].mean(),score_median=v[mod].median(),mean_UMI=v.library_UMI.mean(),MRE11_pct=100*(v.MRE11_UMI>0).mean(),MRE11_CPM=1e6*v.MRE11_UMI.sum()/v.library_UMI.sum()))
a=pd.DataFrame(rows);a.to_csv(out/'tables/sample_scores.tsv',sep='\t',index=False)
rows=[]
for (scope,mod),v in a[a.subtype=='All'].groupby(['scope','module']):
 for rule in [20,1]:
  z=v[v.cells>=rule];e=z[z.group=='Exposed'];c=z[z.group=='Silicosis'];effect=c.score_mean.mean()-e.score_mean.mean()
  rows.append(dict(scope=scope,module=mod,min_cells=rule,n_exposed=len(e),n_silicosis=len(c),exposed_mean=e.score_mean.mean(),silicosis_mean=c.score_mean.mean(),difference=effect,direction='Higher' if effect>0 else 'Lower' if effect<0 else 'Equal',p_status='descriptive_small_sample'))
r=pd.DataFrame(rows);r.to_csv(out/'tables/group_directions.tsv',sep='\t',index=False)
# Composition-standardized score, common weights from eligible Core samples.
base=m[m.core];sizes=base.groupby('sample').size();eligible=sizes[sizes>=20].index
ct=pd.crosstab(base['sample'],base.subtype);frac=ct.div(ct.sum(axis=1),axis=0);weights=frac.loc[eligible].mean();weights.rename('weight').to_csv(out/'tables/standardization_weights.tsv',sep='\t')
std=[]
for (ss,mod),v in a[(a.scope=='Core')&(a.subtype!='All')].groupby(['sample','module']):
 vv=v.set_index('subtype').reindex(weights.index);complete=vv.score_mean.notna().all()
 std.append(dict(sample=ss,group=v.group.iloc[0],module=mod,cells=sizes[ss],eligible20=ss in eligible,score_std=(vv.score_mean*weights).sum() if complete else np.nan,status='complete' if complete else 'missing_subtype'))
std=pd.DataFrame(std);std.to_csv(out/'tables/core_standardized_scores.tsv',sep='\t',index=False)
samples=['GSM5324866','GSM5324867','GSM5324868','GSM5324869','GSM5324870'];short=['Exp1','Exp2','Sil1','Sil2','Sil3']
fig,axs=plt.subplots(2,len(mods),figsize=(3*len(mods),6),squeeze=False)
for j,scope in enumerate(['Core','Expanded']):
 for i,mod in enumerate(mods):
  ax=axs[j,i];v=a[(a.scope==scope)&(a.subtype=='All')&(a.module==mod)].set_index('sample').reindex(samples)
  for k,(_,rr) in enumerate(v.iterrows()):
   color='#31688e' if k<2 else '#c55635';ax.scatter(k,rr.score_mean,facecolors=color if rr.cells>=20 else 'none',edgecolors=color,s=42);ax.annotate('n='+str(int(rr.cells)),(k,rr.score_mean),xytext=(0,7),textcoords='offset points',ha='center',fontsize=7)
  ax.set(title=scope+' / '+mod,xticks=range(5),xticklabels=short,ylabel='UCell');ax.margins(y=.3);ax.spines[['top','right']].set_visible(False)
fig.tight_layout()
for ext in ['png','pdf']:fig.savefig(out/'figures'/('human_scores.'+ext),dpi=350,bbox_inches='tight')
plt.close(fig)
fig,ax=plt.subplots(figsize=(5,3.5));q=r[r.min_cells==20].pivot(index='module',columns='scope',values='difference').reindex(mods);vmax=max(np.abs(q.to_numpy()).max(),1e-8);im=ax.imshow(q.to_numpy(),cmap='RdBu_r',vmin=-vmax,vmax=vmax,aspect='auto');ax.set(xticks=range(2),xticklabels=q.columns,yticks=range(len(q)),yticklabels=q.index)
for i in range(len(q)):
 for j in range(2):ax.text(j,i,f'{q.iloc[i,j]:+.4f}',ha='center',va='center',color='black' if abs(q.iloc[i,j])<vmax*.6 else 'white')
fig.colorbar(im,ax=ax,label='Sil - Exp');fig.tight_layout()
for ext in ['png','pdf']:fig.savefig(out/'figures'/('scope_effects.'+ext),dpi=350,bbox_inches='tight')
plt.close(fig)
(out/'captions.txt').write_text('Human scores: x-axis independent samples Exp1/Exp2=exposed comparators GSM5324866/67, Sil1/Sil2/Sil3=silicosis GSM5324868/69/70. Y-axis mean per-cell UCell score within sample/scope. Filled points have >=20 cells; hollow points have <20 cells; n denotes cell number. Core=F08 clusters 0/1/12, Expanded adds mixed CD4/CD8 clusters 2/3. OxPhos=oxidative phosphorylation; OxStress=response to oxidative stress; Inflam=inflammasome-mediated signaling pathway; PyroResp=pyroptotic inflammatory response. Scope effects: rows modules, columns cell scope, color and printed value are equally weighted sample mean Silicosis minus Exposed for >=20-cell samples. Positive indicates higher score in silicosis. Descriptive differences, not significance. Standardized scores use common Core subtype weights from eligible samples; absent subtypes yield missing values. MemCD4=memory-like CD4, NaiveCD4=naive-like CD4, Treg=regulatory T, Mix2/Mix3=mixed CD4/CD8 clusters.\n')
(out/'methods.txt').write_text('Four GO programs fixed before human scoring, same GO definitions as mouse cohorts; human GOALL symbols derived from org.Hs.eg.db, not capitalization or claimed identical cross-species gene membership. Actual annotated/present/detected genes and annotation versions exported. MRE11 and MRE11A excluded. UCell 2.14.0 maxRank=1500, full-feature raw counts, signature restricted to present features exactly as main-cohort approach. Minimum 3 present features; coverage and detection must be interpreted. Full measured feature background retained. Core >=20 cells primary descriptive direction; Expanded and all-sample summaries sensitivity. No new differential-expression, group-significance or MRE11-correlation tests for this five-subject supplement. All subject points retained; groups weighted equally. Composition-standardization uses fixed Core subtype proportions. Scores describe transcriptional programs; no cross-species absolute-score comparison.\n')
print('GSE174725_SUPPC2_ANALYSIS_PASS',flush=True)
C2_SOURCE
cat > "$RUN/input/expected_samples.tsv" <<'C2_DATA'
dataset	scope	sample	group	day	n_cells	target_positive	target_pct	target_counts	total_counts	target_CPM	log2_CPM1	median_library	eligible20
GSE174725	Core	GSM5324866	Exposed	0	22	2	9.090909090909092	2.0	13064.0	153.09247	7.2676525	570.5	True
GSE174725	Core	GSM5324867	Exposed	0	151	18	11.920529801324504	19.0	722989.0	26.27979	4.7697606	4339.0	True
GSE174725	Core	GSM5324868	Silicosis	0	50	2	4.0	2.0	144338.0	13.856364	3.8930092	2877.5	True
GSE174725	Core	GSM5324869	Silicosis	0	5	1	20.0	1.0	10454.0	95.657166	6.594805	1658.0	False
GSE174725	Core	GSM5324870	Silicosis	0	2730	164	6.007326007326007	174.0	1.4249735e+07	12.210753	3.723641	5091.5	True
GSE174725	Expanded	GSM5324866	Exposed	0	27	2	7.4074074074074066	2.0	15821.0	126.41426	6.993383	566.0	True
GSE174725	Expanded	GSM5324867	Exposed	0	2639	209	7.9196665403561965	229.0	1.4097993e+07	16.243446	4.1079764	5155.0	True
GSE174725	Expanded	GSM5324868	Silicosis	0	50	2	4.0	2.0	144338.0	13.856364	3.8930092	2877.5	True
GSE174725	Expanded	GSM5324869	Silicosis	0	11	1	9.090909090909092	1.0	20290.0	49.285362	5.6520667	1658.0	False
GSE174725	Expanded	GSM5324870	Silicosis	0	2754	166	6.027596223674655	177.0	1.435935e+07	12.326464	3.736222	5091.0	True
GSE183682	Core	GSM5567501	NS	7	271	7	2.5830258302583027	8.0	753900.0	10.611486	3.5374808	2643.0	True
GSE183682	Core	GSM5567502	Silica	7	393	2	0.5089058524173028	2.0	508797.0	3.9308407	2.3018336	1160.0	True
GSE183682	Core	GSM5567503	NS	56	416	22	5.288461538461538	24.0	1.403141e+06	17.10448	4.178275	3205.5	True
GSE183682	Core	GSM5567504	Silica	56	333	9	2.7027027027027026	9.0	1.140086e+06	7.894141	3.1528554	3248.0	True
C2_DATA
"$CONDA" run --no-capture-output -n batch_bench python -m py_compile "$RUN/export.py" "$RUN/summarize.py"
"$CONDA" run --no-capture-output -n R_4.5.2 Rscript -e 'invisible(parse(file=file.path(Sys.getenv("RUN"),"score.R"))); req<-c("Matrix","UCell","BiocParallel","AnnotationDbi","org.Hs.eg.db"); missing<-req[!vapply(req,requireNamespace,logical(1),quietly=TRUE)]; if(length(missing))stop(paste("MISSING_PACKAGES",paste(missing,collapse=","))); cat("R_PARSE_PASS\n")'
echo SYNTAX_CHECK=PASS
"$CONDA" run --no-capture-output -n batch_bench python -u "$RUN/export.py" "$RUN" "$SRC"
"$CONDA" run --no-capture-output -n R_4.5.2 Rscript "$RUN/score.R" "$RUN"
"$CONDA" run --no-capture-output -n batch_bench python -u "$RUN/summarize.py" "$RUN"
printf 'GSE174725_SUPPC2_PASS\n' > "$RUN/STATUS.txt"
PACK="$ROOT/GSE174725_suppC2_return.tar.gz"
tar -C "$RUN" -czf "$PACK.part" tables figures input/cells.tsv.gz input/genes.tsv input/expected_samples.tsv input/provenance.txt export.py score.R summarize.py captions.txt methods.txt STATUS.txt
mv "$PACK.part" "$PACK"
echo "RETURN_PACKAGE=$PACK"
echo GSE174725_SUPPC2_PASS
