from pathlib import Path
import os,json,gzip,csv,tarfile,hashlib,re,subprocess,collections,time
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import numpy as np,pandas as pd
from scipy import sparse,io
from threadpoolctl import threadpool_limits
import anndata as ad
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/GSE174725_reaudit_20261005_v1';R=O/'results';raw=P/'00_raw/GSE174725';ext=P/'01_raw_data/GSE174725/extracted';old=P/'02_analysis/screen6_G174_recover_8E2LLp'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for c in iter(lambda:f.read(8*1024*1024),b''):h.update(c)
 return h.hexdigest()
def wr(d,n):pd.DataFrame(d).to_csv(R/n,sep='\t',index=False,na_rep='NA')
print('START',time.strftime('%FT%T'),flush=True)
features=pd.read_csv(raw/'GSE174725_features.tsv.gz',sep='\t',header=None);bc=pd.read_csv(raw/'GSE174725_barcodes.tsv.gz',header=None)[0].astype(str);suf=bc.str.rsplit('-',n=1).str[-1];base=bc.str.rsplit('-',n=1).str[0];assert bc.is_unique
rscript='/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript'
expr='x<-read.delim(gzfile("'+str(raw/'GSE174725_features.tsv.gz')+'"),header=FALSE);write.table(data.frame(symbol=x[[2]],Rname=make.names(x[[2]],unique=FALSE)),"'+str(R/'pooled_R_make_names.tsv')+'",sep="\\t",quote=FALSE,row.names=FALSE)'
subprocess.run([rscript,'-e',expr],check=True,capture_output=True,text=True);rnames=pd.read_csv(R/'pooled_R_make_names.tsv',sep='\t').Rname.astype(str)
with threadpool_limits(limits=2):
 with gzip.open(raw/'GSE174725_matrix.mtx.gz','rb') as f:pool=io.mmread(f).tocsr().astype('int32')
assert pool.shape==(len(features),len(bc));assert pool.min()>=0
print('POOLED_LOADED',pool.shape,pool.nnz,flush=True)
ids=features[0].astype(str);symbols=features[1].astype(str);species={'features':len(ids),'ENSG':int(ids.str.startswith('ENSG').sum()),'ENSMUSG':int(ids.str.startswith('ENSMUSG').sum()),'unique_ids':bool(ids.is_unique),'duplicate_symbols':int(symbols.duplicated().sum()),'integer_nonnegative':True,'matrix_shape':list(pool.shape),'nnz':pool.nnz,'suffix_counts':suf.value_counts().to_dict()};(R/'pooled_structure_species.json').write_text(json.dumps(species,indent=2))
symbolcounts=collections.Counter(symbols);literal={g:i for i,g in enumerate(symbols) if symbolcounts[g]==1}
hashes=[]
for f in sorted(raw.iterdir()):hashes.append(dict(path=str(f),bytes=f.stat().st_size,sha256=sha(f)))
with tarfile.open(raw/'GSE174725_RAW.tar') as t:
 archive=[]
 for m in t.getmembers():
  if not m.isfile():continue
  b=t.extractfile(m).read();p=ext/Path(m.name).name;assert hashlib.sha256(b).hexdigest()==sha(p);archive.append(dict(member=m.name,bytes=m.size,extracted_identical=True,sha256=sha(p)))
wr(archive,'archive_integrity.tsv');(O/'provenance/raw_sha256.json').write_text(json.dumps(hashes,indent=2))
scope=pd.read_csv(P/'02_analysis/screen2_CD4_target_wz1lfx/GSE174725/scope_membership.tsv.gz',sep='\t');print('SCOPE_COLUMNS',list(scope),flush=True)
core=ad.read_h5ad(P/'02_analysis/screen2_CD4_target_wz1lfx/GSE174725/CD4_core.h5ad');print('CORE',core.shape,list(core.obs),flush=True)
assert core.obs_names.is_unique
SM={'GSM5324866':('EXP1','Exposed','1'),'GSM5324867':('EXP2','Exposed','2'),'GSM5324868':('SIL1','Silicosis','3'),'GSM5324869':('SIL2','Silicosis','5'),'GSM5324870':('SIL3','Silicosis','4')}
MARK='CD3D CD3E CD3G TRAC TRBC1 TRBC2 CD4 CD8A CD8B LYZ LST1 TYROBP FCER1G AIF1 CTSS C1QA C1QB C1QC MARCO MRC1 CD68 CSF1R MS4A1 CD79A EPCAM'.split()
structure=[];comparisons=[];overlaps=[];cellresults=[];availability=[];coreaudit=[];genemap=[];allbarcodes=[]
for path in sorted(ext.glob('*.csv.gz')):
 gsm=path.name.split('_')[0];name,group,expected=SM[gsm];print('CSV_START',gsm,flush=True)
 with gzip.open(path,'rt') as f:header=next(csv.reader(f))
 d=pd.read_csv(path,index_col=0,dtype={c:'float32' for c in header[1:]});blank=d.index.isna()&d.isna().all(axis=1);assert int(blank.sum())==(1 if gsm=='GSM5324868' else 0);d=d.loc[~blank];genes=d.index.astype(str);cells=d.columns.astype(str);assert cells.is_unique
 missing=np.argwhere(~np.isfinite(d.to_numpy()));assert len(missing)==(1 if gsm=='GSM5324866' else 0)
 if len(missing):
  i,j=missing[0];assert genes[i]=='IFNGR2' and cells[j]=='GGCTTTCAGCCAGTAG-1';d=d.drop(columns=[cells[j]])
 genes=d.index.astype(str);cells=d.columns.astype(str);v=d.to_numpy();assert np.isfinite(v).all() and (v>=0).all() and (v==np.rint(v)).all()
 total=v.sum(axis=0,dtype='float64');nf=(v>0).sum(axis=0);mt=np.asarray(genes.str.upper().str.startswith(('MT-','MT.')));mtcount=v[mt].sum(axis=0,dtype='float64');qc=(nf>=200)&(total>=300)&(mtcount/np.maximum(total,1)<.2)
 counts=collections.Counter(genes)
 rcount=collections.Counter(rnames);rn={g:i for i,g in enumerate(rnames) if rcount[g]==1}
 gi=[];pi=[];modes=[]
 for i,g in enumerate(genes):
  if counts[g]!=1:continue
  if g in literal:gi.append(i);pi.append(literal[g]);modes.append('literal')
  elif g in rn:gi.append(i);pi.append(rn[g]);modes.append('R_make_names')
 # Explicitly reject any non-injective symbol transformations.
 pcount=collections.Counter(pi);keep=[k for k,p in enumerate(pi) if pcount[p]==1];gi=np.array(gi)[keep];pi=np.array(pi)[keep];modes=np.array(modes)[keep]
 cb=cells.str.rsplit('-',n=1).str[0];assert cb.is_unique
 for i,p,mode in zip(gi,pi,modes):genemap.append(dict(sample=gsm,csv_row=int(i),csv_symbol=genes[i],pooled_row=int(p),ensembl=ids[p],pooled_symbol=symbols[p],mode=mode))
 dates=[g for g in genes if re.fullmatch(r'\d{1,2}-(Mar|Sep|Dec|Oct|Jan|Feb|Apr|May|Jun|Jul|Aug|Nov)',g)]
 structure.append(dict(sample=gsm,label=name,group=group,input_cells=len(header)-1,retained_cells=len(cells),rows=len(genes),unique_genes=len(set(genes)),duplicate_label_rows=int(genes.duplicated().sum()),date_label_rows=len(dates),literal_mapped=int(sum(modes=='literal')),Rname_mapped=int(sum(modes=='R_make_names')),unmapped_rows=len(genes)-len(gi),missing_cells_excluded=len(missing),total_counts=float(total.sum()),integer_nonnegative=True))
 wr(pd.DataFrame({'row':np.arange(len(genes)),'gene':genes})[lambda q:~q.row.isin(gi)],gsm+'_unmapped_gene_rows.tsv')
 for suffix in sorted(suf.unique()):
  bj=np.where(suf==suffix)[0];lookup={b:i for b,i in zip(base.iloc[bj],bj)};ci=np.array([j for j,b in enumerate(cb) if b in lookup],dtype=int);pj=np.array([lookup[cb[j]] for j in ci],dtype=int);overlaps.append(dict(sample=gsm,suffix=suffix,barcode_overlap=len(ci),csv_retained=len(cells),pooled_cells=len(bj),expected=suffix==expected))
  if not len(ci):continue
  # Limit peak array footprint; evaluate every unambiguous mapped gene, every overlapping cell.
  neq=np.zeros(len(ci),dtype='int64');cgt=np.zeros(len(ci),dtype='int64');clt=np.zeros(len(ci),dtype='int64');union=np.zeros(len(ci),dtype='int64');s1=np.zeros(len(ci),dtype='float64');s2=s1.copy();l1=s1.copy()
  for start in range(0,len(gi),512):
   a=v[np.ix_(gi[start:start+512],ci)];b=pool[pi[start:start+512]][:,pj].toarray();delta=a-b;neq+=(delta!=0).sum(axis=0);cgt+=(delta>0).sum(axis=0);clt+=(delta<0).sum(axis=0);union+=((a>0)|(b>0)).sum(axis=0);s1+=a.sum(axis=0);s2+=b.sum(axis=0);l1+=abs(delta).sum(axis=0)
  comparisons.append(dict(sample=gsm,suffix=suffix,expected=suffix==expected,overlap_cells=len(ci),shared_genes=len(gi),exact_vectors=int((neq==0).sum()),all_csv_le_pool=int((cgt==0).sum()),all_csv_ge_pool=int((clt==0).sum()),mismatch_entries=int(neq.sum()),nonzero_union=int(union.sum()),csv_shared_sum=float(s1.sum()),pooled_shared_sum=float(s2.sum()),absolute_difference=float(l1.sum()),median_cell_nonzero_agreement=float(np.median(1-neq/np.maximum(union,1)))))
  if suffix==expected:
   cellresults.append(pd.DataFrame(dict(sample=gsm,barcode=cells[ci],pooled_barcode=bc.iloc[pj].values,mismatched_genes=neq,csv_gt_pooled=cgt,csv_lt_pooled=clt,csv_shared_sum=s1,pooled_shared_sum=s2,absolute_difference=l1)))
 z={g:v[np.asarray(genes.str.upper()==g)].sum(axis=0) for g in MARK};T=sum(z[g]>0 for g in 'CD3D CD3E CD3G TRAC TRBC1 TRBC2'.split())>=2;my=sum(z[g]>0 for g in 'LYZ LST1 TYROBP FCER1G AIF1 CTSS'.split())>=2;mac=sum(z[g]>0 for g in 'C1QA C1QB C1QC MARCO MRC1 CD68 CSF1R'.split())>=1;compete=T|(z['MS4A1']>0)|(z['CD79A']>0)|(z['EPCAM']>0);conservative=my&mac&~compete
 oldd=pd.read_csv(old/(gsm+'_cells.tsv.gz'),sep='\t').set_index('cell').loc[cells];assert np.array_equal(total,oldd.matrix_total) and np.array_equal(nf,oldd.n_features)
 c=core[core.obs['sample']==gsm];bar=c.obs.barcode.astype(str);assert set(bar)<=set(cells);oldcore=c.layers['counts'].tocsr();diff=0;tested=0
 # Existing objects collapse upper-case duplicate labels by SUM; verify every measured feature on the frozen Core cells.
 cup=genes.str.upper();clook={b:i for i,b in enumerate(cells)};cidx=np.array([clook[b] for b in bar]);gidx={}
 for jj,gg in enumerate(cup):gidx.setdefault(gg,[]).append(jj)
 oldcore=oldcore.tocsc()
 for j,g in enumerate(c.var_names):
  if g not in gidx:continue
  arr=v[np.ix_(gidx[g],cidx)].sum(axis=0);ref=oldcore[:,j].toarray().ravel();diff+=int((arr!=ref).sum());tested+=len(arr)
 coreaudit.append(dict(sample=gsm,core_cells=c.n_obs,features=c.n_vars,compared_entries=tested,mismatch_entries=diff,barcodes_all_in_own_CSV=True));assert diff==0
 ss=scope[scope['sample']==gsm];core_n=int(ss.scope_core.sum());expanded_n=int(ss.scope_expanded.sum());assert core_n==c.n_obs
 availability.append(dict(sample=gsm,label=name,group=group,retained_csv_cells=len(cells),QC_cells=int(qc.sum()),T_marker_supported=int((T&qc).sum()),direct_CD4_supported=int((T&(z['CD4']>0)&(z['CD8A']==0)&(z['CD8B']==0)&qc).sum()),broad_myeloid=int((my&qc).sum()),macrophage_marker_supported=int((conservative&qc).sum()),myeloid_T_overlap=int((my&T&qc).sum()),existing_Core_CD4=core_n,existing_Expanded_CD4=expanded_n,paired_Core20_and_markerMac20=bool(core_n>=20 and (conservative&qc).sum()>=20)))
 md=pd.DataFrame(dict(sample=gsm,barcode=cells,total_counts=total,n_features=nf,mt_pct=100*mtcount/np.maximum(total,1),QC=qc,T_supported=T,myeloid_supported=my,macrophage_marker_supported=conservative));md.to_csv(R/(gsm+'_cell_feasibility_SERVER_ONLY.tsv.gz'),sep='\t',index=False)
 allbarcodes.extend([(gsm,b) for b in cb]);wr(structure,'CSV_structure_summary.tsv');wr(comparisons,'count_vector_comparison_all_suffixes.tsv');wr(availability,'paired_lineage_feasibility.tsv');wr(coreaudit,'existing_Core_count_integrity.tsv');wr(overlaps,'barcode_overlap_reaudit.tsv');wr(genemap,'unambiguous_gene_mapping.tsv.gz');pd.concat(cellresults).to_csv(R/'expected_suffix_cell_integrity_SERVER_ONLY.tsv.gz',sep='\t',index=False);print('CSV_DONE',gsm,json.dumps(availability[-1]),'COUNTS',json.dumps([x for x in comparisons if x['sample']==gsm and x['expected']]),flush=True)
 del d,v
wr(overlaps,'barcode_overlap_reaudit.tsv');wr(genemap,'unambiguous_gene_mapping.tsv.gz');pd.concat(cellresults).to_csv(R/'expected_suffix_cell_integrity_SERVER_ONLY.tsv.gz',sep='\t',index=False)
bcdf=pd.DataFrame(allbarcodes,columns=['sample','base_barcode']);reuse=bcdf.groupby('base_barcode')['sample'].nunique();shared=reuse[reuse>1];wr(bcdf[bcdf.base_barcode.isin(shared.index)],'cross_library_reused_barcodes_SERVER_ONLY.tsv.gz')
result={'status':'COMPLETE','samples':5,'namespace_cell_ids_unique':not bcdf.duplicated(['sample','base_barcode']).any(),'shared_barcodes_across_libraries':len(shared),'all_existing_Core_count_comparisons_exact':all(x['mismatch_entries']==0 for x in coreaudit),'archive_files_identical':len(archive),'new_disease_models':0,'new_communication_runs':0};(O/'provenance/validation_complete.json').write_text(json.dumps(result,indent=2,default=lambda x:bool(x)));print(json.dumps(result,indent=2,default=lambda x:bool(x)),flush=True)
