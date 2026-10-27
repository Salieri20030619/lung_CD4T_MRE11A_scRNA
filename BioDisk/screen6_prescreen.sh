#!/usr/bin/env bash
set -Eeuo pipefail
export PROJECT_ROOT="${PROJECT_ROOT:-/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA}"
mkdir -p "$PROJECT_ROOT/02_analysis" "$PROJECT_ROOT/logs" "$PROJECT_ROOT/tmp"
exec 9>"$PROJECT_ROOT/tmp/screen6_prescreen.lock"
if ! flock -n 9; then echo SCREEN6_PRESCREEN_ALREADY_RUNNING; exit 0; fi
trap 'rc=$?; echo "SCREEN6_PRESCREEN_FAIL EXIT_CODE=$rc"; exit "$rc"' ERR
TASK_DIR=$(mktemp -d "$PROJECT_ROOT/02_analysis/screen6_prescreen_$(date +%Y%m%d_%H%M%S)_XXXXXX")
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
cat > "$TASK_DIR/prescreen.py" <<'SCREEN6_PY'
from pathlib import Path
import os, sys, gzip, csv, io, tarfile, json, traceback
import numpy as np
import pandas as pd
from scipy import io as sio, sparse
import scanpy as sc
ROOT=Path(os.environ.get('PROJECT_ROOT','/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA'))
OUT=Path(sys.argv[1]); OUT.mkdir(parents=True,exist_ok=True)
GENES='CD3D CD3E CD3G TRAC TRBC1 TRBC2 CD247 CD4 CD8A CD8B CD8B1 IL7R LTB TCF7 LEF1 CCR7 SELL FOXP3 IL2RA CTLA4 ICOS NKG7 GNLY NCR1 KLRD1 GZMB CCL5 LYZ LYZ2 LST1 TYROBP EPCAM SFTPC COL1A1 MS4A1 CD79A MRE11 MRE11A ENSMUSG00000031928 ENSG00000020922'.split()
TARGET=['MRE11','MRE11A','ENSMUSG00000031928','ENSG00000020922']
TGENES=['CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2']
reports=[]; audits=[]; statuses=[]
def canon(s):
 s=str(s).strip().upper()
 return s.split('.')[0] if s.startswith('ENS') else s

def metadata(g):
 ans={}; current=None
 for f in (ROOT/'00_raw'/g).glob('*family.soft.gz'):
  with gzip.open(f,'rt') as h:
   for line in h:
    if line.startswith('^SAMPLE = '): current=line.strip().split(' = ',1)[1];ans[current]=[]
    elif current and line.startswith(('!Sample_title =','!Sample_characteristics_ch1 =','!Sample_source_name_ch1 =')):ans[current].append(line.strip())
 return ans

def emit(g,lib,cells,total,nfeat,m,features,integral,extra=None):
 n=len(cells); d=pd.DataFrame({'cell':list(map(str,cells)),'sample':lib,'matrix_total':total,'n_features':nfeat})
 for k,v in m.items():d[k]=v
 if extra is not None:
  for k in extra:d[k]=np.asarray(extra[k])
 if g=='GSE184854':
  d['sample']=[str(x).split('_')[0] for x in cells]
  d['sample']=lib+':'+d['sample']
 d['eligible']=~d.cell.str.contains('doublet|unassigned|undetermined',case=False,regex=True)
 def val(k):return np.asarray(m.get(k,np.zeros(n)))
 d['T_marker_n']=sum((val(k)>0).astype(int) for k in TGENES)
 d['T_candidate']=d.T_marker_n>=2
 if 'productive_TCR' in d:d['T_candidate']=d.T_candidate|d.productive_TCR
 d['CD4_supported']=d.T_candidate & (val('CD4')>0) & ((val('CD8A')+val('CD8B')+val('CD8B1'))==0)
 targets=[k for k in TARGET if k in m]
 d['target_value']=sum((val(k) for k in targets),np.zeros(n)) if targets else np.nan
 d['target_detected']=d.target_value>0
 if 'Cell_Subtype' in d:
  d['author_CD4']=d.Cell_Subtype.isin(['CD4 effector T cell','CD4 helper T cell','Regulatory T cell'])
 tiers=['all','T_candidate','CD4_supported']+(['author_CD4'] if 'author_CD4' in d else [])
 for sample,ds in d.groupby('sample',sort=False):
  for tier in tiers:
   z=ds[ds.eligible & (True if tier=='all' else ds[tier])];count=len(z)
   reports.append(dict(dataset=g,library=lib,sample=sample,tier=tier,n_cells=count,target_features=';'.join(targets),target_positive=int(z.target_detected.sum()) if targets else None,target_pct=100*z.target_detected.mean() if count and targets else None,target_matrix_sum=z.target_value.sum() if targets else None,median_features=z.n_features.median() if count else None,median_matrix_total=z.matrix_total.median() if count else None,integer_nonnegative_matrix=bool(integral)))
 d.to_csv(OUT/f'{g}_{lib}_cells.tsv.gz',sep='\t',index=False)
 audits.append(dict(dataset=g,library=lib,n_cells=n,n_features=features,integer_nonnegative_matrix=bool(integral),markers_present=list(m),target_features=targets))
 print('LIBRARY_DONE',g,lib,n,'target=',targets,flush=True)

def from_sparse(g,lib,X,cells,genes,extra=None):
 X=sparse.csr_matrix(X);assert X.shape==(len(cells),len(genes))
 if not np.isfinite(X.data).all():raise ValueError('nonfinite expression')
 integral=np.all(X.data>=0) and np.allclose(X.data,np.rint(X.data),atol=1e-6,rtol=0)
 names=np.array([canon(x) for x in genes]);m={k:np.asarray(X[:,names==k].sum(axis=1)).ravel() for k in GENES if np.any(names==k)}
 emit(g,lib,cells,np.asarray(X.sum(1)).ravel(),np.asarray((X>0).sum(1)).ravel(),m,len(genes),integral,extra)

def stream(g,lib,handle,orientation,sep=','):
 reader=csv.reader((s for s in handle if s.strip() and not s.lstrip().startswith('#')),delimiter=sep)
 header=next(reader);m={};integral=True;seen=0
 if orientation=='genes':
  cells=header[1:];total=np.zeros(len(cells));nfeat=np.zeros(len(cells),dtype=int)
  for row in reader:
   if len(row)!=len(header):raise ValueError(f'row width {len(row)} != {len(header)} at {seen}')
   a=np.asarray(row[1:],dtype=float)
   if not np.isfinite(a).all():raise ValueError('nonfinite expression')
   integral=integral and bool(np.all(a>=0) and np.allclose(a,np.rint(a),atol=1e-6,rtol=0))
   total+=a;nfeat+=a>0;seen+=1;k=canon(row[0])
   if k in GENES:m[k]=m.get(k,0)+a
  if seen==0:raise ValueError('empty matrix')
  emit(g,lib,cells,total,nfeat,m,seen,integral)
 else:
  genes=[canon(k) for k in header[1:]]
  if not set(genes).intersection(GENES):raise ValueError('No recognizable marker in cell-by-gene header; orientation needs review: '+str(header[:12]))
  idx={k:np.where(np.array(genes)==k)[0] for k in GENES if k in genes}
  vals={k:[] for k in idx};cells=[];total=[];nfeat=[]
  for row in reader:
   if len(row)!=len(header):raise ValueError('cell-by-gene row width mismatch')
   a=np.asarray(row[1:],dtype=float)
   if not np.isfinite(a).all():raise ValueError('nonfinite expression')
   integral=integral and bool(np.all(a>=0) and np.allclose(a,np.rint(a),atol=1e-6,rtol=0))
   cells.append(row[0]);total.append(a.sum());nfeat.append((a>0).sum())
   for k,ii in idx.items():vals[k].append(a[ii].sum())
  emit(g,lib,cells,total,nfeat,{k:np.array(v) for k,v in vals.items()},len(genes),integral)

def run(g):
 raw=ROOT/'00_raw'/g
 (OUT/f'{g}_sample_metadata.json').write_text(json.dumps(metadata(g),indent=2,ensure_ascii=False))
 if g=='GSE184854':
  files=sorted(raw.glob('*matrix_inflection*.txt.gz'));assert len(files)==2
  for f in files:
   with gzip.open(f,'rt') as h:stream(g,f.name.split('_')[0],h,'genes','\t')
 elif g=='GSE174725':
  with tarfile.open(raw/'GSE174725_RAW.tar') as t:
   files=[f for f in t if f.isfile() and f.name.endswith('.csv.gz')];assert len(files)==5
   for f in files:
    with t.extractfile(f) as b,gzip.GzipFile(fileobj=b) as z,io.TextIOWrapper(z) as h:stream(g,Path(f.name).name.split('_')[0],h,'genes')
 elif g=='GSE295463':
  files=sorted(raw.glob('*MolsPerCell.csv.gz'));assert len(files)==5
  for f in files:
   with gzip.open(f,'rt') as h:stream(g,f.name.split('_')[0],h,'cells')
 elif g=='GSE183682':
  mapping={'GSM5567501':'GSM5567509','GSM5567502':'GSM5567510','GSM5567503':'GSM5567511','GSM5567504':'GSM5567512'}
  files=sorted(raw.glob('*matrix.mtx.gz'));assert len(files)==4
  for f in files:
   prefix=f.name[:-len('matrix.mtx.gz')];lib=f.name.split('_')[0]
   cells=pd.read_csv(raw/(prefix+'barcodes.tsv.gz'),header=None,sep='\t')[0].astype(str).values
   genes=pd.read_csv(raw/(prefix+'features.tsv.gz'),header=None,sep='\t')[1].values
   with gzip.open(f,'rb') as h:X=sio.mmread(h).T.tocsr()
   extra={};tc=list(raw.glob(mapping[lib]+'*TCR*filtered_contig_annotations.csv.gz'))
   if tc:
    v=pd.read_csv(tc[0]);good=v.chain.isin(['TRA','TRB'])
    for k in ['productive','high_confidence']:good &= v[k].astype(str).str.lower().eq('true')
    extra['productive_TCR']=np.isin(cells,v.loc[good,'barcode'].values)
   from_sparse(g,lib,X,cells,genes,extra);del X
 elif g=='GSE159354':
  files=sorted(raw.glob('*filtered_feature_bc_matrix.h5'));assert len(files)==10
  for f in files:
   a=sc.read_10x_h5(f);from_sparse(g,f.name.split('_')[0],a.X,a.obs_names,a.var_names);del a
 elif g=='GSE296682':
  a=sc.read_h5ad(raw/'GSE296682_BAL_FINAL.h5ad',backed='r')
  try:
   names=np.array([canon(s) for s in a.var_names]);idx=np.where(np.isin(names,GENES))[0]
   mats=[];tot=[];nf=[];integral=True
   for start in range(0,a.n_obs,1000):
    X=sparse.csr_matrix(a.X[start:start+1000,:]);mats.append(X[:,idx].toarray());tot.extend(np.asarray(X.sum(1)).ravel());nf.extend(np.asarray((X>0).sum(1)).ravel())
    if not np.isfinite(X.data).all():raise ValueError('nonfinite expression')
    integral=integral and bool(np.all(X.data>=0) and np.allclose(X.data,np.rint(X.data),atol=1e-6,rtol=0))
   V=np.concatenate(mats);m={k:V[:,names[idx]==k].sum(1) for k in GENES if np.any(names[idx]==k)}
   extra={k:a.obs[k].astype(str).values for k in ['Cell_Subtype','Cell_Sub_subtype','disease_final','Pool'] if k in a.obs}
   if 'Sample' not in a.obs:raise ValueError('missing donor Sample')
   extra['sample']=a.obs['Sample'].astype(str).values
   emit(g,'BAL_FINAL',a.obs_names,tot,nf,m,a.n_vars,integral,extra)
  finally:a.file.close()

(OUT/'README.txt').write_text('Screening only, not final CD4 annotation. T_candidate: >=2 detected T-lineage genes, or productive high-confidence TRA/TRB where provided. CD4_supported: T_candidate + CD4 detected + no CD8A/B/B1 detection. These gates are depth-dependent; no detection is not proof of absence. No new QC or doublet algorithm was applied. Author CD4 is reported separately. Matrix values may be normalized; integer values alone do not establish raw UMI provenance. TAS-Seq values are not called UMI. Do not compare absolute expression across platforms. No p-values or dataset ranking by target significance. Sample metadata must be reviewed before donor-level inference. GSE184854 sample prefixes retained without unverified time mapping.\n')
for g in ['GSE184854','GSE183682','GSE296682','GSE295463','GSE174725','GSE159354']:
 print('DATASET_START='+g,flush=True)
 try:run(g);statuses.append({'dataset':g,'status':'PASS'})
 except Exception as e:
  traceback.print_exc();statuses.append({'dataset':g,'status':'FAIL','error':str(e)})
 pd.DataFrame(reports).to_csv(OUT/'target_summary.tsv',sep='\t',index=False)
 (OUT/'matrix_audit.json').write_text(json.dumps(audits,indent=2))
 (OUT/'status.json').write_text(json.dumps(statuses,indent=2))
with tarfile.open(ROOT/'tmp/screen6_prescreen_return.tar.gz','w:gz') as t:
 for p in sorted(OUT.iterdir()):
  if p.is_file():t.add(p,arcname=p.name)
passed=all(x['status']=='PASS' for x in statuses)
print('SCREEN6_PRESCREEN_PASS' if passed else 'SCREEN6_PRESCREEN_PARTIAL',flush=True)
sys.exit(0 if passed else 1)
SCREEN6_PY
conda run -n batch_bench python -m py_compile "$TASK_DIR/prescreen.py"
echo SYNTAX_CHECK=PASS
echo "OUTPUT=$TASK_DIR"
conda run --no-capture-output -n batch_bench python -u "$TASK_DIR/prescreen.py" "$TASK_DIR"
