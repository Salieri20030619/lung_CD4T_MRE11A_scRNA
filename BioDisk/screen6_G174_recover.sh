#!/usr/bin/env bash
set -Eeuo pipefail
export PROJECT_ROOT="${PROJECT_ROOT:-/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA}"
mkdir -p "$PROJECT_ROOT/tmp" "$PROJECT_ROOT/02_analysis" "$PROJECT_ROOT/logs"
exec 9>"$PROJECT_ROOT/tmp/screen6_G174_recover.lock"
if ! flock -n 9; then echo G174_RECOVERY_ALREADY_RUNNING; exit 0; fi
trap 'rc=$?; echo "G174_RECOVERY_FAIL EXIT_CODE=$rc"; exit "$rc"' ERR
TASK=$(mktemp -d "$PROJECT_ROOT/02_analysis/screen6_G174_recover_XXXXXX")
cat > "$TASK/recover.py" <<'RECOVERY_PY'
from pathlib import Path
import os,sys,csv,gzip,io,tarfile,json
import numpy as np
import pandas as pd
ROOT=Path(os.environ.get('PROJECT_ROOT','/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA'))
OUT=Path(sys.argv[1]);OUT.mkdir(parents=True,exist_ok=True)
MARKERS='CD3D CD3E CD3G TRAC TRBC1 TRBC2 CD247 CD4 CD8A CD8B IL7R LTB TCF7 LEF1 CCR7 SELL FOXP3 IL2RA CTLA4 ICOS NKG7 GNLY NCR1 KLRD1 GZMB CCL5 LYZ LST1 TYROBP EPCAM SFTPC COL1A1 MS4A1 CD79A MRE11 MRE11A'.split()
summary=[];audit=[]
with tarfile.open(ROOT/'00_raw/GSE174725/GSE174725_RAW.tar') as t:
 members=[m for m in t if m.isfile() and m.name.endswith('.csv.gz')]
 assert len(members)==5
 for member in members:
  lib=Path(member.name).name.split('_')[0]
  with t.extractfile(member) as b,gzip.GzipFile(fileobj=b) as z,io.TextIOWrapper(z) as h:
   reader=csv.reader(h);header=next(reader);cells=np.array(header[1:]);n=len(cells)
   assert len(set(cells))==n and all(cells)
   total=np.zeros(n);nf=np.zeros(n,dtype=int);markers={};blank=0;missing=0;genes=0;integer=True
   excluded=(cells=='GGCTTTCAGCCAGTAG-1') if lib=='GSM5324866' else np.zeros(n,dtype=bool)
   assert excluded.sum()==(1 if lib=='GSM5324866' else 0)
   keep=~excluded
   for line,row in enumerate(reader,2):
    if not row or not any(v.strip() for v in row):blank+=1;continue
    assert len(row)==len(header),f'{lib}: row width mismatch at {line}'
    gene=row[0].strip().upper();assert gene
    values=row[1:]
    empties=[j for j,v in enumerate(values) if not v.strip()]
    if empties:
     assert lib=='GSM5324866' and gene=='IFNGR2' and len(empties)==1 and cells[empties[0]]=='GGCTTTCAGCCAGTAG-1',f'Unexpected missing field: {lib}:{gene}'
     missing+=1
    # Exclude the entire audited cell before numeric conversion; no imputation.
    a=np.array([v for j,v in enumerate(values) if keep[j]],dtype=float)
    assert np.isfinite(a).all() and (a>=0).all(),f'Invalid expression: {lib}:{gene}'
    integer=integer and bool(np.allclose(a,np.rint(a),atol=1e-6,rtol=0))
    total[keep]+=a;nf[keep]+=a>0;genes+=1
    if gene in MARKERS:markers[gene]=markers.get(gene,0)+a
  assert missing==(1 if lib=='GSM5324866' else 0)
  d=pd.DataFrame({'cell':cells[keep],'sample':lib,'matrix_total':total[keep],'n_features':nf[keep]})
  for k,v in markers.items():d[k]=v
  def val(k):return np.asarray(markers.get(k,np.zeros(len(d))))
  d['T_marker_n']=sum((val(k)>0).astype(int) for k in ['CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2'])
  d['T_candidate']=d.T_marker_n>=2
  d['CD4_supported']=d.T_candidate & (val('CD4')>0) & ((val('CD8A')+val('CD8B'))==0)
  targets=[k for k in ['MRE11','MRE11A'] if k in markers];assert targets
  d['target_value']=sum(val(k) for k in targets);d['target_detected']=d.target_value>0
  for tier in ['all','T_candidate','CD4_supported']:
   s=d if tier=='all' else d[d[tier]]
   summary.append(dict(dataset='GSE174725',sample=lib,tier=tier,n_cells=len(s),target_positive=int(s.target_detected.sum()),target_pct=100*s.target_detected.mean() if len(s) else None,target_matrix_sum=s.target_value.sum(),median_matrix_total=s.matrix_total.median() if len(s) else None,median_features=s.n_features.median() if len(s) else None,integer_nonnegative_matrix=integer))
  d.to_csv(OUT/f'{lib}_cells.tsv.gz',sep='\t',index=False)
  record=dict(sample=lib,input_cells=n,retained_cells=len(d),excluded_cells=list(cells[excluded]),blank_rows_skipped=blank,missing_fields=missing,features=genes,markers_present=list(markers),integer_nonnegative_matrix=integer)
  audit.append(record);print('LIBRARY_DONE',json.dumps(record),flush=True)
pd.DataFrame(summary).to_csv(OUT/'GSE174725_target_summary.tsv',sep='\t',index=False)
(OUT/'recovery_audit.json').write_text(json.dumps(audit,indent=2))
(OUT/'README.txt').write_text('One cell with audited IFNGR2 missing expression excluded entirely; no missing expression imputed. Blank rows skipped. Marker gates are preliminary, not final CD4 annotation. No new QC/doublet detection applied. Positive gene detection alone does not establish T lineage; re-clustering required before expression inference.\n')
with tarfile.open(ROOT/'tmp/screen6_G174_recovered_return.tar.gz','w:gz') as t:
 for p in sorted(OUT.iterdir()):
  if p.is_file():t.add(p,arcname=p.name)
print('G174_PRESCREEN_RECOVERED_PASS',flush=True)
RECOVERY_PY
conda run -n batch_bench python -m py_compile "$TASK/recover.py"
echo SYNTAX_CHECK=PASS
conda run --no-capture-output -n batch_bench python -u "$TASK/recover.py" "$TASK"
