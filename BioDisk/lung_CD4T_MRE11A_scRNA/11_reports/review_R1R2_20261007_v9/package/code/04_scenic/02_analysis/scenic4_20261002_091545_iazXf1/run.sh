#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA
PREV1="$ROOT/02_analysis/mechanism_step1_20260927_152218_s4V5eF"
SCENIC1="$ROOT/02_analysis/scenic1_20260929_103513_xRoEnr"
PY="$ROOT/software/pyscenic_0.12.1_py310/bin/python"
export ROOT PREV1 SCENIC1
export WORKERS="${WORKERS:-4}" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
bash -n "$0"
mkdir -p "$ROOT/logs" "$ROOT/tmp" "$ROOT/02_analysis"
exec 9>"$ROOT/tmp/lung_scenic4.lock"
flock -n 9 || { echo SCENIC4_ALREADY_RUNNING; exit 1; }
OUT=$(mktemp -d "$ROOT/02_analysis/scenic4_$(date +%Y%m%d_%H%M%S)_XXXXXX")
export OUT
trap 'rc=$?; if ((rc != 0)); then printf "LUNG_SCENIC4_FAIL EXIT_CODE=%s OUTPUT=%s\n" "$rc" "$OUT"; fi' EXIT
mkdir -p "$OUT/tables" "$OUT/objects" "$OUT/input" "$OUT/figures"
cp -- "$0" "$OUT/run.sh"
printf 'OUTPUT=%s\n' "$OUT"
[[ -x "$PY" ]]
cat > "$OUT/export.py" <<'PY_EXPORT'
import os,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy import sparse
import anndata as ad
root=Path(os.environ['ROOT']);p1=Path(os.environ['PREV1']);out=Path(os.environ['OUT'])
source=root/'02_analysis/screen2_CD4_target_wz1lfx/GSE183682/CD4_core.h5ad'
a=ad.read_h5ad(source);assert a.n_obs==1413 and a.obs_names.is_unique and 'counts' in a.layers
x=sparse.csr_matrix(a.layers['counts'],dtype=float);x.sum_duplicates();x.eliminate_zeros()
assert np.isfinite(x.data).all() and (x.data>=0).all() and np.allclose(x.data,np.rint(x.data))
m=pd.read_csv(p1/'objects/GSE183682_cell_scores.tsv.gz',sep='\t');assert m.cell.is_unique and set(m.cell)==set(a.obs_names)
idx=pd.Index(a.obs_names).get_indexer(m.cell);x=x[idx];assert np.allclose(np.asarray(x.sum(axis=1)).ravel(),m.library)
assert m['sample'].nunique()==4 and set(m.day)=={7,56} and set(m.group)=={'NS','Silica'}
assert m[['sample','group','day','subtype','library','bg']].notna().all().all() and (m.library>0).all()
design={'GSM5567501':('NS',7),'GSM5567502':('Silica',7),'GSM5567503':('NS',56),'GSM5567504':('Silica',56)}
assert set(m['sample'])==set(design)
assert all((r.group,int(r.day))==design[r.sample] for r in m.itertuples())
assert np.array_equal(a.obs['sample'].astype(str).to_numpy()[idx],m['sample'].astype(str).to_numpy())
original=a.var['symbol'].astype('string').fillna('').astype(str).to_numpy() if 'symbol' in a.var else np.array(['']*a.n_vars)
ids=pd.Index(a.var_names.astype(str)).str.replace(r'\.\d+$','',regex=True)
alias={'MRE11A':'MRE11','ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','PARK2':'PRKN'}
genes=[]
for sy,gi in zip(original,ids):
 if not sy.strip() or sy.upper().startswith('ENS'):sy='' if gi.upper().startswith('ENS') else gi
 sy=sy.strip().upper();genes.append(alias.get(sy,sy))
genes=np.asarray(genes);valid=(genes!='')&~np.char.startswith(genes,'ENS');unique=sorted(set(genes[valid]))
pr=sparse.csr_matrix((np.ones(valid.sum()),(np.flatnonzero(valid),pd.Index(unique).get_indexer(genes[valid]))),shape=(len(genes),len(unique)))
x=(x@pr).tocsr();assert len(unique)>5000
sparse.save_npz(out/'input/validation_counts.npz',x)
pd.DataFrame({'gene':unique}).to_csv(out/'input/genes.tsv',sep='\t',index=False)
m.to_csv(out/'input/cells.tsv.gz',sep='\t',index=False)
pd.DataFrame({'feature':ids,'symbol':original,'canonical':genes,'mapped':valid}).to_csv(out/'tables/symbol_mapping.tsv.gz',sep='\t',index=False)
sm=m.groupby(['sample','group','day','subtype'],observed=True).size().reset_index(name='cells');sm.to_csv(out/'tables/input_design.tsv',sep='\t',index=False)
(out/'input_manifest.json').write_text(json.dumps(dict(source=str(source),metadata=str(p1/'objects/GSE183682_cell_scores.tsv.gz'),cells=len(m),genes=len(unique)),indent=2))
print('ANALYSIS_LOGIC_PASS EXTERNAL_COUNTS_FROZEN_METADATA',flush=True)

PY_EXPORT
cat > "$OUT/validate.py" <<'PY_VALIDATE'
import os,json,tarfile,hashlib
from pathlib import Path
import numpy as np,pandas as pd
from scipy import sparse
import loompy
from pyscenic.aucell import aucell
from pyscenic.cli.utils import load_signatures
from ctxcore.genesig import GeneSignature
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(os.environ['ROOT']);OUT=Path(os.environ['OUT']);S1=Path(os.environ['SCENIC1']);P1=Path(os.environ['PREV1'])
ALIAS={'MRE11A':'MRE11','ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','PARK2':'PRKN'}
def canon(s):return ALIAS.get(str(s).upper(),str(s).upper())
def rd(p):return pd.read_csv(p,sep='\t')
def wr(x,n):pd.DataFrame(x).to_csv(OUT/'tables'/n,sep='\t',index=False)
def pairs(table,features):
 rows=[]
 for keys,d in table.groupby(['filter','scope','seed'],observed=True):
  for day in [7,56]:
   v=d[d.day==day];ctrl=v[v.group=='NS'];case=v[v.group=='Silica']
   if len(ctrl)!=1 or len(case)!=1:continue
   for f in features:
    a=float(ctrl[f].iloc[0]);b=float(case[f].iloc[0])
    rows.append(dict(filter=keys[0],scope=keys[1],seed=keys[2],day=day,feature=f,control=a,silica=b,delta=b-a,control_cells=int(ctrl.cells.iloc[0]),silica_cells=int(case.cells.iloc[0]),evidence='one_library_per_condition_time_no_inferential_P'))
 return pd.DataFrame(rows)
def main():
 x=sparse.load_npz(OUT/'input/validation_counts.npz');g=rd(OUT/'input/genes.tsv').gene.to_numpy();m=rd(OUT/'input/cells.tsv.gz')
 mods=rd(P1/'tables/module_definition.tsv');excluded=set(mods.gene.map(canon))|set('IL1R2 PENK TTC39C FOXP1 TNFRSF9 LAG3 MRE11'.split())
 effects=[];coverage=[];activities=[];phen=[];det=[]
 modulefeatures=['score_OXPHOS','score_MitoISR','score_MitoROS','score_Mitophagy','score_Inflam','score_PyroExec','expr_IL1R2']
 for arm in ['All','BG1']:
  take=np.ones(len(m),bool) if arm=='All' else (m.bg<=.01).to_numpy();vm=m.loc[take].reset_index(drop=True);vx=x[take]
  frozen=rd(S1/'Mouse'/arm/'cells.tsv')
  with loompy.connect(str(S1/'Mouse'/arm/'score.loom'),mode='r') as ds:
   dc=np.asarray(ds.ca.CellID).astype(str);dg=np.array([canon(z) for z in ds.ra.Gene]);dv=ds[:,:].T
  assert len(set(dg))==len(dg) and set(dc)==set(frozen.cell)
  pos=pd.Index(dc).get_indexer(frozen.cell);dv=dv[pos];dc=frozen.cell.to_numpy()
  common=np.array(sorted(set(dg)&set(g)));assert len(common)>0
  dj=pd.Index(dg).get_indexer(common);vj=pd.Index(g).get_indexer(common)
  vv=vx[:,vj].multiply(1e4/vm.library.to_numpy()[:,None]).toarray();vv=np.log1p(vv)
  expression=pd.DataFrame(np.vstack([dv[:,dj],vv]),columns=common,index=['D::'+z for z in dc]+['V::'+z for z in vm.cell])
  assert expression.index.is_unique and np.isfinite(expression.to_numpy()).all()
  for seed in [17,29,43]:
   regs=load_signatures(str(S1/'Mouse'/arm/f'seed{seed}/regulons.csv'));sig=[]
   for tf in ['IKZF2','STAT1']:
    match=[r for r in regs if canon(r.transcription_factor)==tf];assert len(match)==1
    r=match[0];original={canon(z):w for z,w in r.gene2weight.items()};retained={z:w for z,w in original.items() if z not in excluded|{tf}};kept={z:w for z,w in retained.items() if z in set(common)}
    coverage.append(dict(filter=arm,seed=seed,tf=tf,original_genes=len(original),DropOverlap_genes=len(retained),shared_genes=len(kept),shared_fraction=len(kept)/len(retained) if retained else np.nan,ranking_genes=len(common),discovery_ranking_genes=len(dg)))
    assert kept,'Empty external regulon: '+tf
    sig.append(GeneSignature(name=tf,gene2weight=kept))
   scores=aucell(expression,sig,seed=seed,num_workers=int(os.environ.get('WORKERS','4')),noweights=True)
   assert scores.index.is_unique and set(scores.index)==set(expression.index)
   assert np.isfinite(scores.to_numpy()).all(), 'Nonfinite AUCell scores'
   scores.to_csv(OUT/'objects'/f'{arm}_seed{seed}_shared_auc.csv.gz',compression='gzip')
   for cohort,meta,prefix in [('Discovery',frozen,'D::'),('External',vm,'V::')]:
    aa=scores.loc[[prefix+z for z in meta.cell]]
    for scope in ['Core']+sorted(meta.subtype.unique()):
     sel=np.ones(len(meta),bool) if scope=='Core' else meta.subtype.eq(scope).to_numpy()
     for ss,sub in meta.loc[sel].groupby('sample',observed=True):
      if len(sub)<20:continue
      ii=meta.index.get_indexer(sub.index);mu=aa.iloc[ii].mean()
      activities.append(dict(cohort=cohort,filter=arm,scope=scope,seed=seed,sample=ss,group=sub.group.iloc[0],day=int(sub.day.iloc[0]) if 'day' in sub else -1,cells=len(sub),IKZF2=mu['IKZF2'],STAT1=mu['STAT1']))
   print('TRANSFER_DONE='+arm+' seed='+str(seed),flush=True)
  # Phenotype comparisons use frozen scores; keep them independent of regulon definition.
  for scope in ['Core']+sorted(vm.subtype.unique()):
   z=vm if scope=='Core' else vm[vm.subtype==scope]
   for ss,sub in z.groupby('sample',observed=True):
    if len(sub)<20:continue
    row=dict(filter=arm,scope=scope,seed=0,sample=ss,group=sub.group.iloc[0],day=int(sub.day.iloc[0]),cells=len(sub))
    for f in modulefeatures:
     if f=='expr_IL1R2':row[f]=np.mean(sub['count_IL1R2']/sub.library*1e4)
     else:row[f]=sub[f].mean()
    phen.append(row)
    for gene in ['IKZF2','STAT1','IL1R2','IL18R1','IL18RAP','CXCR6','IL2RB','IL2RG']:
     jj=np.flatnonzero(g==gene);ix=vm.index.get_indexer(sub.index)
     val=vx[ix][:,jj[0]].toarray().ravel() if len(jj) else np.full(len(sub),np.nan)
     det.append(dict(**row,gene=gene,mean_CP10k=np.mean(val/sub.library.to_numpy()*1e4),pct=100*np.mean(val>0) if len(jj) else np.nan))
 a=pd.DataFrame(activities);wr(a,'frozen_regulon_sample_activity.tsv');wr(coverage,'transfer_coverage.tsv')
 e=pairs(a[a.cohort=='External'],['IKZF2','STAT1']);wr(e,'external_time_matched_TF_deltas.tsv')
 ph=pd.DataFrame(phen);wr(ph,'external_phenotype_sample.tsv');wr(pairs(ph,modulefeatures),'external_time_matched_phenotype_deltas.tsv');wr(det,'external_gene_detection.tsv')
 for (arm,scope,seed),z in a[a.cohort=='Discovery'].groupby(['filter','scope','seed'],observed=True):
  for tf in ['IKZF2','STAT1']:
   c=z[z.group=='Saline'][tf];s=z[z.group=='Silica'][tf]
   if len(c) and len(s):effects.append(dict(filter=arm,scope=scope,seed=seed,tf=tf,control_mean=c.mean(),silica_mean=s.mean(),delta=s.mean()-c.mean()))
 wr(effects,'discovery_shared_universe_deltas.tsv')
 # Transfer audit: compare discovery directions after common-universe rescoring with external day-specific directions.
 dd=pd.DataFrame(effects);validation=e.merge(dd.rename(columns={'tf':'feature','delta':'discovery_delta'})[['filter','scope','seed','feature','discovery_delta']],on=['filter','scope','seed','feature'],how='left')
 validation['direction_matches_discovery']=np.where(validation.discovery_delta.notna(),np.sign(validation.delta)==np.sign(validation.discovery_delta),np.nan)
 wr(validation,'direction_replication.tsv')
 fig,axes=plt.subplots(1,2,figsize=(8,4),sharey=False)
 for ax,tf in zip(axes,['IKZF2','STAT1']):
  z=e[(e.scope=='ActCD4')&(e.feature==tf)]
  for arm,col,off in [('All','#176B9B',-.08),('BG1','#C75327',.08)]:
   zz=z[z['filter']==arm]
   for seed,ss in zz.groupby('seed',observed=True):
    xx=ss.day.map({7:0,56:1}).to_numpy()+off+(seed-29)*.001
    ax.plot(xx,ss.delta,marker='o',ms=4,lw=.8,c=col,alpha=.85,label=arm if seed==17 else None)
  ax.axhline(0,c='gray',lw=.8,ls='--');ax.set_xticks([0,1],['D7','D56']);ax.set_title(tf);ax.set_ylabel('AUC diff.');ax.legend(frameon=False)
 fig.tight_layout()
 for ext in ['png','pdf']:fig.savefig(OUT/'figures'/f'ExtTF.{ext}',dpi=300)
 plt.close(fig)
 (OUT/'figures/captions.txt').write_text('ExtTF: external GSE183682 ActCD4, Silica minus NS per time point. D7/D56: days7/56. All/BG1: frozen filtering arms. Three lines per arm are frozen regulons from random seeds17/29/43, not biological replicates. Each time-condition has one library; no inferential P. AUC diff.: AUCell mean difference. Same common measured gene universe used to rescore discovery and external cells.\n')
 (OUT/'methods.txt').write_text('External transfer validation, not de novo network inference. GSE183682 four libraries: NS/Silica at days7/56. Frozen IKZF2 and STAT1 regulons from GSE314943 seeds17/29/43, All/BG1. Remove previous candidates, TF itself and frozen module genes (DropOverlap). Canonical aliases fixed as previous pipeline. Intersection of discovery ranking genes and external measured features; no external detection, disease or P-based gene selection. Rescore discovery and external on exactly same measured universe using unweighted AUCell; report shared target coverage without a new cutoff. Original full-library CP10k. Frozen CD4/subtype/barcode/background definitions; >=20 cells per sample/subtype retained. Frozen phenotype scores reused. One library per condition/time: report day-specific deltas only, no cell-as-replicate tests or time-point-as-biological-replicate P. Matching TF activity directions is not evidence for macrophage origin or causal mitochondrial damage. No external macrophage compartment analyzed in this step.\n')
 (OUT/'status.json').write_text(json.dumps(dict(status='PASS',stage='lung_scenic4_external_transfer',output=str(OUT),dataset='GSE183682',assay_libraries=4,inferential_p=False,tf=['IKZF2','STAT1']),indent=2))
 ret=ROOT/'tmp/lung_scenic4_return.tar.gz';part=Path(str(ret)+'.part')
 with tarfile.open(part,'w:gz') as t:
  for p in sorted(OUT.rglob('*')):
   if p.is_file() and not {'input','objects','__pycache__'}&set(p.relative_to(OUT).parts):t.add(p,arcname=str(p.relative_to(OUT)))
 with tarfile.open(part) as t:
  for member in t:
   if member.isfile():
    with t.extractfile(member) as f:
     while f.read(2**20):pass
 part.replace(ret);print('RETURN='+str(ret),flush=True);print('RETURN_SHA256='+hashlib.sha256(ret.read_bytes()).hexdigest(),flush=True);print('LUNG_SCENIC4_PASS',flush=True)
if __name__=='__main__':main()

PY_VALIDATE
"$PY" -m py_compile "$OUT/export.py" "$OUT/validate.py"
conda run --no-capture-output -n batch_bench python -m py_compile "$OUT/export.py"
echo SYNTAX_CHECK=PASS
"$PY" - <<'PY_PREFLIGHT'
import os,json
from pathlib import Path
r=Path(os.environ['ROOT']);p=Path(os.environ['PREV1']);s=Path(os.environ['SCENIC1'])
required=[p/'objects/GSE183682_cell_scores.tsv.gz',p/'tables/module_definition.tsv',s/'status.json',r/'02_analysis/screen2_CD4_target_wz1lfx/GSE183682/CD4_core.h5ad']
for arm in ['All','BG1']:
 required += [s/'Mouse'/arm/'score.loom',s/'Mouse'/arm/'cells.tsv']
 required += [s/'Mouse'/arm/f'seed{seed}/regulons.csv' for seed in [17,29,43]]
missing=[str(f) for f in required if not f.is_file()]
assert not missing, 'Missing frozen inputs: '+json.dumps(missing)
st=json.loads((s/'status.json').read_text())
assert st['status']=='PASS' and st['stage']=='lung_scenic1_network'
import loompy
from pyscenic.aucell import aucell
from pyscenic.cli.utils import load_signatures
print('INPUT_PREFLIGHT=PASS',flush=True)
PY_PREFLIGHT
conda run --no-capture-output -n batch_bench python -u "$OUT/export.py"
"$PY" -u "$OUT/validate.py"
