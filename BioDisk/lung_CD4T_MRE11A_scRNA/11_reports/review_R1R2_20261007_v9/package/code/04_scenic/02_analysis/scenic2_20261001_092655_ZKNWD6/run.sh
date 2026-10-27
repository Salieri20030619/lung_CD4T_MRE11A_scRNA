#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA
PREV="$ROOT/02_analysis/scenic1_20260929_103513_xRoEnr"
PREV1="$ROOT/02_analysis/mechanism_step1_20260927_152218_s4V5eF"
PREV2="$ROOT/02_analysis/mechanism_step2_20260927_154919_EsN4gk"
PY="$ROOT/software/pyscenic_0.12.1_py310/bin/python"
export ROOT PREV PREV1 PREV2
export WORKERS="${WORKERS:-4}" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
bash -n "$0"
[[ -x "$PY" && -f "$PREV/status.json" ]]
mkdir -p "$ROOT/logs" "$ROOT/tmp" "$ROOT/02_analysis"
exec 9>"$ROOT/tmp/lung_scenic2.lock"
flock -n 9 || { echo SCENIC2_ALREADY_RUNNING; exit 1; }
OUT=$(mktemp -d "$ROOT/02_analysis/scenic2_$(date +%Y%m%d_%H%M%S)_XXXXXX")
export OUT
trap 'rc=$?; if ((rc != 0)); then printf "LUNG_SCENIC2_FAIL EXIT_CODE=%s OUTPUT=%s\n" "$rc" "$OUT"; fi' EXIT
mkdir -p "$OUT/tables" "$OUT/objects"
cp -- "$0" "$OUT/run.sh"
printf 'OUTPUT=%s\n' "$OUT"
cat > "$OUT/stats.py" <<'PY_STATS'
import numpy as np
import pandas as pd
from scipy import stats

def canon(g):
 a={'MRE11A':'MRE11','ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','PARK2':'PRKN'}
 s=str(g).upper();return a.get(s,s)

def design(d,model,rank_depth=False):
 cols=[np.ones(len(d)),d.case.to_numpy(float)]
 if model in ['Source','SourceDep','SourceComp']:
  for v in sorted(d.series.unique())[1:]:cols.append((d.series==v).to_numpy(float))
 if model in ['Depth','SourceDep','SourceComp','Comp']:
  x=d.library.to_numpy(float);cols.append(stats.rankdata(x) if rank_depth else np.log1p(x))
 if model in ['SourceComp','Comp']:
  comps=[c for c in d.columns if c.startswith('prop_')]
  for c in comps[:-1]:cols.append(d[c].to_numpy(float))
 return np.column_stack(cols)

def effect(y,d,model):
 y=np.asarray(y,float);assert len(y)==len(d)
 nc=int((d.case==0).sum());nt=int((d.case==1).sum())
 r=dict(n_control=nc,n_case=nt,delta=np.nan,coef=np.nan,p=np.nan,df=np.nan,reason='')
 if min(nc,nt)<2:r['reason']='FEWER_THAN_TWO_PER_GROUP';return r
 r['delta']=float(y[d.case==1].mean()-y[d.case==0].mean())
 X=design(d,model);df=len(y)-np.linalg.matrix_rank(X);r['df']=df
 if np.linalg.matrix_rank(X)<X.shape[1]:r['reason']='RANK_DEFICIENT';return r
 if df<=0:r['reason']='NO_RESIDUAL_DF';return r
 if not np.isfinite(y).all():r['reason']='NONFINITE';return r
 b=np.linalg.lstsq(X,y,rcond=None)[0];r['coef']=b[1]
 s2=float(np.sum((y-X@b)**2)/df);se=np.sqrt(max(0,np.linalg.inv(X.T@X)[1,1]*s2))
 if se>1e-15:r['p']=2*stats.t.sf(abs(b[1]/se),df)
 else:r['reason']='ZERO_RESIDUAL_VARIANCE'
 return r

def association(x,y,d,model):
 x=np.asarray(x,float);y=np.asarray(y,float);assert len(x)==len(y)==len(d)
 r=dict(n=len(x),rho=np.nan,p=np.nan,df=np.nan,reason='')
 if not np.isfinite(x).all() or not np.isfinite(y).all():r['reason']='NONFINITE';return r
 if len(np.unique(x))<2 or len(np.unique(y))<2:r['reason']='CONSTANT';return r
 if model=='Spearman':
  z=stats.spearmanr(x,y);r['rho']=z.statistic
  if len(x)>=10:r['p']=z.pvalue
  r['df']=len(x)-2
 else:
  X=design(d,model,rank_depth=True);df=len(x)-np.linalg.matrix_rank(X)-1;r['df']=df
  if np.linalg.matrix_rank(X)<X.shape[1] or df<=0:r['reason']='UNIDENTIFIABLE';return r
  rx=stats.rankdata(x);ry=stats.rankdata(y);rx=rx-X@np.linalg.lstsq(X,rx,rcond=None)[0];ry=ry-X@np.linalg.lstsq(X,ry,rcond=None)[0]
  if np.linalg.norm(rx)<1e-10 or np.linalg.norm(ry)<1e-10:r['reason']='CONSTANT_RESIDUAL';return r
  rho=float(np.clip(np.corrcoef(rx,ry)[0,1],-1,1));r['rho']=rho
  if len(x)>=10:r['p']=2*stats.t.sf(abs(rho)*np.sqrt(df/max(1-rho*rho,1e-15)),df)
 assert not np.isfinite(r['rho']) or abs(r['rho'])<=1
 return r

def bh(p):
 p=np.asarray(p,float);o=np.argsort(p);q=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];z=np.empty(len(p));z[o]=np.minimum(q,1);return z

def adjust(table,keys,kind):
 table['FDR']=np.nan
 for _,ii in table.groupby(keys,observed=True,dropna=False).groups.items():
  rows=table.loc[ii];ok=rows.p.notna()
  if kind=='effect':ok &= (rows.n_control>=10)&(rows.n_case>=10)
  else:ok &= rows.n>=10
  idx=rows.index[ok]
  if len(idx):table.loc[idx,'FDR']=bh(table.loc[idx,'p'])
 return table

PY_STATS
cat > "$OUT/analyze.py" <<'PY_ANALYZE'
import os, json, tarfile, hashlib, gc
from pathlib import Path
import numpy as np
import pandas as pd
from stats import canon,effect,association,adjust
ROOT=Path(os.environ['ROOT']);PREV=Path(os.environ['PREV']);P1=Path(os.environ['PREV1']);P2=Path(os.environ['PREV2']);OUT=Path(os.environ['OUT'])
TARGETS={'IL1R2','PENK','TTC39C','FOXP1','TNFRSF9','LAG3','MRE11'}
PRIOR={'STAT1','STAT3','STAT5A','STAT5B','NFKB1','RELA','TBX21','BATF'}
WORKERS=int(os.environ.get('WORKERS','4'))
def rd(p):return pd.read_csv(p,sep='\t')
def wr(d,name):pd.DataFrame(d).to_csv(OUT/'tables'/name,sep='\t',index=False)
def looeffect(y,d,model):
 z=[]
 for j in range(len(d)):
  keep=np.arange(len(d))!=j;q=effect(np.asarray(y)[keep],d.iloc[keep],model)
  if np.isfinite(q['coef']):z.append(q['coef'])
 return dict(loo_valid=len(z),loo_min=min(z) if z else np.nan,loo_max=max(z) if z else np.nan)
def looassoc(x,y,d,model):
 z=[]
 for j in range(len(d)):
  keep=np.arange(len(d))!=j;q=association(np.asarray(x)[keep],np.asarray(y)[keep],d.iloc[keep],model)
  if np.isfinite(q['rho']):z.append(q['rho'])
 return dict(loo_valid=len(z),loo_min=min(z) if z else np.nan,loo_max=max(z) if z else np.nan)
def sample_design(m,q,species,scope):
 z=m if scope=='Core' else m[m.subtype==scope]
 n=z.groupby('sample',observed=True).size();subjects=sorted(n[n>=20].index)
 d=z[['sample','group','series']].drop_duplicates();assert not d['sample'].duplicated().any()
 d=d.set_index('sample').loc[subjects].reset_index();d['cells']=d['sample'].map(n);d['library']=d['sample'].map(z.groupby('sample',observed=True).library.mean());d['case']=d.group.isin(['Silica','CWP']).astype(int)
 prop=pd.crosstab(z['sample'],z.subtype,normalize='index').reindex(subjects)
 for j,c in enumerate(prop.columns):d[f'prop_{j}']=prop[c].to_numpy(float)
 prevscope=('All' if scope=='Core' else scope) if species=='Mouse' else ('Core' if scope=='Core' else 'Subtype:'+scope)
 qq=q[q.scope==prevscope].set_index('sample').reindex(subjects)
 assert qq.group.notna().all() and np.array_equal(qq.group.to_numpy(),d.group.to_numpy()),'Previous score metadata mismatch'
 assert np.array_equal(qq.cells.to_numpy(),d.cells.to_numpy()),'Previous score cell counts mismatch'
 assert np.allclose(qq.library.to_numpy(),d.library.to_numpy()),'Depth mismatch'
 for c in qq.columns:
  if c.startswith('score_') or (c.startswith('expr_') and c[5:] in TARGETS):d[c]=qq[c].to_numpy()
 return d,z

def main():
 import loompy
 from ctxcore.genesig import GeneSignature
 from ctxcore.rnkdb import FeatherRankingDatabase
 from pyscenic.cli.utils import load_signatures
 from pyscenic.aucell import aucell
 status=json.loads((PREV/'status.json').read_text());assert status['status']=='PASS' and status['stage']=='lung_scenic1_network' and status['runs']==12
 mods=rd(P1/'tables/module_definition.tsv');modgenes={canon(x) for x in mods.gene}
 cov=[];effects=[];assocs=[];acts=[];edge_reports=[];summaries=[];aliasaudit=[]
 for sp in ['Mouse','Human']:
  for arm in ['All','BG1']:
   p=PREV/sp/arm;md=rd(p/'cells.tsv');assert md.cell.is_unique
   q=rd((P1/'tables/sample_scores.tsv') if sp=='Mouse' else (P2/'tables/human_sample_scores.tsv'))
   if sp=='Mouse':q=q[q.dataset=='GSE314943']
   q=q[q['filter']==arm]
   e=rd(p/'candidate_edges.tsv');e['self_edge']=e.tf==e.target
   e['species']=sp;e['filter']=arm;edge_reports.append(e)
   focus={canon(t) for t in e.loc[~e.self_edge,'tf']}|PRIOR
   gf=rd(p/'gene_filter.tsv');res=rd(PREV/'tables'/f'{sp}_resources.tsv')
   dbpaths=[Path(f) for f in res.file if str(f).endswith('.feather')]
   dbsets=[set(FeatherRankingDatabase(str(f),name=str(i)).genes) for i,f in enumerate(dbpaths)]
   assert len(dbsets)==2
   for symbol in ['MRE11','MRE11A','Mre11','Mre11a','IL1R2','Il1r2']:
    hit=gf[gf.gene==symbol]
    aliasaudit.append(dict(species=sp,filter=arm,symbol=symbol,input_present=len(hit)>0,retained=bool(hit.retained.iloc[0]) if len(hit) else False,db_10kb=symbol in dbsets[0],db_500bp=symbol in dbsets[1]))
   with loompy.connect(str(p/'score.loom'),mode='r') as ds:
    cells=np.asarray(ds.ca.CellID).astype(str);genes=np.asarray(ds.ra.Gene).astype(str);values=ds[:,:].T
   assert set(cells)==set(md.cell) and len(set(genes))==len(genes)
   expression=pd.DataFrame(values,index=cells,columns=genes).loc[md.cell]
   assert np.isfinite(expression.to_numpy()).all()
   designs={sc:sample_design(md,q,sp,sc) for sc in ['Core']+sorted(md.subtype.unique())}
   print(f'ANALYSIS_LOGIC_PASS INPUT_SCOPE_AND_SCORE_ALIGNMENT {sp} {arm}',flush=True)
   for seed in [17,29,43]:
    regs=load_signatures(str(p/f'seed{seed}'/'regulons.csv'))
    original=pd.read_csv(p/f'seed{seed}'/'auc.csv',index_col=0);original.index=original.index.astype(str)
    assert original.index.is_unique and set(original.index)==set(md.cell)
    original=original.loc[md.cell]
    for variant in ['Original','DropTarget','DropOverlap']:
     print(f'SCORE_START={sp} {arm} seed={seed} variant={variant}',flush=True)
     if variant=='Original':a=original
     else:
      signatures=[]
      for r in regs:
       exclusion=TARGETS|{canon(r.transcription_factor)}
       if variant=='DropOverlap':exclusion |= modgenes
       kept={g:w for g,w in r.gene2weight.items() if canon(g) not in exclusion}
       cov.append(dict(species=sp,filter=arm,seed=seed,variant=variant,regulon=r.name,original_genes=len(r.gene2weight),retained_genes=len(kept),status='SCORED' if kept else 'EMPTY_AFTER_EXCLUSION'))
       if kept:signatures.append(GeneSignature(name=r.name,gene2weight=kept))
      if not signatures:continue
      a=aucell(expression,signatures,seed=seed,num_workers=WORKERS,noweights=True)
      assert a.index.is_unique and set(a.index)==set(md.cell)
      a=a.loc[md.cell]
      assert np.isfinite(a.to_numpy()).all()
      a.to_csv(OUT/'objects'/f'{sp}_{arm}_{seed}_{variant}_auc.csv.gz',compression='gzip')
     for scope,(d,z) in designs.items():
      grouped=a.loc[z.cell].groupby(z.set_index('cell')['sample'],observed=True).mean().reindex(d['sample'])
      assert grouped.notna().all().all() and np.array_equal(grouped.index,d['sample'])
      models=(['Group','Depth']+(['Comp'] if scope=='Core' else [])) if sp=='Mouse' else ['Source','SourceDep']+(['SourceComp'] if scope=='Core' else [])+['WithinG227']
      for reg in grouped.columns:
       tf=reg.rsplit('(',1)[0];isfocus=canon(tf) in focus;y=grouped[reg].to_numpy()
       for i,ss in enumerate(d['sample']):acts.append(dict(species=sp,filter=arm,seed=seed,variant=variant,scope=scope,sample=ss,group=d.group.iloc[i],regulon=reg,mean_auc=y[i],cells=int(d.cells.iloc[i])))
       key=dict(species=sp,filter=arm,seed=seed,variant=variant,scope=scope,regulon=reg,focused=isfocus)
       for model in models:
        take=(d.series=='G227').to_numpy() if model=='WithinG227' else np.ones(len(d),bool)
        dd=d.loc[take].reset_index(drop=True);yy=y[take];fitmodel='Group' if model=='WithinG227' else model
        r=effect(yy,dd,fitmodel)
        if isfocus:r.update(looeffect(yy,dd,fitmodel))
        effects.append(dict(**key,model=model,**r))
       if not isfocus:continue
       for feature in [c for c in d if c.startswith(('score_','expr_'))]:
        x=d[feature].to_numpy();ok=np.isfinite(x)&np.isfinite(y);dd=d.loc[ok].reset_index(drop=True);xx=x[ok];yy=y[ok]
        for model in (['Spearman','Group']+(['Comp'] if scope=='Core' else []) if sp=='Mouse' else ['Spearman','SourceDep']+(['SourceComp'] if scope=='Core' else [])):
         r=association(xx,yy,dd,model);r.update(looassoc(xx,yy,dd,model))
         assocs.append(dict(**key,feature=feature,model=model,**r))
     print(f'SCORE_DONE={sp} {arm} seed={seed} variant={variant}',flush=True)
    summaries.append(dict(species=sp,filter=arm,seed=seed,regulons=len(regs),cells=len(md)))
   del expression,values;gc.collect()
 ef=adjust(pd.DataFrame(effects),['species','filter','seed','variant','scope','model'],'effect')
 ac=adjust(pd.DataFrame(assocs),['species','filter','seed','variant','scope','model'],'association')
 wr(ef,'TF_group_effects.tsv');wr(ac,'TF_associations.tsv');wr(cov,'regulon_exclusion_coverage.tsv');wr(aliasaudit,'gene_alias_audit.tsv');wr(pd.concat(edge_reports),'candidate_edge_audit.tsv');wr(summaries,'run_summary.tsv')
 pd.DataFrame(acts).to_csv(OUT/'tables/TF_sample_activity.tsv.gz',sep='\t',index=False,compression='gzip')
 for name,tab,value in [('effect_seed_summary.tsv',ef,'coef'),('association_seed_summary.tsv',ac,'rho')]:
  keys=['species','filter','variant','scope','regulon','model']+(['feature'] if value=='rho' else [])
  rows=[]
  for k,v in tab.groupby(keys,observed=True):
   vv=v[np.isfinite(v[value])];rows.append(dict(zip(keys,k),seeds_evaluable=vv.seed.nunique(),positive_seeds=int((vv[value]>0).sum()),negative_seeds=int((vv[value]<0).sum()),estimate_min=vv[value].min(),estimate_max=vv[value].max(),nominal_p_min=vv.p.min(),nominal_p_max=vv.p.max(),FDR_min=vv.FDR.min(),FDR_max=vv.FDR.max()))
  wr(rows,name)
 (OUT/'methods.txt').write_text('Original scores reuse SCENIC1 outputs. DropTarget removes the seven candidate genes and each TF itself from its regulon. DropOverlap additionally removes the union of frozen functional-module genes, recalculating unweighted AUCell as in SCENIC1 in the same expression/ranking universe and random seed. Empty signatures are explicitly recorded; all nonempty signatures are retained with sizes reported, not filtered by a new size threshold. All regulons receive sample-level group models; associations and leave-one-subject-out effects cover all nonself upstream TFs of the predefined candidates plus eight previous TFs, selected without reference to group p-values. Samples/subtypes use the existing >=20-cell convention. Mouse: Group, Depth and Core Comp effects (unidentifiable fits are reported explicitly); Human: Source, SourceDep, Core SourceComp, and within-G227. Correlation models additionally adjust for case group where indicated. No inferential correlation P for n<10; group BH only if both groups >=10; association BH within species/filter/seed/variant/scope/model across TF-feature pairs if n>=10. Three random seeds are not independent biological replicates. Self-target edges are flagged as algorithmically included TF genes and are not counted as regulatory evidence. MRE11 gene-alias audit does not silently alter or rerun frozen networks.\n')
 (OUT/'status.json').write_text(json.dumps(dict(status='PASS',stage='lung_scenic2_activity_sensitivity',previous=str(PREV),output=str(OUT),pending=['cross_species_orthology_validation','integration_with_macrophage_ligands','final_mechanism_literature_review']),indent=2))
 ret=ROOT/'tmp/lung_scenic2_return.tar.gz';part=Path(str(ret)+'.part')
 with tarfile.open(part,'w:gz') as t:
  for p in sorted(OUT.rglob('*')):
   if p.is_file() and 'objects' not in p.relative_to(OUT).parts and '__pycache__' not in p.parts:t.add(p,arcname=str(p.relative_to(OUT)))
 with tarfile.open(part) as t:
  for member in t:
   if member.isfile():
    with t.extractfile(member) as f:
     for b in iter(lambda:f.read(2**20),b''):pass
 part.replace(ret);print('RETURN='+str(ret),flush=True);print('RETURN_BYTES='+str(ret.stat().st_size),flush=True);print('RETURN_SHA256='+hashlib.sha256(ret.read_bytes()).hexdigest(),flush=True);print('LUNG_SCENIC2_PASS',flush=True)
if __name__=='__main__':main()

PY_ANALYZE
"$PY" -m py_compile "$OUT/stats.py" "$OUT/analyze.py"
echo SYNTAX_CHECK=PASS
"$PY" -u "$OUT/analyze.py"
