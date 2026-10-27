import os,json,gzip,hashlib,itertools
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from shared import TARGETS,SETS,canonical,bgmask,scores,effect

root=Path(os.environ['ROOT']);out=Path(os.environ['OUT']);prev=Path(os.environ['PREV'])
def wr(df,name):df.to_csv(out/'tables'/name,sep='\t',index=False,na_rep='NA')
def matrix(path):
    with gzip.open(path,'rb') if str(path).endswith('.gz') else open(path,'rb') as h:return sparse.csr_matrix(mmread(h).T,dtype=float)
def check(x,m):
    x.sum_duplicates();x.eliminate_zeros()
    assert x.shape[0]==len(m) and m.cell.is_unique
    assert np.isfinite(x.data).all() and (x.data>=0).all() and np.allclose(x.data,np.rint(x.data),atol=1e-6,rtol=0)
    assert (m.library>0).all()
    assert (np.asarray(x.sum(axis=1)).ravel()<=m.library.to_numpy()+1e-6).all()

# Human counts retain original full-library denominators and exact subject/source labels.
hm=pd.read_csv(out/'input/human_cells.tsv',sep='\t');raw=matrix(out/'input/human_counts.mtx');check(raw,hm)
hm['core']=hm.core.astype(str).str.upper().map({'TRUE':True,'FALSE':False});assert hm.core.notna().all()
original=pd.read_csv(out/'input/human_genes.tsv',sep='\t').symbol.astype('string').fillna('')
sy=np.array([canonical(s) for s in original]);valid=(sy!='')&~np.char.startswith(sy,'ENS')
symbols=np.array(sorted(set(sy[valid])));assert len(symbols)>5000
proj=sparse.csr_matrix((np.ones(valid.sum()),(np.flatnonzero(valid),pd.Index(symbols).get_indexer(sy[valid]))),shape=(len(sy),len(symbols)))
hx=(raw@proj).tocsr();lib=hm.library.to_numpy()
wr(pd.DataFrame({'original':original,'canonical':sy,'retained':valid}),'human_gene_mapping.tsv.gz')
hm['bg']=np.asarray(hx[:,bgmask(symbols)].sum(axis=1)).ravel()/lib
hm['mt']=np.asarray(hx[:,np.char.startswith(symbols,'MT-')].sum(axis=1)).ravel()/lib
hm['nfeature']=np.asarray((hx>0).sum(axis=1)).ravel()
y=hx.multiply((1e4/lib)[:,None]).tocsr();y.data=np.log1p(y.data)
vals,det,cov,ctrl=scores(y,hx,symbols,hm['sample'].to_numpy())
wr(pd.DataFrame(cov),'human_module_coverage.tsv');wr(pd.DataFrame(ctrl),'human_score_controls.tsv.gz')
for mod in SETS:hm['score_'+mod]=vals[mod];hm['detect_'+mod]=det[mod]
for gene in TARGETS:
    j=np.flatnonzero(symbols==gene);a=np.asarray(hx[:,j].sum(axis=1)).ravel() if len(j) else np.full(len(hm),np.nan)
    hm['expr_'+gene]=a/lib*1e4;hm['positive_'+gene]=a>0 if len(j) else np.nan
rows=[]
for filt,cut in [('All',1.),('BG1',.01)]:
    for scope in ['Primary','Core']+['Subtype:'+s for s in sorted(hm.subtype.unique())]:
        mask=np.ones(len(hm),bool) if scope=='Primary' else hm.core.to_numpy() if scope=='Core' else (hm.subtype==scope.split(':',1)[1]).to_numpy()
        for ss in sorted(hm['sample'].unique()):
            allm=hm[hm['sample']==ss];z=hm[mask&(hm['sample']==ss)&(hm.bg<=cut)];n=len(z)
            d=dict(scope=scope,filter=filt,sample=ss,group=allm.group.iloc[0],series=allm.series.iloc[0],cells=n,original_cells=int(np.sum(mask&(hm['sample']==ss))),eligible20=n>=20,library=z.library.mean(),bg=z.bg.mean(),mt=z.mt.mean())
            for k in SETS:d['score_'+k]=z['score_'+k].mean();d['detect_'+k]=z['detect_'+k].mean()
            for k in TARGETS:d['expr_'+k]=z['expr_'+k].mean();d['pct_'+k]=z['positive_'+k].mean()*100
            rows.append(d)
ht=pd.DataFrame(rows);wr(ht,'human_sample_scores.tsv')
hm.to_csv(out/'objects/human_cell_scores.tsv.gz',sep='\t',index=False,na_rep='NA')
# Fixed subtype weights shared across groups/filters within each human scope.
std=[]
for sc in ['Primary','Core']:
    b=hm if sc=='Primary' else hm[hm.core]
    ct=pd.crosstab(b['sample'],b.subtype);weights=ct.div(ct.sum(axis=1),axis=0).mean()
    for filt in ['All','BG1']:
        for ss in sorted(b['sample'].unique()):
            z=ht[(ht['filter']==filt)&(ht['sample']==ss)&ht.scope.str.startswith('Subtype:')].copy();z['type']=z.scope.str.replace('Subtype:','',regex=False);z=z.set_index('type').reindex(weights.index)
            for metric in ['score_'+k for k in SETS]+['expr_'+g for g in TARGETS]:
                std.append(dict(scope=sc,filter=filt,sample=ss,group=z.group.iloc[0],series=z.series.iloc[0],metric=metric,value=float(z[metric]@weights) if z[metric].notna().all() else np.nan,all_subtypes_eligible20=bool((z.cells>=20).all()),weights=json.dumps(weights.to_dict())))
wr(pd.DataFrame(std),'human_fixed_composition.tsv')

# Decompose main mouse and human panel scores into measured gene contributions.
mm=pd.read_csv(prev/'comm_input/CD4_cells.tsv',sep='\t');mx=matrix(prev/'comm_input/CD4_counts.mtx.gz');check(mx,mm)
ms=pd.read_csv(prev/'comm_input/CD4_genes.tsv',sep='\t').canonical.astype(str).to_numpy();assert mx.shape[1]==len(ms)
mm['core']=True;mm['series']='G314';mm['bg']=np.asarray(mx[:,bgmask(ms)].sum(axis=1)).ravel()/mm.library.to_numpy()
rawrows=[];contrib=[];joint=[];lomo=[]
extra=['GSDME','CASP3','CASP1','GSDMD','PYCARD','IL1B','IL18','NLRP3']
for ds,x,meta,genes,ct in [('GSE314943',mx,mm,ms,pd.read_csv(prev/'tables/score_controls.tsv.gz',sep='\t').query("dataset == 'GSE314943'")),('HumanLung',hx,hm,symbols,pd.DataFrame(ctrl))]:
    log=x.multiply((1e4/meta.library.to_numpy())[:,None]).tocsr();log.data=np.log1p(log.data)
    targets=sorted(set(TARGETS+extra+list(itertools.chain.from_iterable(SETS.values()))));counts={}
    for g in targets:
        j=np.flatnonzero(genes==g);counts[g]=np.asarray(x[:,j].sum(axis=1)).ravel() if len(j) else np.full(len(meta),np.nan)
    available={g:np.any(genes==g) for g in targets}
    ctrls={k:np.flatnonzero(np.isin(genes,ct.loc[ct.module==k,'gene'])) for k in SETS}
    ctrlvalues={k:np.asarray(log[:,j].mean(axis=1)).ravel() if len(j) else np.full(len(meta),np.nan) for k,j in ctrls.items()}
    subtypes=sorted(meta.subtype.unique());scopes=['Primary','Core']+['Subtype:'+k for k in subtypes] if ds=='HumanLung' else ['Core']+['Subtype:'+k for k in subtypes]
    for filt,cut in [('All',1.),('BG1',.01)]:
        for sc in scopes:
            sm=np.ones(len(meta),bool) if sc=='Primary' else meta.core.to_numpy() if sc=='Core' else (meta.subtype==sc.split(':',1)[1]).to_numpy()
            for ss in sorted(meta['sample'].unique()):
                ix=np.flatnonzero(sm&(meta['sample']==ss)&(meta.bg<=cut));n=len(ix);base=meta[meta['sample']==ss].iloc[0]
                key=dict(dataset=ds,scope=sc,filter=filt,sample=ss,group=base.group,series=base.series,cells=n,eligible20=n>=20)
                for g in targets:
                    v=counts[g][ix];ok=available[g] and n>0
                    rawrows.append(dict(**key,gene=g,present=available[g],positive=int(np.sum(v>0)) if ok else np.nan,pct=100*np.mean(v>0) if ok else np.nan,CP10k=np.mean(v/meta.library.to_numpy()[ix]*1e4) if ok else np.nan,mean_log=np.mean(np.log1p(v/meta.library.to_numpy()[ix]*1e4)) if ok else np.nan))
                for label,gs in [('CaspGsd',['CASP1','GSDMD']),('Pyro3',['CASP1','GSDMD','PYCARD']),('AltPyro',['CASP3','GSDME'])]:
                    ok=all(available[g] for g in gs) and n>0
                    pos=np.logical_and.reduce([counts[g][ix]>0 for g in gs]) if ok else None
                    joint.append(dict(**key,pair=label,genes=';'.join(gs),positive=int(pos.sum()) if ok else np.nan,pct=100*pos.mean() if ok else np.nan))
                for mod in ['MitoISR','MitoROS','PyroExec','Inflam','OXPHOS']:
                    gs=[g for g in SETS[mod] if available[g]];nmod=len(gs)
                    if not n or not nmod:continue
                    bg=float(np.mean(ctrlvalues[mod][ix]));means={g:np.mean(np.log1p(counts[g][ix]/meta.library.to_numpy()[ix]*1e4)) for g in gs};full=float(np.mean(list(means.values()))-bg)
                    for g,v in means.items():
                        contrib.append(dict(**key,module=mod,gene=g,mean_log=v,control_mean=bg,contribution=(v-bg)/nmod,module_score=full,module_genes=nmod))
                        lomo.append(dict(**key,module=mod,removed_gene=g,full_score=full,reduced_score=(sum(means.values())-v)/(nmod-1)-bg if nmod>1 else np.nan))
wr(pd.DataFrame(rawrows),'panel_gene_expression.tsv');wr(pd.DataFrame(contrib),'module_gene_contributions.tsv');wr(pd.DataFrame(joint),'pyro_codetection.tsv');wr(pd.DataFrame(lomo),'module_leave_one_gene.tsv')
# Exact additive decomposition audit against previous mouse panel scores.
c=pd.DataFrame(contrib);prevtab=pd.read_csv(prev/'tables/sample_scores.tsv',sep='\t');errors=[]
for (scope,filt,ss,mod),v in c[c.dataset=='GSE314943'].groupby(['scope','filter','sample','module'],observed=True):
    old_scope='All' if scope=='Core' else scope.split(':',1)[1]
    q=prevtab[(prevtab.dataset=='GSE314943')&(prevtab.scope==old_scope)&(prevtab['filter']==filt)&(prevtab['sample']==ss)]
    assert len(q)==1
    observed=float(v.contribution.sum());expected=float(q.iloc[0]['score_'+mod])
    errors.append(dict(scope=scope,filter=filt,sample=ss,module=mod,absolute_error=abs(observed-expected)))
assert max(z['absolute_error'] for z in errors)<1e-8,'Score decomposition mismatch'
wr(pd.DataFrame(errors),'decomposition_audit.tsv')
# Correlation leverage in the six mouse samples: direction and range after deletion.
from scipy.stats import rankdata

def residual_corr(a,b,group):
    a=np.asarray(a,float);b=np.asarray(b,float);group=np.asarray(group,float)
    ok=np.isfinite(a)&np.isfinite(b)&np.isfinite(group);a=a[ok];b=b[ok];group=group[ok]
    if len(a)<5 or np.ptp(a)==0 or np.ptp(b)==0 or len(set(group))!=2:return np.nan
    d=np.column_stack([np.ones(len(a)),group]);a=rankdata(a);b=rankdata(b);a-=d@np.linalg.lstsq(d,a,rcond=None)[0];b-=d@np.linalg.lstsq(d,b,rcond=None)[0]
    if np.std(a)<1e-10 or np.std(b)<1e-10:return np.nan
    return float(np.clip(np.corrcoef(a,b)[0,1],-1,1))
leverage=[]
for (scope,filt),v in prevtab[(prevtab.dataset=='GSE314943')&prevtab.eligible20].groupby(['scope','filter'],observed=True):
    for g in TARGETS:
        for mod in SETS:
            z=v[['sample','group','expr_'+g,'score_'+mod]].dropna()
            if len(z)!=6:continue
            full=residual_corr(z['expr_'+g],z['score_'+mod],(z.group=='Silica').astype(int))
            for ss in z['sample']:
                q=z[z['sample']!=ss];val=residual_corr(q['expr_'+g],q['score_'+mod],(q.group=='Silica').astype(int))
                leverage.append(dict(scope=scope,filter=filt,gene=g,module=mod,removed_sample=ss,full_rho=full,loo_rho=val,same_sign=bool(np.sign(val)==np.sign(full)) if np.isfinite(val) and np.isfinite(full) and full!=0 else np.nan))
wr(pd.DataFrame(leverage),'mouse_correlation_leverage.tsv')
# Audit known NicheNet cache files without network or package changes.
cache=[]
for directory in [root/'00_raw',root/'02_analysis',root/'08_MRE11A',root/'resources']:
    if directory.exists():
        for dp,dirs,files in os.walk(directory):
            dirs[:]=[d for d in dirs if d not in ['comm_input','node_modules','.git']]
            for name in files:
                if name.endswith(('.rds','.RDS')) and any(k in name.lower() for k in ['ligand_target','lr_network','weighted_network','nichenet']):
                    p=Path(dp)/name;cache.append(dict(path=str(p),size_bytes=p.stat().st_size))
wr(pd.DataFrame(cache,columns=['path','size_bytes']),'nichenet_cache_inventory.tsv')
(out/'python_status.json').write_text(json.dumps({'status':'PASS','human_cells':len(hm),'human_subjects':hm['sample'].nunique(),'mouse_cells':len(mm),'max_decomposition_error':max(z['absolute_error'] for z in errors)},indent=2))
print('ANALYSIS_LOGIC_PASS=HUMAN_SCOPE_SOURCE_COUNTS_AND_GENE_DECOMPOSITION',flush=True)

