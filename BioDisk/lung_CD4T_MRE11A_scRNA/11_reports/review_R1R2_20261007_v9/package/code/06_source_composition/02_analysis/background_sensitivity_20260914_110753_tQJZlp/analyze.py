import os, re, json, tarfile, itertools
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from scipy.stats import rankdata

root=Path(os.environ['ROOT']); out=Path(os.environ['OUT'])
path=root/'02_analysis/GSE314943/stage5_20260905_222457_a2dApT/objects/CD4_core.h5ad'
a=ad.read_h5ad(path)
assert a.n_obs==2311 and a.obs_names.is_unique
assert {'sample','group','Subtype'} <= set(a.obs.columns) and 'counts' in a.layers
x=sparse.csr_matrix(a.layers['counts']); x.sum_duplicates(); x.eliminate_zeros()
assert np.isfinite(x.data).all() and (x.data>=0).all() and np.allclose(x.data,np.rint(x.data),atol=1e-6,rtol=0)
symbols=a.var['symbol'].astype('string').fillna('').str.strip().str.upper().to_numpy(dtype=str)
targets=['IL1R2','PENK','TTC39C','FOXP1','TNFRSF9','LAG3']
ig=np.array([bool(re.match(r'^(IGH[VDJMAEG]|IGK[VJC]|IGL[VJC])',s)) for s in symbols])
epi=np.isin(symbols,['SFTPA1','SFTPA2','SFTPB','SFTPC','SFTPD','SCGB1A1','SCGB3A1','SCGB3A2','ABCA3','SLC34A2','NAPSA'])
myeloid=np.isin(symbols,['LYZ1','LYZ2','S100A8','S100A9','CSF1R','LST1','TYROBP','FCER1G'])
panels={'Ig':ig,'Epi':epi,'Myeloid':myeloid,'Total':ig|epi|myeloid}
assert not np.any(panels['Total'] & np.isin(symbols,targets))
lib=np.asarray(x.sum(axis=1)).ravel(); assert (lib>0).all()
m=a.obs[['sample','group','Subtype']].astype(str).reset_index(drop=True).rename(columns={'Subtype':'subtype'})
m['cell']=a.obs_names.astype(str);m['library_UMI']=lib
assert set(m.group)=={'Saline','Silica'} and set(m.subtype)=={'Naive','ActCD4','Treg'}
samples=['Saline_1','Saline_2','Saline_3','Silica_1','Silica_2','Silica_3']
assert set(m['sample'])==set(samples)
assert (m.group==m['sample'].str.split('_').str[0]).all()
audit=[]
for panel,mask in panels.items():
    m[panel+'_UMI']=np.asarray(x[:,mask].sum(axis=1)).ravel()
    m[panel+'_frac']=m[panel+'_UMI']/lib
    for j in np.flatnonzero(mask):audit.append({'panel':panel,'gene_id':str(a.var_names[j]),'symbol':symbols[j]})
m['nonBG_UMI']=lib-m.Total_UMI
assert (m.nonBG_UMI>0).all(), 'A cell has no non-panel counts; explicit handling required'
for g in targets:
    ix=np.flatnonzero(symbols==g); assert len(ix)>0,g
    m[g]=np.asarray(x[:,ix].sum(axis=1)).ravel()
    for j in ix:audit.append({'panel':'Target','gene_id':str(a.var_names[j]),'symbol':symbols[j]})
def write(df,name):df.to_csv(out/name,sep='\t',index=False,na_rep='NA')
write(pd.DataFrame(audit),'gene_panel.tsv');write(m,'cell_burden.tsv.gz')
print('ANALYSIS_LOGIC_PASS=RAW_COUNTS_PANEL_UNION_TARGETS_EXCLUDED',flush=True)

burden=[];rows=[];corr=[]
def rho(x,y):
    if len(x)<3 or np.ptp(x)==0 or np.ptp(y)==0:return np.nan
    return float(np.corrcoef(rankdata(x),rankdata(y))[0,1])
def partial_rho(x,y,depth):
    if len(x)<20 or np.ptp(x)==0 or np.ptp(y)==0:return np.nan
    d=np.column_stack([np.ones(len(depth)),rankdata(depth)])
    rx=rankdata(x);ry=rankdata(y)
    ex=rx-d@np.linalg.lstsq(d,rx,rcond=None)[0]
    ey=ry-d@np.linalg.lstsq(d,ry,rcond=None)[0]
    return float(np.corrcoef(ex,ey)[0,1]) if np.std(ex)>1e-10 and np.std(ey)>1e-10 else np.nan
thresholds=[('All',None)]+[(f'BG{int(v*100):02d}',v) for v in [.20,.10,.05,.02,.01]]
for scope in ['Core','Naive','ActCD4','Treg']:
    z=m if scope=='Core' else m[m.subtype==scope]
    for sample in samples:
        v=z[z['sample']==sample];n=len(v);assert n>0
        group=v.group.iloc[0]
        for panel in panels:
            values=v[panel+'_frac'].to_numpy()
            burden.append(dict(scope=scope,sample=sample,group=group,panel=panel,cells=n,
                pooled_fraction=v[panel+'_UMI'].sum()/v.library_UMI.sum(),mean=values.mean(),
                median=np.median(values),p90=np.quantile(values,.9),p95=np.quantile(values,.95)))
        for g in targets:
            vals=v[g].to_numpy();dep=v.library_UMI.to_numpy()
            for panel in panels:
                b=v[panel+'_frac'].to_numpy()
                corr.append(dict(scope=scope,sample=sample,group=group,gene=g,panel=panel,cells=n,
                    rho_raw=rho(vals,b),rho_CP10k=rho(vals/dep,b),
                    rho_nonBG=rho(vals/v.nonBG_UMI.to_numpy(),b),
                    partial_rank_rho_depth=partial_rho(vals,b,dep)))
        for label,cut in thresholds:
            w=v if cut is None else v[v.Total_frac<=cut]
            k=len(w)
            for g in targets:
                count=w[g].to_numpy();den=w.library_UMI.to_numpy();nb=w.nonBG_UMI.to_numpy()
                rows.append(dict(scope=scope,sample=sample,group=group,filter=label,gene=g,
                    original_cells=n,retained=k,retention=k/n,eligible20=k>=20,
                    positive=int((count>0).sum()),target_UMI=float(count.sum()),
                    pct=100*np.mean(count>0) if k else np.nan,
                    CPM=1e6*count.sum()/den.sum() if k else np.nan,
                    CPM_nonBG=1e6*count.sum()/nb.sum() if k else np.nan,
                    CP10k=np.mean(count/den*1e4) if k else np.nan,
                    CP10k_nonBG=np.mean(count/nb*1e4) if k else np.nan))
tab=pd.DataFrame(rows)
write(pd.DataFrame(burden),'burden_by_sample.tsv');write(pd.DataFrame(corr),'within_sample_correlations.tsv')
write(tab,'sample_sensitivity.tsv')
effects=[]
for (scope,filt,gene),v in tab.groupby(['scope','filter','gene'],observed=True):
    for mincells in [1,20]:
        w=v[v.retained>=mincells]
        for metric in ['pct','CPM','CPM_nonBG','CP10k','CP10k_nonBG']:
            b=w.loc[w.group=='Saline',metric].dropna().to_numpy();c=w.loc[w.group=='Silica',metric].dropna().to_numpy()
            valid=len(b)>0 and len(c)>0;diff=c.mean()-b.mean() if valid else np.nan
            loo=[]
            if len(b)>=2 and len(c)>=2:
                loo=[c.mean()-np.delete(b,i).mean() for i in range(len(b))]+[np.delete(c,i).mean()-b.mean() for i in range(len(c))]
            effects.append(dict(scope=scope,filter=filt,gene=gene,min_cells=mincells,metric=metric,
                n_control=len(b),n_case=len(c),complete_3v3=len(b)==3 and len(c)==3,
                control_mean=b.mean() if len(b) else np.nan,case_mean=c.mean() if len(c) else np.nan,
                difference=diff,loo_same_sign_fraction=np.mean(np.sign(loo)==np.sign(diff)) if loo and diff!=0 else np.nan,
                min_retained=int(v.retained.min()),status='COMPLETE' if len(b)==3 and len(c)==3 else 'PARTIAL' if valid else 'NOT_COMPARABLE'))
ef=pd.DataFrame(effects)
baseline=ef[ef['filter']=='All'][['scope','gene','min_cells','metric','difference']].rename(columns={'difference':'baseline_difference'})
ef=ef.merge(baseline,on=['scope','gene','min_cells','metric'],how='left',validate='many_to_one')
ef['same_direction_as_all']=np.where(ef.difference.notna()&ef.baseline_difference.notna(),np.sign(ef.difference)==np.sign(ef.baseline_difference),np.nan)
write(ef,'effect_sensitivity.tsv')

# Original equal-sample subtype weights remain fixed at every threshold.
ct=m.groupby(['sample','subtype'],observed=True).size().unstack(fill_value=0).reindex(samples)
weights=ct.div(ct.sum(axis=1),axis=0).mean(axis=0)
standard=[]
for (filt,gene),v in tab[tab.scope!='Core'].groupby(['filter','gene'],observed=True):
    for sample in samples:
        z=v[v['sample']==sample].set_index('scope').reindex(weights.index)
        for metric in ['pct','CP10k','CP10k_nonBG']:
            ok=z[metric].notna().all()
            standard.append(dict(filter=filt,gene=gene,sample=sample,group=sample.split('_')[0],metric=metric,
                value=float(z[metric]@weights) if ok else np.nan,complete_subtypes=bool(ok),
                all_subtypes_eligible20=bool((z.retained>=20).all()),weights=json.dumps(weights.to_dict())))
write(pd.DataFrame(standard),'fixed_composition_sensitivity.tsv')
(out/'methods.txt').write_text(
 'Background-panel burden is a transcript fraction, not an estimated contamination probability.\n'
 'Panels: immunoglobulin family; epithelial surfactant/club markers; myeloid markers. Union avoids double counting.\n'
 'All cells retained in baseline. Fixed total-panel cutoffs <=20%,10%,5%,2%,1%; no threshold selected for significance.\n'
 'Targets are excluded from panels. nonBG denominators remove panel counts but do not correct target counts.\n'
 'Per-sample/per-subtype correlations are descriptive, no cell-level p-values. Partial rank correlations control library UMI ranks.\n'
 'All comparisons retain min_cells=1 and min_cells=20 outputs. Incomplete 3+3 comparisons are flagged, not treated as equivalent evidence.\n'
 'Leave-one-out recomputes mean differences, not edgeR models. No new DE P values or FDR are computed.\n'
 'Fixed-composition weights derive from original samples and remain unchanged across filters. Missing subtype produces NA.\n'
 'Panel genes may reflect true biology; high burden alone does not establish a doublet or justify deleting a cell from the source object.\n',encoding='utf-8')
(out/'status.json').write_text(json.dumps({'status':'PASS','source':str(path),'cells':len(m),'samples':samples,'targets':targets,'thresholds':thresholds,'analysis':'background burden descriptive sensitivity; no contamination correction or DEG refit'},indent=2))
dest=root/'tmp/lung_background_sensitivity_return.tar.gz';tmp=dest.with_name(dest.name+'.part')
with tarfile.open(tmp,'w:gz') as t:
    for p in sorted(out.iterdir()):
        if p.is_file():t.add(p,arcname=p.name)
assert tmp.stat().st_size<3*1024*1024, 'Return package exceeds 3 MiB'
os.replace(tmp,dest)
print('RETURN_PACKAGE='+str(dest),flush=True)
print('RETURN_BYTES='+str(dest.stat().st_size),flush=True)
print('BACKGROUND_SENSITIVITY_PASS',flush=True)
