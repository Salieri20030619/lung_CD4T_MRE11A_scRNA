import os, json, itertools, tarfile, shutil
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

root = Path(os.environ['ROOT']); out = Path(os.environ['OUT']); base = root/'02_analysis'
targets = ['IL1R2','PENK','TNFRSF9','LAG3','TTC39C','FOXP1']
markers = ['CD3D','CD3E','TRAC','FOXP3','LYZ2','LYZ','TYROBP','CSF1R','IGKC','JCHAIN','IGHA','SFTPC','EPCAM']
genes = targets + markers
spec = [
 ('GSE314943','CD4','GSE314943/stage5_20260905_222457_a2dApT/objects/CD4_core.h5ad',2311),
 ('GSE314943','Mac','GSE314943/stage11_20260906_100314_pALwAB/objects/macrophages_core.h5ad',1796),
 ('GSE183682','CD4','screen2_CD4_target_wz1lfx/GSE183682/CD4_core.h5ad',1413),
 ('GSE274372','TLike','GSE274372/stage9_20260908_140721_Yj56oq/objects/CD4_author_reclustered.h5ad',888),
 ('GSE174725','CD4','screen2_Trefine_nuzVQu/GSE174725/T_refined.h5ad',10286)]
assert all((base/p).is_file() for _,_,p,_ in spec), 'Missing fixed input object'
def write(df, name):
    df.to_csv(out/name, sep='\t', index=False, na_rep='NA')
sample_rows=[]; co_rows=[]; audit=[]; mouse_ids={}; cell_inventory=[]
for dataset, lineage, rel, expected in spec:
    print('INPUT_START='+dataset+'_'+lineage, flush=True)
    a = ad.read_h5ad(base/rel)
    assert a.n_obs == expected and a.obs_names.is_unique
    assert 'counts' in a.layers and 'sample' in a.obs
    x = sparse.csr_matrix(a.layers['counts']); x.sum_duplicates(); x.eliminate_zeros()
    assert np.isfinite(x.data).all() and (x.data>=0).all()
    assert np.allclose(x.data,np.rint(x.data),atol=1e-6,rtol=0)
    library = np.asarray(x.sum(axis=1)).ravel()
    assert (library>0).all()
    ids = pd.Index(a.var_names.astype(str)).str.replace(r'\.\d+$','',regex=True)
    sy = a.var['symbol'].astype('string').fillna('').str.strip().str.upper().to_numpy(dtype=str) if 'symbol' in a.var else ids.str.upper().to_numpy()
    values = {}; found = {}
    for g in genes:
        j = np.flatnonzero(sy==g)
        if dataset=='GSE274372' and mouse_ids.get(g):
            j = np.union1d(j,np.flatnonzero(ids.isin(mouse_ids[g])))
        if dataset=='GSE314943' and lineage=='CD4':
            mouse_ids[g] = list(ids[j])
        found[g] = len(j)>0
        values[g] = np.asarray(x[:,j].sum(axis=1)).ravel() if len(j) else np.full(a.n_obs,np.nan)
        audit.append(dict(dataset=dataset,lineage=lineage,gene=g,features=';'.join(ids[j]),found=found[g],path=str(base/rel)))
    m = pd.DataFrame({'cell':a.obs_names.astype(str),'sample':a.obs['sample'].astype(str).to_numpy(),'library':library,**values})
    if dataset=='GSE314943':
        m['subtype']=a.obs['Subtype' if lineage=='CD4' else 'MacType'].astype(str).to_numpy()
        m['group']=a.obs['group'].astype(str).to_numpy(); m['day']=-1
        assert set(m.group)=={'Saline','Silica'} and m['sample'].nunique()==6
    elif dataset=='GSE183682':
        mapping={'GSM5567501':('NS',7),'GSM5567502':('Silica',7),'GSM5567503':('NS',56),'GSM5567504':('Silica',56)}
        m['group']=m['sample'].map({k:v[0] for k,v in mapping.items()})
        m['day']=m['sample'].map({k:v[1] for k,v in mapping.items()})
        m['subtype']=a.obs.F08.astype(str).map({'0':'NaiveCD4','3':'ActCD4','5':'Treg'}).to_numpy()
        assert m.subtype.notna().all()
    elif dataset=='GSE274372':
        m['day']=pd.to_numeric(a.obs.day).to_numpy()
        m['group']='D'+m.day.astype(int).astype(str)
        m['subtype']=a.obs.R08.astype(str).map({'0':'ActReg','1':'TLike1'}).to_numpy()
        m=m[m.subtype.notna()].copy(); assert len(m)==494
    else:
        m['subtype']=a.obs.F08.astype(str).map({'0':'MemCD4','1':'NaiveCD4','12':'Treg'}).to_numpy()
        m=m[m.subtype.notna()].copy(); assert len(m)==2958
        m['group']=m['sample'].map({'GSM5324866':'Exposed','GSM5324867':'Exposed','GSM5324868':'Silicosis','GSM5324869':'Silicosis','GSM5324870':'Silicosis'})
        m['day']=-1
    assert m[['sample','group','subtype','day']].notna().all().all()
    for sample, v in m.groupby('sample',observed=True):
        for subtype in ['All']+sorted(m.subtype.unique()):
            z=v if subtype=='All' else v[v.subtype==subtype]
            n=len(z); denom=z.library.sum()
            common=dict(dataset=dataset,lineage=lineage,subtype=subtype,sample=sample,group=v.group.iloc[0],day=v.day.iloc[0],cells=n)
            cell_inventory.append(common)
            for gene in targets:
                arr=z[gene].to_numpy(); ok=found[gene] and n>0
                sample_rows.append(dict(**common,gene=gene,found=found[gene],positive=int((arr>0).sum()) if ok else np.nan,
                    pct=100*np.mean(arr>0) if ok else np.nan,target_UMI=arr.sum() if ok else np.nan,
                    CPM=1e6*arr.sum()/denom if ok else np.nan,
                    CP10k=np.mean(arr/z.library.to_numpy()*1e4) if ok else np.nan,
                    mean_library=denom/n if n else np.nan,eligible20=n>=20))
                if not ok: continue
                pos=z[gene]>0
                for marker in markers:
                    measured=found[marker]
                    co_rows.append(dict(**common,gene=gene,marker=marker,marker_found=measured,
                        target_positive=int(pos.sum()),both_positive=int((pos & (z[marker]>0)).sum()) if measured else np.nan,
                        marker_pct_in_targetpos=100*np.mean(z.loc[pos,marker]>0) if measured and pos.any() else np.nan,
                        marker_pct_in_targetneg=100*np.mean(z.loc[~pos,marker]>0) if measured and (~pos).any() else np.nan))
    del a,x,m
    print('INPUT_DONE='+dataset+'_'+lineage,flush=True)
tab=pd.DataFrame(sample_rows)
assert not tab.duplicated(['dataset','lineage','subtype','sample','gene']).any()
write(tab,'sample_expression.tsv'); write(pd.DataFrame(co_rows),'marker_codetection.tsv')
write(pd.DataFrame(audit),'feature_audit.tsv');write(pd.DataFrame(cell_inventory),'cell_inventory.tsv')

# Each comparison gives equal weight to measured independent samples/libraries.
# Leave-one-out statistics describe direction robustness and are not new tests.
comparisons=[]
for keys,v in tab.groupby(['dataset','lineage','subtype','gene'],observed=True):
    dataset,lineage,subtype,gene=keys
    pairs=[('Saline','Silica',-1)] if dataset=='GSE314943' else [('NS','Silica',7),('NS','Silica',56)] if dataset=='GSE183682' else [('D0','D7',-1),('D0','D28',-1),('D0','D56',-1)] if dataset=='GSE274372' else [('Exposed','Silicosis',-1)]
    for control,case,day in pairs:
        for threshold in [1,20]:
            z=v[(v.cells>=threshold)&v.found]
            if day!=-1:z=z[z.day==day]
            for metric in ['pct','CPM','CP10k']:
                b=z.loc[z.group==control,metric].dropna().to_numpy(); c=z.loc[z.group==case,metric].dropna().to_numpy()
                effect=c.mean()-b.mean() if len(b) and len(c) else np.nan
                loo=[]
                if len(b)>=2 and len(c)>=2:
                    loo=[c.mean()-np.delete(b,i).mean() for i in range(len(b))]+[np.delete(c,i).mean()-b.mean() for i in range(len(c))]
                comparisons.append(dict(dataset=dataset,lineage=lineage,subtype=subtype,gene=gene,control=control,case=case,day=day,min_cells=threshold,metric=metric,
                    n_control=len(b),n_case=len(c),control_mean=b.mean() if len(b) else np.nan,case_mean=c.mean() if len(c) else np.nan,
                    difference=effect,loo_same_sign_fraction=np.mean(np.sign(loo)==np.sign(effect)) if loo and effect!=0 else np.nan,
                    loo_min=min(loo) if loo else np.nan,loo_max=max(loo) if loo else np.nan,
                    pair_case_higher=np.mean(c[:,None]>b[None,:]) if len(b) and len(c) else np.nan,
                    max_case_expression_share=c.max()/c.sum() if len(c) and c.sum()>0 else np.nan))
write(pd.DataFrame(comparisons),'sample_consistency.tsv')

# Standardize main CD4 subtype composition using fixed equal-sample weights.
std=[]
q=tab[(tab.dataset=='GSE314943')&(tab.lineage=='CD4')&(tab.subtype!='All')]
for gene,v in q.groupby('gene',observed=True):
    ct=v.pivot(index='sample',columns='subtype',values='cells')
    w=ct.div(ct.sum(axis=1),axis=0).mean(axis=0)
    for metric in ['pct','CP10k']:
        mat=v.pivot(index='sample',columns='subtype',values=metric).reindex(columns=w.index)
        for sample,row in mat.iterrows():
            std.append(dict(gene=gene,sample=sample,group=v.loc[v['sample']==sample,'group'].iloc[0],metric=metric,
                value=float(row@w) if row.notna().all() else np.nan,weights=json.dumps(w.to_dict())))
write(pd.DataFrame(std),'main_fixed_composition.tsv')

# Extract existing sample-level TMM logCPM and normalization metadata; no refit.
pb=base/'GSE314943/stage6_20260905_231807_L44XdT/tables'
pbrows=[];pb_status=[]
id_to_gene={i:g for g in targets for i in mouse_ids.get(g,[])}
for scope in ['Core','ActCD4','Naive','Treg']:
    for arm in ['Full','NoBG']:
        p=pb/f'{scope}_{arm}_logCPM.tsv'
        assert p.is_file(),str(p)
        z=pd.read_csv(p,sep='\t',dtype={'gene_id':str})
        z['gene']=z.gene_id.map(id_to_gene)
        for gene in targets:
            pb_status.append(dict(scope=scope,arm=arm,gene=gene,tested=bool((z.gene==gene).any())))
        z=z[z.gene.notna()]
        z=z.melt(id_vars=['gene_id','gene'],var_name='sample',value_name='TMM_logCPM')
        z['scope']=scope;z['arm']=arm;pbrows.append(z)
        n=pb/f'{scope}_{arm}_normalization.tsv'
        if n.is_file():shutil.copy2(n,out/n.name)
write(pd.concat(pbrows,ignore_index=True),'main_TMM_logCPM.tsv')
write(pd.DataFrame(pb_status),'main_TMM_eligibility.tsv')
(out/'methods.txt').write_text(
    'Six prespecified targets. Raw counts; duplicate gene symbols sum; missing features remain NA.\n'
    'Independent samples/libraries are weighted equally. No cell-level significance tests.\n'
    'Min_cells=1 is descriptive; min_cells=20 is sensitivity, retaining the original all-sample output.\n'
    'GSE183682 has one library per condition/time. GSE274372 is TLike and has no late matched controls.\n'
    'GSE174725 uses exposed rather than healthy comparators. Human and mouse are never pooled.\n'
    'Marker co-detection describes cellular provenance; absence of T markers can reflect dropout. No automatic cell deletion.\n'
    'Mac reference is limited to main-cohort macrophages, not a whole-lung specificity test.\n'
    'CPM=pooled target counts / pooled total counts * 1e6. CP10k=mean per-cell normalized counts.\n'
    'Leave-one-out evaluates direction under removal of each sample without refitting edgeR.\n'
    'TMM logCPM is extracted from the original edgeR output, with original normalization metadata.\n',encoding='utf-8')
print('ANALYSIS_LOGIC_PASS=COUNTS_SCOPES_SAMPLE_UNIT_MISSING_FEATURES',flush=True)
archive=root/'tmp/lung_candidate6_verify_return.tar.gz';tmp=archive.with_name(archive.name+'.part')
with tarfile.open(tmp,'w:gz') as t:
    for p in sorted(out.iterdir()):
        if p.is_file():t.add(p,arcname=p.name)
assert tmp.stat().st_size<3*1024*1024,'Return package exceeds 3 MiB; keep local output'
os.replace(tmp,archive)
print('RETURN_PACKAGE='+str(archive),flush=True)
print('RETURN_BYTES='+str(archive.stat().st_size),flush=True)
print('CANDIDATE6_VERIFY_PASS',flush=True)
