import os, re, json, gzip, tarfile, hashlib, platform, itertools
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata
from scipy.io import mmwrite

TARGETS = ['IL1R2','PENK','TTC39C','FOXP1','TNFRSF9','LAG3','MRE11']
SETS = {
 'OXPHOS':'NDUFA1 NDUFA2 NDUFA9 NDUFB8 NDUFS1 NDUFS2 NDUFS3 NDUFV1 SDHA SDHB SDHC UQCRC1 UQCRC2 UQCRFS1 CYC1 COX4I1 COX5A COX5B COX6A1 ATP5F1A ATP5F1B ATP5F1C ATP5MC1 ATP5PO',
 'MitoISR':'ATF4 ATF5 DDIT3 DDIT4 CHAC1 ASNS TRIB3 PPP1R15A SESN2 DELE1 EIF2AK1 HSPA9 HSPD1 HSPE1 LONP1 CLPP OMA1',
 'MitoROS':'SOD2 PRDX3 PRDX5 GPX4 TXN2 TXNRD2 HMOX1 NQO1 GCLC GCLM SRXN1 FTH1',
 'Mitophagy':'PINK1 PRKN BNIP3 BNIP3L FUNDC1 SQSTM1 OPTN TBK1 ULK1 BECN1 MAP1LC3B',
 'Inflam':'NLRP3 AIM2 NLRC4 PYCARD CASP1',
 'PyroExec':'PYCARD CASP1 GSDMD IL1B IL18',
 'ApoExec':'BAX BAK1 APAF1 CASP9 CASP3 CASP7 CYCS',
 'TAct':'CD69 IL2RA ICOS TNFRSF4 BATF IRF4 NR4A1 NR4A2 EGR1',
 'TReg':'FOXP3 IL2RA CTLA4 IKZF2 IKZF4 LRRC32 ENTPD1 TIGIT',
 'Inhib':'PDCD1 CTLA4 TIGIT HAVCR2 TOX TOX2 ENTPD1',
 'Naive':'CCR7 SELL TCF7 LEF1 KLF2 IL7R MAL',
 'Cycle':'MKI67 TOP2A PCNA MCM2 MCM3 MCM4 MCM5 MCM6 MCM7 UBE2C'
}
SETS={k:v.split() for k,v in SETS.items()}
EPI='SFTPA1 SFTPA2 SFTPB SFTPC SFTPD SCGB1A1 SCGB3A1 SCGB3A2 ABCA3 SLC34A2 NAPSA'.split()
MYE='LYZ LYZ1 LYZ2 S100A8 S100A9 CSF1R LST1 TYROBP FCER1G'.split()
SOURCE_GENES='IL1A IL1B IL18 CXCL16 TNF IL6 IL10 TGFB1 IL12A IL12B IL23A CCL3 CCL4 CCL5 SPP1 LGALS9 CD80 CD86 CD274 PDCD1LG2 TNFSF9'.split()
RECEPTORS='IL1R1 IL1RAP IL1R2 IL18R1 IL18RAP CXCR6 TNFRSF1A TNFRSF1B IL6R IL6ST IL10RA IL10RB TGFBR1 TGFBR2 IL12RB1 IL12RB2 IL23R CCR5 CD44 ITGAV ITGB1 ITGB3 HAVCR2 CD28 CTLA4 PDCD1 TNFRSF9'.split()
SPECS=[
 ('GSE314943','CD4','GSE314943/stage5_20260905_222457_a2dApT/objects/CD4_core.h5ad',2311),
 ('GSE314943','Mac','GSE314943/stage11_20260906_100314_pALwAB/objects/macrophages_core.h5ad',1796),
 ('GSE183682','CD4','screen2_CD4_target_wz1lfx/GSE183682/CD4_core.h5ad',1413),
 ('GSE274372','TLike','GSE274372/stage9_20260908_140721_Yj56oq/objects/CD4_author_reclustered.h5ad',888),
 ('GSE174725','CD4','screen2_Trefine_nuzVQu/GSE174725/T_refined.h5ad',10286)]

def canonical(s):
    s=str(s).strip().upper()
    return {'MRE11A':'MRE11','ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','PARK2':'PRKN'}.get(s,s)

def rho(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float);ok=np.isfinite(a)&np.isfinite(b);a=a[ok];b=b[ok]
    if len(a)<3 or np.ptp(a)==0 or np.ptp(b)==0:return np.nan
    return float(np.clip(np.corrcoef(rankdata(a),rankdata(b))[0,1],-1,1))

def partial(a,b,z):
    a=np.asarray(a,float);b=np.asarray(b,float);z=np.asarray(z,float)
    if z.ndim==1:z=z[:,None]
    ok=np.isfinite(a)&np.isfinite(b)&np.isfinite(z).all(axis=1);a=a[ok];b=b[ok];z=z[ok]
    if len(a)<6 or np.ptp(a)==0 or np.ptp(b)==0:return np.nan
    d=np.column_stack([np.ones(len(a))]+[rankdata(z[:,j]) for j in range(z.shape[1])])
    if len(a)-np.linalg.matrix_rank(d)<3:return np.nan
    u=rankdata(a);v=rankdata(b);u-=d@np.linalg.lstsq(d,u,rcond=None)[0];v-=d@np.linalg.lstsq(d,v,rcond=None)[0]
    if np.std(u)<1e-10 or np.std(v)<1e-10:return np.nan
    return float(np.clip(np.corrcoef(u,v)[0,1],-1,1))

def effect(b,c):
    b=np.asarray(b,float);c=np.asarray(c,float);b=b[np.isfinite(b)];c=c[np.isfinite(c)]
    d=float(c.mean()-b.mean()) if len(b) and len(c) else np.nan
    loo=[]
    if min(len(b),len(c))>=2:
        loo=[c.mean()-np.delete(b,i).mean() for i in range(len(b))]+[np.delete(c,i).mean()-b.mean() for i in range(len(c))]
    return dict(n_control=len(b),n_case=len(c),control_mean=b.mean() if len(b) else np.nan,case_mean=c.mean() if len(c) else np.nan,
      difference=d,loo_same=np.mean(np.sign(loo)==np.sign(d)) if loo and d!=0 else np.nan,
      loo_min=min(loo) if loo else np.nan,loo_max=max(loo) if loo else np.nan)

def bgmask(symbols):
    return np.array([bool(re.match(r'^(IGH[VDJMAEG]|IGK[VJC]|IGL[VJC])',s)) or s in EPI or s in MYE for s in symbols])

def scores(y, x, symbols, samples):
    # Equal-sample reference expression; controls fixed before group comparisons.
    reference=np.mean([np.asarray(y[np.asarray(samples)==s].mean(axis=0)).ravel() for s in sorted(set(samples))],axis=0)
    excluded=set(TARGETS)|set(itertools.chain.from_iterable(SETS.values()))|set(EPI)|set(MYE)
    pool=np.array([s not in excluded and not re.match(r'^(MT-|RPL|RPS|MRPL|MRPS|IGH|IGK|IGL)',s) for s in symbols])&(reference>0)
    edges=np.unique(np.quantile(reference[reference>0],np.linspace(0,1,26)))
    bins=np.digitize(reference,edges[1:-1]);coverage=[];controls=[];values={};members={}
    for name,expected in SETS.items():
        j=np.flatnonzero(np.isin(symbols,expected));observed=j[np.asarray(x[:,j].sum(axis=0)).ravel()>0]
        present=list(symbols[j]);eligible=len(j)>=3 and len(j)/len(expected)>=.6 and len(observed)>=3
        selected=set()
        for ix in j:
            same=np.flatnonzero(pool & (bins==bins[ix]))
            options=same if len(same)>=20 else np.flatnonzero(pool)
            order=np.argsort(np.abs(reference[options]-reference[ix]),kind='stable')
            selected.update(options[order[:20]].tolist())
        cc=np.array(sorted(selected),dtype=int)
        eligible=eligible and len(cc)>=20
        coverage.append(dict(module=name,expected=len(expected),present=len(j),expressed=len(observed),coverage=len(j)/len(expected),eligible=eligible,
          genes=';'.join(present),missing=';'.join(sorted(set(expected)-set(present))),controls=len(cc)))
        for ix in cc:controls.append(dict(module=name,gene=symbols[ix],reference_mean=reference[ix],bin=int(bins[ix])))
        mean=np.asarray(y[:,j].mean(axis=1)).ravel() if len(j) else np.full(y.shape[0],np.nan)
        values[name]=mean-np.asarray(y[:,cc].mean(axis=1)).ravel() if eligible else np.full(y.shape[0],np.nan)
        members[name]=np.asarray((x[:,j]>0).sum(axis=1)).ravel() if len(j) else np.zeros(y.shape[0])
    return values,members,coverage,controls

def main():
    import anndata as ad
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    root=Path(os.environ['ROOT']);out=Path(os.environ['OUT']);base=root/'02_analysis'
    for sub in ['tables','figures','objects','comm_input']:(out/sub).mkdir(exist_ok=True)
    def wr(df,name):df.to_csv(out/'tables'/name,sep='\t',index=False,na_rep='NA')
    missing=[str(base/p) for _,_,p,_ in SPECS if not (base/p).is_file()]
    if missing:raise FileNotFoundError('Missing frozen inputs: '+json.dumps(missing))
    mouse_map={};mouse_original={};mouse_names={};main_symbols=None;inventory=[];coverage_all=[];control_all=[];mapping_all=[];sample_rows=[];within=[];source=[];paircounts=[];exports=[];frames={}
    for dataset,lineage,rel,expected in SPECS:
        print('DATASET_START='+dataset+'_'+lineage,flush=True)
        a=ad.read_h5ad(base/rel)
        assert a.n_obs==expected and a.obs_names.is_unique,(rel,a.shape)
        assert 'counts' in a.layers and 'sample' in a.obs
        raw=sparse.csr_matrix(a.layers['counts'],dtype=float);raw.sum_duplicates();raw.eliminate_zeros()
        assert np.isfinite(raw.data).all() and (raw.data>=0).all() and np.allclose(raw.data,np.rint(raw.data),rtol=0,atol=1e-6)
        lib=np.asarray(raw.sum(axis=1)).ravel();assert (lib>0).all()
        ids=pd.Index(a.var_names.astype(str)).str.replace(r'\.\d+$','',regex=True)
        original=a.var['symbol'].astype('string').fillna('').astype(str).to_numpy() if 'symbol' in a.var else np.array(['']*a.n_vars)
        if dataset=='GSE314943' and lineage=='CD4':
            for gi,sy in zip(ids,original):
                if sy.strip() and not sy.upper().startswith('ENS'):
                    if gi in mouse_map:assert mouse_map[gi]==canonical(sy)
                    mouse_map[gi]=canonical(sy);mouse_original[gi]=sy;mouse_names.setdefault(canonical(sy),sy)
        sy=[]
        for gi,s in zip(ids,original):
            if dataset=='GSE274372' and gi in mouse_map:s=mouse_map[gi]
            elif not s.strip() or s.upper().startswith('ENS'):s='' if gi.upper().startswith('ENS') else gi
            sy.append(canonical(s))
        sy=np.array(sy);valid=(sy!='')&~np.char.startswith(sy,'ENS');symbols=np.array(sorted(set(sy[valid])))
        assert len(symbols)>5000,(dataset,'insufficient mapped genes',len(symbols))
        projection=sparse.csr_matrix((np.ones(valid.sum()),(np.flatnonzero(valid),pd.Index(symbols).get_indexer(sy[valid]))),shape=(a.n_vars,len(symbols)))
        x=(raw@projection).tocsr();x.eliminate_zeros()
        m=pd.DataFrame(dict(cell=a.obs_names.astype(str),sample=a.obs['sample'].astype(str).to_numpy(),library=lib))
        m['day']=-1
        if dataset=='GSE314943':
            m['subtype']=a.obs['Subtype' if lineage=='CD4' else 'MacType'].astype(str).to_numpy();m['group']=a.obs['group'].astype(str).to_numpy()
            assert set(m.group)=={'Saline','Silica'} and m['sample'].nunique()==6
            assert (m.group==m['sample'].str.split('_').str[0]).all()
        elif dataset=='GSE183682':
            mp={'GSM5567501':('NS',7),'GSM5567502':('Silica',7),'GSM5567503':('NS',56),'GSM5567504':('Silica',56)}
            m['group']=m['sample'].map({k:v[0] for k,v in mp.items()});m['day']=m['sample'].map({k:v[1] for k,v in mp.items()})
            m['subtype']=a.obs.F08.astype(str).map({'0':'NaiveCD4','3':'ActCD4','5':'Treg'}).to_numpy()
        elif dataset=='GSE274372':
            m['day']=pd.to_numeric(a.obs.day).to_numpy();m['group']='D'+m.day.astype(int).astype(str)
            m['subtype']=a.obs.R08.astype(str).map({'0':'ActReg','1':'TLike1'}).to_numpy()
        else:
            m['group']=m['sample'].map({'GSM5324866':'Exposed','GSM5324867':'Exposed','GSM5324868':'Silicosis','GSM5324869':'Silicosis','GSM5324870':'Silicosis'})
            m['subtype']=a.obs.F08.astype(str).map({'0':'MemCD4','1':'NaiveCD4','12':'Treg'}).to_numpy()
        take=m.subtype.notna().to_numpy();m=m.loc[take].reset_index(drop=True);x=x[take];lib=lib[take]
        assert m[['sample','group','subtype','day']].notna().all().all()
        if dataset=='GSE274372':assert len(m)==494
        if dataset=='GSE174725':assert len(m)==2958
        mapped=np.asarray(x.sum(axis=1)).ravel()
        inventory.append(dict(dataset=dataset,lineage=lineage,path=str(base/rel),input_cells=a.n_obs,retained_cells=len(m),features=a.n_vars,mapped_symbols=len(symbols),mapped_UMI_fraction=mapped.sum()/lib.sum(),samples=m['sample'].nunique()))
        for gi,orig,cs,keep in zip(ids,original,sy,valid):mapping_all.append(dict(dataset=dataset,lineage=lineage,gene_id=gi,original=orig,canonical=cs,mapped=keep))
        m['bg']=np.asarray(x[:,bgmask(symbols)].sum(axis=1)).ravel()/lib
        m['mt']=np.asarray(x[:,np.char.startswith(symbols,'MT-')].sum(axis=1)).ravel()/lib
        m['nfeature']=np.asarray((x>0).sum(axis=1)).ravel()
        for g in set(TARGETS+SOURCE_GENES+RECEPTORS):
            ix=np.flatnonzero(symbols==g);m['count_'+g]=np.asarray(x[:,ix].sum(axis=1)).ravel() if len(ix) else np.nan
        # Save aligned raw counts and metadata on server only for the subsequent ligand-target step.
        if dataset=='GSE314943':
            if lineage=='CD4':main_symbols=symbols.copy()
            else:assert np.array_equal(main_symbols,symbols),'Main Mac/CD4 gene universe mismatch'
            prefix=out/'comm_input'/lineage
            with gzip.open(str(prefix)+'_counts.mtx.gz','wb') as handle:mmwrite(handle,x.T.astype(np.int64))
            pd.DataFrame({'canonical':symbols,'mouse_symbol':[mouse_names.get(g,g) for g in symbols]}).to_csv(str(prefix)+'_genes.tsv',sep='\t',index=False)
            m.assign(lineage=lineage)[['cell','sample','group','lineage','subtype','library']].to_csv(str(prefix)+'_cells.tsv',sep='\t',index=False)
            exports.append(dict(lineage=lineage,cells=len(m),genes=len(symbols),counts=str(prefix)+'_counts.mtx.gz'))
            for (ss,ty),v in m.groupby(['sample','subtype'],observed=True):paircounts.append(dict(sample=ss,group=v.group.iloc[0],lineage=lineage,subtype=ty,cells=len(v),eligible20=len(v)>=20))
        if lineage=='Mac':
            for ty in ['All']+sorted(m.subtype.unique()):
                z=m if ty=='All' else m[m.subtype==ty]
                for ss,v in z.groupby('sample',observed=True):
                    for g in SOURCE_GENES+TARGETS:
                        q=v['count_'+g].to_numpy();source.append(dict(sample=ss,group=v.group.iloc[0],subtype=ty,gene=g,cells=len(v),positive=np.sum(q>0) if np.isfinite(q).any() else np.nan,pct=100*np.mean(q>0) if np.isfinite(q).any() else np.nan,CP10k=np.mean(q/v.library.to_numpy()*1e4)))
            del a,raw,x,m
            continue
        y=x.multiply((1e4/lib)[:,None]).tocsr();y.data=np.log1p(y.data)
        values,members,cov,ctrl=scores(y,x,symbols,m['sample'].to_numpy())
        coverage_all.extend([dict(dataset=dataset,**v) for v in cov]);control_all.extend([dict(dataset=dataset,**v) for v in ctrl])
        for mod in SETS:m['score_'+mod]=values[mod];m['detect_'+mod]=members[mod]
        frames[dataset]=m
        m.to_csv(out/'objects'/f'{dataset}_cell_scores.tsv.gz',sep='\t',index=False,na_rep='NA')
        for filt,cut in [('All',1.0),('BG1',.01)]:
            z=m[m.bg<=cut]
            for scope in ['All']+sorted(m.subtype.unique()):
                v=z if scope=='All' else z[z.subtype==scope]
                for ss in sorted(m['sample'].unique()):
                    group=m.loc[m['sample']==ss,'group'].iloc[0];day=int(m.loc[m['sample']==ss,'day'].iloc[0]);w=v[v['sample']==ss];n=len(w)
                    original_n=len(m[(m['sample']==ss)&((m.subtype==scope) if scope!='All' else True)])
                    common=dict(dataset=dataset,scope=scope,filter=filt,sample=ss,group=group,day=day,cells=n,original_cells=original_n,eligible20=n>=20)
                    stats={'bg':w.bg.mean(),'mt':w.mt.mean(),'library':w.library.mean(),'nfeature':w.nfeature.mean()}
                    for mod in SETS:stats['score_'+mod]=w['score_'+mod].mean();stats['detect_'+mod]=w['detect_'+mod].mean()
                    for g in sorted(set(TARGETS+RECEPTORS)):
                        q=w['count_'+g].to_numpy();ok=n>0 and np.isfinite(q).all()
                        stats['expr_'+g]=np.mean(q/w.library.to_numpy()*1e4) if ok else np.nan
                        stats['pct_'+g]=100*np.mean(q>0) if ok else np.nan
                    sample_rows.append(dict(**common,**stats))
                    # Within-sample, within-subtype only; descriptive partial rank correlations.
                    if scope!='All':
                        for g in TARGETS:
                            q=w['count_'+g].to_numpy();npos=int(np.sum(q>0));eligible=n>=50 and npos>=10 and n-npos>=10
                            for mod in SETS:
                                yy=w['score_'+mod].to_numpy();ok=eligible and np.isfinite(q).all() and np.isfinite(yy).all()
                                within.append(dict(**common,gene=g,module=mod,positive=npos,eligible=bool(ok),rho=rho(q/w.library.to_numpy(),yy) if ok else np.nan,rho_depth_bg=partial(q/w.library.to_numpy(),yy,w[['library','bg']].to_numpy()) if ok else np.nan))
        print('DATASET_DONE='+dataset+'_'+lineage,flush=True)
        del a,raw,x,y
    tab=pd.DataFrame(sample_rows);assert not tab.duplicated(['dataset','scope','filter','sample']).any()
    wr(pd.DataFrame(inventory),'input_inventory.tsv');wr(pd.DataFrame(mapping_all),'gene_mapping.tsv.gz');wr(pd.DataFrame(coverage_all),'module_coverage.tsv');wr(pd.DataFrame(control_all),'score_controls.tsv.gz')
    wr(tab,'sample_scores.tsv');wr(pd.DataFrame(within),'within_sample_correlations.tsv');wr(pd.DataFrame(source),'mac_source_expression.tsv');wr(pd.DataFrame(paircounts),'communication_eligibility.tsv')
    wr(pd.DataFrame([dict(module=k,gene=g) for k,gs in SETS.items() for g in gs]),'module_definition.tsv')
    effects=[];assoc=[];standard=[]
    for (ds,scope,filt),v in tab.groupby(['dataset','scope','filter'],observed=True):
        pairs=[('Saline','Silica',-1)] if ds=='GSE314943' else [('NS','Silica',7),('NS','Silica',56)] if ds=='GSE183682' else [('D0','D7',-1),('D0','D28',-1),('D0','D56',-1)] if ds=='GSE274372' else [('Exposed','Silicosis',-1)]
        for c0,c1,day in pairs:
            w=v[v.group.isin([c0,c1])&(v.day.eq(day) if day!=-1 else True)]
            for minimum in [1,20]:
                z=w[w.cells>=minimum]
                for metric in ['score_'+k for k in SETS]+['expr_'+g for g in TARGETS]+['pct_'+g for g in TARGETS]:
                    d=effect(z.loc[z.group==c0,metric],z.loc[z.group==c1,metric]);effects.append(dict(dataset=ds,scope=scope,filter=filt,control=c0,case=c1,day=day,min_cells=minimum,metric=metric,**d))
                for g in TARGETS:
                    for mod in SETS:
                        q=z[['group','expr_'+g,'score_'+mod]].dropna();ok=len(q)>=6 and q.groupby('group',observed=True).size().min()>=2
                        assoc.append(dict(dataset=ds,scope=scope,filter=filt,control=c0,case=c1,day=day,min_cells=minimum,gene=g,module=mod,n_samples=len(q),eligible=bool(ok),rho_pooled=rho(q['expr_'+g],q['score_'+mod]) if ok else np.nan,rho_group=partial(q['expr_'+g],q['score_'+mod],(q.group==c1).astype(float)) if ok else np.nan))
    ef=pd.DataFrame(effects);wr(ef,'group_effects.tsv');wr(pd.DataFrame(assoc),'sample_associations.tsv')
    # Fixed baseline subtype proportions, calculated with equal sample weights.
    for ds in tab.dataset.unique():
        baseline=tab[(tab.dataset==ds)&(tab['filter']=='All')&(tab.scope!='All')]
        counts=baseline.pivot(index='sample',columns='scope',values='cells');weights=counts.div(counts.sum(axis=1),axis=0).mean()
        for (filt,ss),v in tab[(tab.dataset==ds)&(tab.scope!='All')].groupby(['filter','sample'],observed=True):
            z=v.set_index('scope').reindex(weights.index)
            for metric in ['score_'+k for k in SETS]+['expr_'+g for g in TARGETS]:
                standard.append(dict(dataset=ds,filter=filt,sample=ss,group=v.group.iloc[0],metric=metric,value=float(z[metric]@weights) if z[metric].notna().all() else np.nan,all_subtypes_eligible20=bool((z.cells>=20).all()),weights=json.dumps(weights.to_dict())))
    wr(pd.DataFrame(standard),'fixed_composition.tsv')
    # Descriptive matched-sample source/receiver associations, no invented LR edges.
    cross=[];src=pd.DataFrame(source)
    for scope in sorted(tab.loc[tab.dataset=='GSE314943','scope'].unique()):
        t=tab[(tab.dataset=='GSE314943')&(tab.scope==scope)&(tab['filter']=='All')&(tab.cells>=20)]
        for ty in sorted(src.subtype.unique()):
            for gene in SOURCE_GENES:
                s=src[(src.subtype==ty)&(src.gene==gene)&(src.cells>=20)][['sample','group','CP10k']]
                z=t.merge(s,on=['sample','group'],how='inner',validate='one_to_one')
                for metric in ['score_'+k for k in SETS]+['expr_'+g for g in TARGETS]:
                    q=z[['group','CP10k',metric]].dropna();ok=len(q)>=6 and q.group.nunique()==2
                    cross.append(dict(mac_type=ty,CD4_type=scope,source_gene=gene,receiver_metric=metric,n_samples=len(q),rho=rho(q.CP10k,q[metric]) if ok else np.nan,rho_group=partial(q.CP10k,q[metric],(q.group=='Silica').astype(float)) if ok else np.nan))
    wr(pd.DataFrame(cross),'mac_CD4_associations.tsv')
    human=root/'06_CD4T_analysis/CD4_final/GSE288887_CD4_final.rds'
    wr(pd.DataFrame([dict(dataset='HumanLung',path=str(human),exists=human.is_file(),size_bytes=human.stat().st_size if human.is_file() else 0,status='NEXT_STAGE_SOURCE_ADJUSTED_EXPORT')]),'human_lung_next_input.tsv')
    print('ANALYSIS_LOGIC_PASS=COUNTS_SCOPE_SAMPLE_UNIT_TARGET_EXCLUSION_FIXED_CONTROLS',flush=True)
    # English labels <=10 characters, arbitrary gene score units; cohorts are never pooled.
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(3,4,figsize=(11,7.5));v=tab[(tab.dataset=='GSE314943')&(tab.scope=='All')&(tab['filter']=='All')]
    for ax,mod in zip(axs.flat,SETS):
        for j,gr in enumerate(['Saline','Silica']):
            y=v.loc[(v.group==gr)&v.eligible20,'score_'+mod].dropna().to_numpy()
            ax.scatter(j+np.linspace(-.06,.06,len(y)),y,color=['#176B9B','#C44232'][j],s=26)
            if len(y):ax.hlines(y.mean(),j-.12,j+.12,color='black',lw=1)
        ax.set(title=mod,xticks=[0,1],xticklabels=['Saline','Silica'],ylabel='Score')
    fig.tight_layout();fig.savefig(out/'figures/module_scores.png',dpi=300);plt.close(fig)
    q=pd.DataFrame(assoc);q=q[(q.dataset=='GSE314943')&(q.scope=='All')&(q['filter']=='All')&(q.min_cells==20)]
    mat=q.pivot(index='gene',columns='module',values='rho_group').reindex(index=TARGETS,columns=list(SETS))
    fig,ax=plt.subplots(figsize=(10,4));im=ax.imshow(np.ma.masked_invalid(mat.to_numpy()),vmin=-1,vmax=1,cmap='RdBu_r',aspect='auto');ax.set(xticks=range(len(SETS)),xticklabels=list(SETS),yticks=range(len(TARGETS)),yticklabels=TARGETS);ax.tick_params(axis='x',rotation=45)
    for i in range(len(TARGETS)):
        for j in range(len(SETS)):
            if np.isfinite(mat.iloc[i,j]):ax.text(j,i,f'{mat.iloc[i,j]:.2f}',ha='center',va='center',fontsize=8)
    fig.colorbar(im,ax=ax,label='Adj rho');fig.tight_layout();fig.savefig(out/'figures/candidate_modules.png',dpi=300);plt.close(fig)
    (out/'methods.txt').write_text('Exploratory stage1. Modules are frozen exploratory panels: MitoISR/MitoROS/Mitophagy/PyroExec follow previous project definitions (14_GSE288887_CD4_damage_pyroptosis_modules.R); other lists are explicit manually defined transcript panels, not official pathway sets. No candidate belongs to a signature.\nIndependent samples are the reporting unit; no cell-level P values. Small cohorts: effect sizes, direction and leave-one-sample-out mean differences; no new P/FDR in this stage.\nScores: mean log1p(CP10k) of fixed signature genes minus mean matched control gene expression. This is a custom matched-expression score, not UCell or pathway activation. Controls use equal-sample baseline mean expression, 25 quantile bins and up to20 nearest genes per signature gene; all targets/signatures, background, mitochondrial-encoded and ribosomal genes excluded from control pool. Same controls apply to baseline/BG1. Absolute scores are not compared across cohorts.\nCoverage requires >=3 measured and >=3 expressed genes, >=60% measured membership, >=20 control genes; otherwise NA. Coverage/detection exported. OXPHOS uses nuclear genes. MitoROS is antioxidant response; PyroExec is transcriptional expression, not protein cleavage. Gene overlap between modules is retained in module_definition.tsv.\nRaw duplicate symbols summed; aliases explicitly canonicalized; GSE274372 mouse Ensembl mapping comes from the existing main mouse object. Full original raw-count library remains the CP10k denominator.\nMain analysis All; BG1 is fixed panel-UMI fraction<=1% sensitivity on CD4 only. Core and subtype rows separate. Min_cells1 descriptive; min_cells20 eligibility. Both original and retained sample counts exported. Fixed-composition weights use original equal-sample proportions.\nSample associations require>=6 samples and>=2 per group; adjusted rho residualizes rank-transformed expression/score on condition. Within-sample subtype correlations require>=50 cells,>=10 positive and>=10 negative cells; descriptive rho adjusts depth and panel burden. No correlation establishes direction or mediation.\nGSE183682 has one library per time/condition; GSE274372 is TLike with baseline comparisons; GSE174725 is human BALF with exposed controls; human lung RDS is inventoried for next source-adjusted stage, not analyzed here.\nMac source expression and matched-sample associations are screening summaries, not ligand activity or a curated ligand-receptor prediction. Raw expression for all mapped genes is exported to comm_input on server. NicheNet/MultiNicheNet and human lung fit remain next stages.\n',encoding='utf-8')
    (out/'captions.txt').write_text('module_scores.png: x=Saline/Silica; y=custom matched-expression module score; one point per sample, bars=equal-sample means; main Core scope, no background filtering. OXPHOS=oxidative phosphorylation; MitoISR=integrated mitochondrial stress response; MitoROS=antioxidant response; Mitophagy=mitophagy-associated genes; Inflam=inflammasome-associated genes; PyroExec=pyroptosis-associated execution transcripts; ApoExec=apoptosis-associated transcripts; TAct=T-cell activation; TReg=regulatory T-associated genes; Inhib=inhibitory-receptor-associated program; Naive=naive-associated genes; Cycle=cell-cycle genes.\ncandidate_modules.png: x=module; y=candidate; color/number=sample-level partial rank correlation controlling condition (Adj rho). No significance stars. Missing estimates are blank. MRE11 includes MRE11A alias. All labels are English and <=10 letters.\n')
    status=dict(status='PASS',stage='mechanism_step1',python=platform.python_version(),anndata=ad.__version__,inputs=inventory,exports=exports,pending=['human_lung_source_adjusted_scoring','ligand_target_inference','candidate_regulatory_network'])
    (out/'status.json').write_text(json.dumps(status,indent=2))
    # Only figures, compact tables and executable provenance are returned; large raw matrices remain on server.
    dest=root/'tmp/lung_mechanism_step1_return.tar.gz';tmp=dest.with_name(dest.name+'.part')
    with tarfile.open(tmp,'w:gz') as tf:
        for folder in ['tables','figures']:
            for p in sorted((out/folder).iterdir()):tf.add(p,arcname=str(p.relative_to(out)))
        for name in ['analyze.py','methods.txt','captions.txt','status.json']:tf.add(out/name,arcname=name)
    assert tmp.stat().st_size<25*1024*1024,'Return package exceeds25MiB; outputs remain on server'
    os.replace(tmp,dest)
    (root/'tmp/lung_mechanism_step1_latest.txt').write_text(str(out)+'\n')
    print('RETURN_PACKAGE='+str(dest),flush=True);print('RETURN_BYTES='+str(dest.stat().st_size),flush=True)
    print('LUNG_MECHANISM_STEP1_PASS',flush=True)

if __name__=='__main__':main()

