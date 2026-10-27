import os,json,gzip,tarfile,hashlib,sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from scipy.stats import rankdata

TARGETS=['IL1R2','PENK','TTC39C','FOXP1','TNFRSF9','LAG3','MRE11']
# Species-specific mouse Il6ra is a symbol alias for the IL6R comparison here.
ALIASES={'MRE11A':'MRE11','IL6RA':'IL6R'}
BGPANEL=set('SFTPA1 SFTPA2 SFTPB SFTPC SFTPD SCGB1A1 SCGB3A1 SCGB3A2 ABCA3 SLC34A2 NAPSA LYZ LYZ1 LYZ2 S100A8 S100A9 CSF1R LST1 TYROBP FCER1G'.split())
AXES={
 'IL15_TP':(['IL15','IL15RA'],['IL2RB','IL2RG'],'IL15/IL15RA donor and IL2RB/IL2RG receiver; RNA co-detection does not establish surface trans-presentation'),
 'TNF_R1':(['TNF'],['TNFRSF1A'],'TNF receptor1'),
 'TNF_R2':(['TNF'],['TNFRSF1B'],'TNF receptor2'),
 'TGFB':(['TGFB1'],['TGFBR1','TGFBR2'],'Latent ligand activation is not measured'),
 'IL1_SIG':(['IL1B'],['IL1R1','IL1RAP'],'Signaling receptor complex; cytokine maturation is not measured'),
 'IL1_DECOY':(['IL1B'],['IL1R2'],'Decoy binding; not a signaling receptor axis'),
 'CD86_CD28':(['CD86'],['CD28'],'Costimulation; ligand-target prior is not evidence of a direct path'),
 'CD86_CTLA4':(['CD86'],['CTLA4'],'Inhibitory receptor context'),
 'IL18':(['IL18'],['IL18R1','IL18RAP'],'Cytokine maturation is not measured'),
 'CXCL16':(['CXCL16'],['CXCR6'],'Membrane/soluble forms not distinguished'),
 'VEGFA_KDR':(['VEGFA'],['KDR','NRP1'],'Canonical receptor plus co-receptor audit'),
 'VEGFA_FLT':(['VEGFA'],['FLT1','NRP1'],'Canonical receptor plus co-receptor audit')}

def canonical(s):
    s=str(s).strip().upper();return ALIASES.get(s,s)

def load(countpath,genepath,metapath,column):
    with (gzip.open(countpath,'rb') if str(countpath).endswith('.gz') else open(countpath,'rb')) as f:x=sparse.csr_matrix(mmread(f).T,dtype=float)
    g=pd.read_csv(genepath,sep='\t');m=pd.read_csv(metapath,sep='\t')
    assert x.shape==(len(m),len(g)) and m.cell.is_unique
    assert np.isfinite(x.data).all() and (x.data>=0).all() and np.allclose(x.data,np.rint(x.data),rtol=0,atol=1e-6)
    assert (m.library>0).all()
    names=np.array([canonical(s) for s in g[column].fillna('')]);keep=(names!='')&~np.char.startswith(names,'ENS')
    unique=np.array(sorted(set(names[keep])))
    pr=sparse.csr_matrix((np.ones(keep.sum()),(np.flatnonzero(keep),pd.Index(unique).get_indexer(names[keep]))),shape=(len(names),len(unique)))
    x=(x@pr).tocsr();x.eliminate_zeros()
    assert (np.asarray(x.sum(axis=1)).ravel()<=m.library.to_numpy()+1e-6).all()
    import re
    b=np.array([s in BGPANEL or bool(re.match(r'^(IGH[VDJMAEG]|IGK[VJC]|IGL[VJC])',s)) for s in unique])
    m['bg']=np.asarray(x[:,b].sum(axis=1)).ravel()/m.library.to_numpy()
    return x,unique,m,pd.DataFrame({'input':g[column],'canonical':names,'retained':keep})

def partial(x,y,group):
    x=np.asarray(x,float);y=np.asarray(y,float);group=np.asarray(group)
    ok=np.isfinite(x)&np.isfinite(y);x=x[ok];y=y[ok];group=group[ok]
    if len(x)<5 or len(set(group))!=2 or min(pd.Series(group).value_counts())<2:return np.nan
    z=np.column_stack([np.ones(len(x)),group==sorted(set(group))[1]])
    u=rankdata(x);v=rankdata(y);u=u-z@np.linalg.lstsq(z,u,rcond=None)[0];v=v-z@np.linalg.lstsq(z,v,rcond=None)[0]
    if np.std(u)<1e-10 or np.std(v)<1e-10:return np.nan
    return float(np.clip(np.corrcoef(u,v)[0,1],-1,1))

def audit_scope(x,genes,m,dataset,scope,filt,lineage):
    # x and m are already aligned and filtered for one scope.
    expected=sorted(set(TARGETS)|set(g for v in AXES.values() for gs in v[:2] for g in gs))
    lookup={g:i for i,g in enumerate(genes)}
    vals={g:np.asarray(x[:,lookup[g]].toarray()).ravel() if g in lookup else None for g in expected}
    comp=[];single=[]
    for sample,indices in m.groupby('sample',observed=True,sort=True).indices.items():
        q=np.asarray(indices);md=m.iloc[q];n=len(q)
        base=dict(dataset=dataset,lineage=lineage,scope=scope,filter=filt,sample=sample,group=md.group.iloc[0],series=md.series.iloc[0] if 'series' in md else dataset,cells=n,eligible20=n>=20)
        for g,v in vals.items():
            present=v is not None
            a=v[q] if present else None
            single.append(dict(**base,gene=g,present=present,positive=int((a>0).sum()) if present else np.nan,pct=100*float(np.mean(a>0)) if present else np.nan,mean_CP10k=float(np.mean(a/md.library.to_numpy()*1e4)) if present else np.nan))
        for axis,(sender,receiver,note) in AXES.items():
            members=sender if lineage=='Mac' else receiver
            present=all(vals[g] is not None for g in members)
            co=np.all(np.column_stack([vals[g][q]>0 for g in members]),axis=1) if present else None
            row=dict(**base,axis=axis,members=';'.join(members),all_features_present=present,missing=';'.join(g for g in members if vals[g] is None),codetected=int(co.sum()) if present else np.nan,pct=100*float(co.mean()) if present else np.nan,interpretation=note)
            comp.append(row)
    return comp,single

def main():
    root=Path(os.environ['ROOT']);out=Path(os.environ['OUT']);p1=Path(os.environ['PREV1']);p2=Path(os.environ['PREV2']);p3=Path(os.environ['PREV3'])
    for p,stage in [(p1,'mechanism_step1'),(p2,'mechanism_step2_human_and_driver_audit'),(p3,'mechanism_step3_nichenet')]:
        s=json.loads((p/'status.json').read_text());assert s['status']=='PASS' and s['stage']==stage
    def wr(df,n):df.to_csv(out/'tables'/n,sep='\t',index=False,na_rep='NA')
    comp=[];single=[];inventory=[]
    for lineage,expected in [('CD4',2311),('Mac',1796),('HumanCD4',6523)]:
        if lineage=='HumanCD4':
            x,g,m,mapping=load(p2/'input/human_counts.mtx',p2/'input/human_genes.tsv',p2/'input/human_cells.tsv','symbol')
            dataset='HumanLung';line='CD4'
            core=m.core.astype(str).str.lower().map({'true':True,'false':False,'1':True,'0':False});assert core.notna().all();m['core']=core.astype(bool)
        else:
            x,g,m,mapping=load(p1/'comm_input'/f'{lineage}_counts.mtx.gz',p1/'comm_input'/f'{lineage}_genes.tsv',p1/'comm_input'/f'{lineage}_cells.tsv','canonical')
            dataset='GSE314943';line=lineage;m['core']=True
        assert len(m)==expected,(lineage,len(m))
        wr(mapping,f'{lineage}_mapping.tsv')
        scopes=['Core']+sorted(m.loc[m.core,'subtype'].unique()) if line=='CD4' else ['AllMac']+sorted(m.subtype.unique())
        for scope in scopes:
            mask=m.core.to_numpy() if scope=='Core' else np.ones(len(m),bool) if scope=='AllMac' else (m.subtype==scope).to_numpy()
            for filt in (['All','BG1'] if line=='CD4' else ['All']):
                take=mask & ((m.bg<=.01).to_numpy() if filt=='BG1' else True)
                v=m.loc[take].reset_index(drop=True)
                co,ss=audit_scope(x[take],g,v,dataset,scope,filt,line)
                comp.extend(co);single.extend(ss)
                for samp in sorted(m['sample'].unique()):
                    inventory.append(dict(dataset=dataset,lineage=line,scope=scope,filter=filt,sample=samp,cells=int((v['sample']==samp).sum())))
        print('RAW_COMPLEX_AUDIT_DONE='+dataset+'_'+line,flush=True)
    c=pd.DataFrame(comp);e=pd.DataFrame(single)
    wr(c,'complex_codetection.tsv');wr(e,'axis_gene_expression.tsv');wr(pd.DataFrame(inventory),'cell_inventory.tsv')
    summary=[]
    for key,v in c.groupby(['dataset','lineage','scope','filter','axis','group'],observed=True):
        z=v[v.eligible20];d=dict(zip(['dataset','lineage','scope','filter','axis','group'],key))
        d.update(samples=len(v),eligible_samples=len(z),features_present=bool(v.all_features_present.all()),mean_pct=z.pct.mean(),min_pct=z.pct.min(),max_pct=z.pct.max(),positive_cells=z.codetected.sum(min_count=1))
        for cut in [1,5,10]:d['samples_ge'+str(cut)+'pct']=int((z.pct>=cut).sum())
        summary.append(d)
    wr(pd.DataFrame(summary),'complex_group_summary.tsv')
    # Match sample-specific receiver detection and sender detection. No product called probability of communication.
    pairs=[]
    mc=c[(c.dataset=='GSE314943')&(c.lineage=='Mac')];tc=c[(c.dataset=='GSE314943')&(c.lineage=='CD4')]
    for (sender,axis),v in mc.groupby(['scope','axis'],observed=True):
        for (receiver,filt),w in tc[tc.axis==axis].groupby(['scope','filter'],observed=True):
            z=v.merge(w,on=['sample','group','axis'],suffixes=('_sender','_receiver'),validate='one_to_one')
            for _,r in z.iterrows():
                pairs.append(dict(sender=sender,receiver=receiver,filter=filt,axis=axis,sample=r['sample'],group=r['group'],sender_cells=r.cells_sender,receiver_cells=r.cells_receiver,sender_pct=r.pct_sender,receiver_pct=r.pct_receiver,eligible20=bool(r.eligible20_sender and r.eligible20_receiver),all_features_present=bool(r.all_features_present_sender and r.all_features_present_receiver)))
    pair=pd.DataFrame(pairs);wr(pair,'paired_complex_support.tsv')
    # ActCD4 Il1r2/Penk can fail the 10% target-expression gate despite real DE: report thresholds, no relabeling of the original model.
    gates=[]
    for key,v in e[(e.dataset=='GSE314943')&(e.lineage=='CD4')&(e.gene.isin(TARGETS))].groupby(['scope','filter','gene'],observed=True):
        for cut in [1,5,10]:
            n=[int(((v.group==gr)&v.eligible20&(v.pct>=cut)).sum()) for gr in ['Saline','Silica']]
            gates.append(dict(scope=key[0],filter=key[1],gene=key[2],cutoff_pct=cut,control_supported=n[0],case_supported=n[1],either_group_ge2=max(n)>=2))
    wr(pd.DataFrame(gates),'candidate_expression_gate_audit.tsv')
    # Matched-sample correlations now include ActCD4 exploratory candidates too, not only FDR NicheNet sets.
    assoc=[]
    sources={'IL15_TP':'IL15','TNF_R1':'TNF','TNF_R2':'TNF','TGFB':'TGFB1','IL1_SIG':'IL1B','IL1_DECOY':'IL1B','CD86_CD28':'CD86','CD86_CTLA4':'CD86','IL18':'IL18','CXCL16':'CXCL16','VEGFA_KDR':'VEGFA','VEGFA_FLT':'VEGFA'}
    me=e[(e.dataset=='GSE314943')&(e.lineage=='Mac')];te=e[(e.dataset=='GSE314943')&(e.lineage=='CD4')]
    for (sender,receiver,filt,axis),v in pair.groupby(['sender','receiver','filter','axis'],observed=True):
        v=v[v.eligible20]
        ligand=sources[axis];lv=me[(me.scope==sender)&(me.gene==ligand)][['sample','mean_CP10k']].rename(columns={'mean_CP10k':'ligand_expr'})
        for gene in TARGETS:
            tv=te[(te.scope==receiver)&(te['filter']==filt)&(te.gene==gene)][['sample','mean_CP10k']].rename(columns={'mean_CP10k':'target_expr'})
            z=v.merge(lv,on='sample',validate='one_to_one').merge(tv,on='sample',validate='one_to_one')
            for metric in ['ligand_expr','sender_pct','receiver_pct']:
                zz=z[np.isfinite(z[metric])&np.isfinite(z.target_expr)].copy()
                nc=int((zz.group=='Saline').sum());nd=int((zz.group=='Silica').sum())
                r=partial(zz[metric],zz.target_expr,zz.group) if nc>=3 and nd>=3 else np.nan
                ll=[partial(zz.drop(i)[metric],zz.drop(i).target_expr,zz.drop(i).group) for i in zz.index] if np.isfinite(r) else []
                ll=np.asarray(ll);ll=ll[np.isfinite(ll)]
                assoc.append(dict(sender=sender,receiver=receiver,filter=filt,axis=axis,ligand=ligand,target=gene,metric=metric,n_control=nc,n_case=nd,rho=r,loo_min=float(ll.min()) if len(ll) else np.nan,loo_max=float(ll.max()) if len(ll) else np.nan,loo_same=float(np.mean(np.sign(ll)==np.sign(r))) if len(ll) and abs(r)>1e-10 else np.nan,status='DESCRIPTIVE_3V3' if nc>=3 and nd>=3 else 'INSUFFICIENT_3V3'))
    wr(pd.DataFrame(assoc),'paired_axis_candidate_correlations.tsv')
    # Save original method outcomes and their gene-set denominators with the receptor audit.
    for name in ['geneset_status.tsv','candidate_prior_coverage.tsv','ligand_rank_sensitivity.tsv']:
        wr(pd.read_csv(p3/'tables'/name,sep='\t'),'step3_'+name)
    methods='''Step4 audits receptor biology and sample coverage using frozen raw-count exports. No new cohort, embedding, or cluster selection. Independent sample is the summary unit. Core and each original subtype, All and <=1% CD4 background burden are separate; macrophages are not filtered by a myeloid background panel.
Canonical mapping: mouse IL6RA to IL6R; MRE11A to MRE11; missing features are NA, measured zero remains zero. The original full raw-library total is the expression denominator. Complex detection requires nonzero counts for every listed gene in the same cell. Thresholds 1%,5%,10% are reported together; no threshold is selected for significance. These are RNA expression feasibility checks, not evidence of protein complex assembly or spatial contact.
IL15 trans-presentation requires ligand/chaperone context on the donor and beta/gamma receptor context on the receiver; IL2RA alone is not accepted as the IL15 signaling receptor. IL1R2 is classified as a decoy. VEGFA co-receptor NRP1 alone is not called a complete canonical receptor complex.
Paired correlations are descriptive, adjust condition by residualized ranks, and require3 control and3 case samples. Leave-one-out requires at least2 per condition. They do not prove causality or pass a predeclared significance test. Disease-enriched Spp1Mac with no eligible control samples cannot provide a balanced DE comparison. Human lung provides receiver expression only, not matched human macrophage communication validation.
ActCD4 target-expression gate audit reports 1/5/10% support explicitly. It does not overwrite or relabel step3 FDR/nominal results. Downregulated gene-set activity in NicheNet does not assign inhibitory direction to a ligand.
IL15 primary reference: https://pubmed.ncbi.nlm.nih.gov/34654688/
IL1R2 curated mouse protein: https://www.uniprot.org/uniprotkb/P27931/entry
'''
    (out/'methods.txt').write_text(methods)
    (out/'status.json').write_text(json.dumps(dict(status='PASS',stage='mechanism_step4_receptor_audit',previous=str(p3),output=str(out),complex_rows=len(c),gene_rows=len(e),note='No new NicheNet ranking; targeted raw-count receptor feasibility and threshold audit.'),indent=2))
    dest=root/'tmp/lung_mechanism_step4_return.tar.gz';part=dest.with_suffix('.tar.gz.part')
    with tarfile.open(part,'w:gz') as t:
        for p in sorted(out.rglob('*')):
            if p.is_file() and '__pycache__' not in p.parts:t.add(p,arcname=str(p.relative_to(out)))
    with tarfile.open(part) as t:assert any(x.name=='status.json' for x in t)
    part.replace(dest)
    print('ANALYSIS_LOGIC_PASS RAW_COUNTS_COMPLEXES_SAMPLE_ALIGNMENT',flush=True)
    print('RETURN='+str(dest),flush=True);print('RETURN_BYTES='+str(dest.stat().st_size),flush=True)
    print('RETURN_SHA256='+hashlib.sha256(dest.read_bytes()).hexdigest(),flush=True)
    print('LUNG_MECHANISM_STEP4_PASS',flush=True)

if __name__=='__main__':main()

