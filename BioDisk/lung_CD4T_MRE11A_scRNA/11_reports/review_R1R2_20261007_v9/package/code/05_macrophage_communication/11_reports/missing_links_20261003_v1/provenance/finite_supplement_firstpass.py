"""Finite follow-up only. Inputs immutable; no GRN, edgeR or legacy model rerun."""
from pathlib import Path
import os, sys, json, hashlib, itertools, fcntl, datetime, re
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.io import mmread
import anndata as ad

P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
O=P/'11_reports/missing_links_20261003_v1'
A=P/'02_analysis'
M1=A/'mechanism_step1_20260927_152218_s4V5eF'
M2=A/'mechanism_step2_20260927_154919_EsN4gk'
S1=A/'scenic1_20260929_103513_xRoEnr'
S2=A/'scenic2_20261001_092655_ZKNWD6'
S3=A/'scenic3_20261001_123452_Kvp9h5'
RV=P/'11_reports/technical_delivery_20261002_ikzf2_stat1_v1/revision_v2/final_release/prior_review_repaired/scenic3_review/results'
lock=open(O/'analysis.lock','a+')
fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
lock.seek(0);lock.truncate();lock.write(str(os.getpid()));lock.flush()
if (O/'provenance/complete.json').exists():raise RuntimeError('Completed supplement exists; refusing rerun')
manifest={}
def track(p):
 p=Path(p)
 if str(p) not in manifest:
  h=hashlib.sha256()
  with p.open('rb') as f:
   for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
  manifest[str(p)]=dict(path=str(p),bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns,sha256=h.hexdigest())
 return p
def rd(p,**kw):return pd.read_csv(track(p),sep='\t',**kw)
def wr(d,name):pd.DataFrame(d).to_csv(O/'results'/name,sep='\t',index=False,na_rep='NA')
def canon(g):
 s=str(g).upper()
 return {'ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','MRE11A':'MRE11'}.get(s,s)
ligands='IL15 IL15RA TNF IL18 CXCL16 IL27 EBI3 IFNB1 IFNG'.split()
complexes={'IL15':['IL2RB','IL2RG'],'TNFR1':['TNFRSF1A'],'TNFR2':['TNFRSF1B'],'IL18':['IL18R1','IL18RAP'],'CXCL16':['CXCR6'],'IL27':['IL27RA','IL6ST'],'IFNB1':['IFNAR1','IFNAR2'],'IFNG':['IFNGR1','IFNGR2']}
modules={
 'STAT5_proxy':'CISH SOCS2 BCL2 IL2RA PIM1',
 'NFkB_proxy':'NFKBIA TNFAIP3 NFKBIZ',
 'RELB_components':'RELB NFKB2',
 'AP1_components':'JUN JUNB FOS FOSB ATF3',
 'IFN_state':'IFIT1 IFIT2 IFIT3 ISG15 MX1 OAS1 IRF1 GBP1',
 'TregNoTF':'FOXP3 IL2RA CTLA4 IKZF4 LRRC32 ENTPD1 TIGIT',
 'Activation':'CD69 IL2RA ICOS TNFRSF4 BATF IRF4 NR4A1 NR4A2 EGR1',
 'Cycle':'MKI67 TOP2A PCNA MCM2 MCM3 MCM4 MCM5 MCM6 MCM7 UBE2C'}
modules={k:v.split() for k,v in modules.items()}
oldmods=rd(M1/'tables/module_definition.tsv')
modules['OXPHOS_nuclear']=oldmods.loc[oldmods.module.eq('OXPHOS'),'gene'].map(canon).tolist()
assert len(modules['OXPHOS_nuclear'])==24 and not any(g.startswith('MT-') for g in modules['OXPHOS_nuclear'])
wr([dict(module=k,gene=g,definition='frozen original nuclear OXPHOS' if k=='OXPHOS_nuclear' else 'prespecified exploratory RNA proxy; not a validated activity assay') for k,gs in modules.items() for g in gs],'module_definition.tsv')
wr([dict(axis=k,receiver_components='+'.join(v),sender_components='IL15+IL15RA' if k=='IL15' else 'IL27+EBI3' if k=='IL27' else 'TNF' if k.startswith('TNFR') else k) for k,v in complexes.items()],'topology.tsv')
coverage=[];gene_audit=[];pair_audit=[];overlap=[];receiver=[];sender=[]

def load(prefix,gcol):
 prefix=str(prefix);fp=Path(prefix+'_counts.mtx.gz')
 if not fp.exists():fp=Path(prefix+'_counts.mtx')
 x=mmread(track(fp)).T.tocsr();g=rd(prefix+'_genes.tsv')[gcol].map(canon).to_numpy();m=rd(prefix+'_cells.tsv')
 assert x.shape==(len(m),len(g)) and m.cell.is_unique and (m.library>0).all()
 assert np.isfinite(x.data).all() and (x.data>=0).all() and np.allclose(x.data,np.rint(x.data))
 un,ix=np.unique(g,return_inverse=True)
 x=x@sparse.csr_matrix((np.ones(len(g)),(np.arange(len(g)),ix)),shape=(len(g),len(un)))
 return x.tocsr(),un,m

def counts(x,g,genes):
 return {v:np.asarray(x[:,np.flatnonzero(g==v)].sum(axis=1)).ravel() if np.any(g==v) else np.full(x.shape[0],np.nan) for v in genes}

def aggregate_sender(sp,source,m,val):
 m=m.copy().reset_index(drop=True)
 for sample,z in m.groupby('sample',sort=True):
  ids=z.index.to_numpy();lib=z.library.to_numpy();r=dict(species=sp,sender=source,sample=sample,group=z.group.iloc[0],sender_cells=len(z),sender_library=lib.mean())
  for gene in ligands:
   v=val[gene][ids];r['L_'+gene]=np.mean(v/lib*1e4)
   gene_audit.append(dict(species=sp,role='sender',scope=source,sample=sample,gene=gene,cells=len(z),present=np.isfinite(v).all(),positive=int(np.sum(v>0)) if np.isfinite(v).all() else np.nan,pct=100*np.mean(v>0) if np.isfinite(v).all() else np.nan,mean_CP10k=r['L_'+gene],raw_sum=v.sum()))
  for name,gs in [('IL15',['IL15','IL15RA']),('IL27',['IL27','EBI3'])]:
   a=np.column_stack([val[g][ids] for g in gs]);r['L_'+name+'_complex']=min(r['L_'+g] for g in gs) if np.isfinite(a).all() else np.nan
   r['joint_'+name]=100*np.mean((a>0).all(axis=1)) if np.isfinite(a).all() else np.nan
  sender.append(r)

def process_receiver(sp,x,g,m):
 m=m.copy().reset_index(drop=True)
 if sp=='Human':m=m[m.core.astype(str).str.lower().eq('true')].copy();x=x[m.index];m=m.reset_index(drop=True)
 else:m['series']='GSE314943'
 genes=sorted(set(ligands+['IKZF2','STAT1']+list(itertools.chain.from_iterable(complexes.values()))+list(itertools.chain.from_iterable(modules.values()))))
 raw=counts(x,g,genes);lib=m.library.to_numpy()
 # Reference scale is frozen over all original Core cells, before sample filtering or tests.
 expr={k:np.log1p(v/lib*1e4) for k,v in raw.items()}
 zval={k:(v-v.mean())/v.std() if np.isfinite(v).all() and v.std()>0 else np.full(len(v),np.nan) for k,v in expr.items()}
 edge=rd(S1/sp/'All/regulon_edges.tsv');edge.tf=edge.tf.map(canon);edge.target=edge.target.map(canon)
 ex=set(oldmods.gene.map(canon))|set('IL1R2 PENK TTC39C FOXP1 TNFRSF9 LAG3 MRE11'.split())
 union=set(edge.loc[edge.tf.isin(['STAT1','IKZF2']),'target'])
 for tf in ['STAT1','IKZF2']:
  for seed in [17,29,43]:
   original=set(edge.loc[edge.tf.eq(tf)&edge.seed.eq(seed),'target']);drop=original-ex-{tf}
   for name,gs in modules.items():
    overlap.append(dict(species=sp,tf=tf,seed=seed,module=name,original_overlap=';'.join(sorted(original&set(gs))),DropOverlap_overlap=';'.join(sorted(drop&set(gs))),module_size=len(gs),original_overlap_n=len(original&set(gs)),DropOverlap_overlap_n=len(drop&set(gs))))
    if name=='OXPHOS_nuclear':assert not drop&set(gs)
 scores={}
 for name,gs in modules.items():
  for variant,expected in [('full',gs),('NoRegUnion',[q for q in gs if q not in union])]:
   used=[q for q in expected if q in set(g) and np.isfinite(zval[q]).all()]
   eligible=len(used)>=2 and len(used)>=.6*len(expected)
   scores[name+'_'+variant]=np.mean([zval[q] for q in used],axis=0) if eligible else np.full(len(m),np.nan)
   coverage.append(dict(species=sp,module=name,variant=variant,original_expected=len(gs),after_overlap=len(expected),used=len(used),eligible=eligible,used_genes=';'.join(used),removed_regulon_genes=';'.join(sorted(set(gs)-set(expected))),missing_or_constant=';'.join(sorted(set(expected)-set(used)))))
 for scope in ['Core','CD4-IL7R' if sp=='Human' else 'ActCD4']:
  take=np.ones(len(m),bool) if scope=='Core' else m.subtype.eq(scope).to_numpy()
  for sample,v in m[take].groupby('sample',sort=True):
   ids=v.index.to_numpy();r=dict(species=sp,scope=scope,sample=sample,group=v.group.iloc[0],series=v.series.iloc[0],receiver_cells=len(v),receiver_library=v.library.mean(),Treg_fraction=float(v.subtype.str.contains('Treg',case=False).mean()))
   for tf in ['IKZF2','STAT1']:r[tf+'_RNA']=np.mean(raw[tf][ids]/lib[ids]*1e4)
   for name,score in scores.items():r[name]=score[ids].mean()
   r['IFNG_RNA']=np.mean(raw['IFNG'][ids]/lib[ids]*1e4)
   for ax,gs in complexes.items():
    a=np.column_stack([raw[q][ids] for q in gs]);r['R_'+ax]=100*np.mean((a>0).all(axis=1)) if np.isfinite(a).all() else np.nan
    for gene in gs:
     q=raw[gene][ids];gene_audit.append(dict(species=sp,role='receiver',scope=scope,sample=sample,gene=gene,cells=len(v),present=np.isfinite(q).all(),positive=int(np.sum(q>0)) if np.isfinite(q).all() else np.nan,pct=100*np.mean(q>0) if np.isfinite(q).all() else np.nan,mean_CP10k=np.mean(q/lib[ids]*1e4),raw_sum=q.sum()))
   receiver.append(r)
 return set(m.cell)

print('EXTRACT_RECEIVERS',flush=True)
mouse_x,mouse_g,mouse_m=load(M1/'comm_input/CD4','canonical')
mouse_cd4=process_receiver('Mouse',mouse_x,mouse_g,mouse_m)
human_x,human_g,human_m=load(M2/'input/human','symbol')
human_all_cd4=set(human_m.cell)
human_cd4=process_receiver('Human',human_x,human_g,human_m)
print('EXTRACT_SENDERS',flush=True)
x,g,m=load(M1/'comm_input/Mac','canonical');take=m.subtype.eq('AMlike');aggregate_sender('Mouse','AM',m[take],counts(x[take],g,ligands))
x,g,m=load(S3/'input/human_mac','gene');review=rd(RV/'reviewed_cells.tsv.gz');assert set(review.cell)==set(m.cell)
label=review.set_index('cell').review_scope.reindex(m.cell).to_numpy();take=(label=='MacAM')&m['sample'].notna().to_numpy()&m['sample'].ne('Unknown').to_numpy()
aggregate_sender('Human','AM',m[take],counts(x[take],g,ligands));del x
t=rd(O/'inputs/human_TNK_selected_counts.tsv');take=~t.cell.isin(human_all_cd4);tt=t[take].copy();aggregate_sender('Human','TNK_other',tt,{v:tt[v].to_numpy() for v in ligands})
p=track(A/'GSE314943/stage4_20260905_215942_3XPis4/objects/TNK_harmony.h5ad');a=ad.read_h5ad(p)
take=~a.obs_names.isin(mouse_cd4);a=a[take].copy();x=sparse.csr_matrix(a.layers['counts']);g=a.var.symbol.map(canon).to_numpy();m=pd.DataFrame(dict(cell=a.obs_names,sample=a.obs['sample'].astype(str).to_numpy(),group=a.obs.group.astype(str).to_numpy(),library=np.asarray(x.sum(axis=1)).ravel()))
aggregate_sender('Mouse','TNK_other',m,counts(x,g,ligands));del a,x
rec=pd.DataFrame(receiver);sen=pd.DataFrame(sender)
assert not rec.duplicated(['species','scope','sample']).any() and not sen.duplicated(['species','sender','sample']).any()
wr(rec,'receiver_raw_all_samples.tsv');wr(sen,'sender_all_samples.tsv');wr(gene_audit,'gene_coverage_by_sample.tsv');wr(coverage,'module_coverage.tsv');wr(overlap,'regulon_overlap_audit.tsv')
# Reuse only frozen AUC, never impute missing upstream excluded donor AUC.
old=rd(P/'11_reports/stepwise_diagnostic_20261002_v1/inputs/receiver_all_positive_cell_samples.tsv')
old=old[old['filter'].eq('All')]
wide=old.pivot(index=['species','scope','sample'],columns=['tf','seed'],values='mean_auc');wide.columns=[f'{tf}_AUC_{seed}' for tf,seed in wide.columns]
rec=rec.merge(wide.reset_index(),on=['species','scope','sample'],how='left',validate='one_to_one')
for row in rec.itertuples():
 pair_audit.append(dict(species=row.species,scope=row.scope,sample=row.sample,kind='RNA_vs_AUC',receiver_cells=row.receiver_cells,AUC_available=np.isfinite(row.IKZF2_AUC_17),reason='' if np.isfinite(row.IKZF2_AUC_17) else 'not in frozen upstream SCENIC roster; RNA retained, AUC not imputed'))
wr(rec,'receiver_paired_modalities.tsv')
paired=[]
for (sp,scope),r in rec.groupby(['species','scope']):
 for source in ['AM','TNK_other']:
  s=sen[sen.species.eq(sp)&sen.sender.eq(source)]
  q=r.merge(s.drop(columns='species'),on='sample',how='outer',suffixes=('','_sender'),indicator=True,validate='one_to_one')
  for _,v in q.iterrows():pair_audit.append(dict(species=sp,scope=scope,sample=v['sample'],kind=source,receiver_cells=v.receiver_cells,sender_cells=v.sender_cells,join=v['_merge'],reason='' if v['_merge']=='both' else str(v['_merge'])))
  q=q[q['_merge'].eq('both')].copy();assert q.group.eq(q.group_sender).all();paired.append(q.drop(columns=['_merge','group_sender']))
pair=pd.concat(paired,ignore_index=True);wr(pair,'sender_receiver_paired_inputs.tsv');wr(pair_audit,'pairing_audit.tsv')
print('COVERAGE_AND_PAIRING_WRITTEN',len(rec),len(sen),len(pair),flush=True)
(O/'provenance/extraction_complete.json').write_text(json.dumps(dict(receiver_rows=len(rec),sender_rows=len(sen),paired_rows=len(pair),time=datetime.datetime.now().isoformat()),indent=2))
pd.DataFrame(manifest.values()).to_csv(O/'provenance/input_manifest.tsv',sep='\t',index=False)

if '--extract-only' in sys.argv:sys.exit(0)
