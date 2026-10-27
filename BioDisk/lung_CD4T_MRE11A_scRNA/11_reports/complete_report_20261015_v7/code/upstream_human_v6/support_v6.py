from pathlib import Path
import sys,json,hashlib,itertools
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from scipy.stats import permutation_test,false_discovery_control
O=Path(sys.argv[1]);P=O.parents[1];V=P/'11_reports/integrated_macrophage_CD4_20261015_v4';H=P/'11_reports/GSE174725_annotation_validation_20261005_v1/results';G=['STAT1','GSDMD','CASP1','C3','C3AR1','CXCL16','CXCR6','IL27','EBI3','IL27RA','IL6ST','IFNG','IFNGR1','IFNGR2','IFNB1','IFNAR1','IFNAR2','IL12A','IL12B','IL12RB1','IL12RB2','TNF','TNFRSF1A','TNFRSF1B','IL15','IL15RA','IL2RB','IL2RG','IL18','IL18R1','IL18RAP'];pairs={'CXCL16_CXCR6':['CXCR6'],'C3a_C3AR1':['C3AR1'],'IL27_receptor':['IL27RA','IL6ST'],'IFNG_receptor':['IFNGR1','IFNGR2'],'IFNB_receptor':['IFNAR1','IFNAR2'],'IL12_receptor':['IL12RB1','IL12RB2'],'IL15_receptor':['IL2RB','IL2RG'],'IL18_receptor':['IL18R1','IL18RAP'],'IL27_ligand':['IL27','EBI3'],'IL12_ligand':['IL12A','IL12B'],'IL15_presentation':['IL15','IL15RA']}
(O/'provenance/specification.json').write_text(json.dumps({'genes':G,'complexes':pairs,'primary':'GSE174725 fixed reviewed states, donor-level expression and codetection; mouse original complete-family upstream statistics','human_AUC_excluded':'cached Human SCENIC file is CWP, not GSE174725','sensitivity':'leave-one-donor descriptive effects; frozen doublet/mixed flag exclusion; descriptive exact composition decomposition','GSE296682':'patient localization only, no group tests or 2-donor correlation','no_cell_threshold':True},indent=2),encoding='utf-8')
def wr(d,n):d.to_csv(O/'tables'/n,sep='\t',index=False,na_rep='NA')
def bh(v):
 p=np.array(v,float);q=p.copy();ok=np.isfinite(p);q[ok]=false_discovery_control(p[ok]) if ok.any() else [];return q
hf=H/'broad_annotation_review_SERVER_ONLY.h5ad';mf=H/'frozen_annotation_cells_SERVER_ONLY.tsv.gz';a=ad.read_h5ad(hf);m=pd.read_csv(mf,sep='\t',index_col=0).loc[a.obs_names];x=sparse.csr_matrix(a.layers['counts']);gi=pd.Index(a.var_names.str.upper());samples=sorted(m['sample'].unique());rows=[];co=[]
def flag(c):return m[c].astype(str).str.lower().isin(['true','1'])
exclude=flag('predicted_doublet')|flag('lineage_mixed_flag');scopes={'Core':m.old_CD4.ne('NotInCore'),'MemCD4':m.old_CD4.eq('MemCD4'),'NaiveCD4':m.old_CD4.eq('NaiveCD4'),'Treg':m.old_CD4.eq('Treg'),'BALF_AM':m.review_label.eq('AM_like'),'FCN1_Mono_comparator':m.review_label.eq('FCN1_Mono'),'FCGR3A_Mono_comparator':m.review_label.eq('FCGR3A_Mono')}
for variant in ['Reviewed','Exclude_frozen_flags']:
 for scope,base in scopes.items():
  for sample in samples:
   mask=(base&m['sample'].eq(sample)&(~exclude if variant!='Reviewed' else True)).to_numpy();n=int(mask.sum());group=m.loc[m['sample']==sample,'group'].iloc[0]
   for g in G:
    j=gi.get_indexer([g])[0];z=x[mask,j].toarray().ravel() if n and j>=0 else None;rows.append(dict(variant=variant,scope=scope,sample=sample,group=group,gene=g,cells=n,feature_present=j>=0,CP10k=float(np.mean(z*1e4/m.loc[mask,'total_counts'])) if z is not None else np.nan,positive_cells=int((z>0).sum()) if z is not None else np.nan,pct=float((z>0).mean()*100) if z is not None else np.nan))
   for name,gs in pairs.items():
    ix=gi.get_indexer(gs);ok=(ix>=0).all();pos=(x[mask][:,ix].toarray()>0).all(axis=1) if ok and n else None;co.append(dict(variant=variant,scope=scope,sample=sample,group=group,complex=name,components=';'.join(gs),cells=n,measured=ok,positive_cells=int(pos.sum()) if pos is not None else np.nan,pct=float(pos.mean()*100) if pos is not None else np.nan))
rows=pd.DataFrame(rows);wr(rows,'human_donor_expression.tsv');wr(pd.DataFrame(co),'human_complex_codetection.tsv');stats=[];loo=[]
for (variant,scope,g),b in rows.groupby(['variant','scope','gene']):
 c=b[b.group=='Exposed'].CP10k.dropna().to_numpy();s=b[b.group=='Silicosis'].CP10k.dropna().to_numpy();pv=np.nan
 if min(len(c),len(s))>=2 and np.ptp(np.r_[c,s])>0:pv=permutation_test((c,s),lambda c,s:abs(s.mean()-c.mean()),n_resamples=np.inf,alternative='greater').pvalue
 delta=s.mean()-c.mean() if len(c) and len(s) else np.nan;stats.append(dict(variant=variant,scope=scope,gene=g,n_control=len(c),n_case=len(s),difference=delta,P=pv,family=f'human31_{variant}_{scope}'))
 for drop in samples:
  bb=b[b['sample']!=drop];cc=bb[bb.group=='Exposed'].CP10k.dropna();ss=bb[bb.group=='Silicosis'].CP10k.dropna();loo.append(dict(variant=variant,scope=scope,gene=g,dropped_sample=drop,n_control=len(cc),n_case=len(ss),difference=ss.mean()-cc.mean(),inferential_P_not_used=True))
stats=pd.DataFrame(stats);stats['BH']=stats.groupby('family').P.transform(bh);wr(stats,'human_donor_group_effects.tsv');wr(pd.DataFrame(loo),'human_leave_one_donor_effects.tsv')
# Exact descriptive decomposition of equal-donor Core means; no invented empty-state expression.
dec=[]
for g in G:
 groups={}
 for group in ['Exposed','Silicosis']:
  b=rows[(rows.variant=='Reviewed')&(rows.group==group)&(rows.gene==g)];core=b[b.scope=='Core'].set_index('sample');vals=[]
  for state in ['MemCD4','NaiveCD4','Treg']:
   z=b[b.scope==state].set_index('sample').reindex(core.index);p=z.cells/core.cells;contrib=(p*z.CP10k).where(p>0,0);mp=p.mean();mu=contrib.mean()/mp if mp>0 else np.nan;vals.append((state,mp,mu))
  groups[group]=pd.DataFrame(vals,columns=['state','p','mu']).set_index('state')
 for state in ['MemCD4','NaiveCD4','Treg']:
  c=groups['Exposed'].loc[state];s=groups['Silicosis'].loc[state];dec.append(dict(gene=g,state=state,control_fraction=c.p,case_fraction=s.p,composition_component=(s.p-c.p)*(s.mu+c.mu)/2,within_state_component=(s.mu-c.mu)*(s.p+c.p)/2))
dec=pd.DataFrame(dec);wr(dec,'human_Core_descriptive_decomposition.tsv')
for g,b in dec.groupby('gene'):
 expected=stats[(stats.variant=='Reviewed')&(stats.scope=='Core')&(stats.gene==g)].difference.iloc[0]
 if np.isfinite(expected):assert np.isclose((b.composition_component+b.within_state_component).sum(),expected,atol=1e-9)
# Mouse source effects and strict biological receptor requirements, no pooled composite-receptor substitution.
send=pd.read_csv(V/'02_sender/all_ligand_sample_expression.tsv.gz',sep='\t');rec=pd.read_csv(V/'03_receiver/receptor_candidate_sample_expression.tsv.gz',sep='\t');up=pd.read_csv(V/'04_chains/all_sender_receiver_associations.tsv.gz',sep='\t');out=[]
for source in ['AMlike','AllMac']:
 de=pd.read_csv(V/'02_sender'/f'{source}_DE_all_genes.tsv.gz',sep='\t')
 for lig in ['TNF','IL15','IL18','C3','CXCL16','IL27','IFNG','IL12A','IFNAB']:
  z=de[de.canonical==lig];u=up[(up.sender==source)&(up.x==lig)&(up.y=='AUC::Stat1(+)')].copy()
  for _,r in u.iterrows():
   q=r.to_dict();q.update(ligand_logFC=z.logFC.iloc[0] if len(z) else np.nan,ligand_P=z.PValue.iloc[0] if len(z) else np.nan,ligand_BH=z.FDR.iloc[0] if len(z) else np.nan);out.append(q)
wr(pd.DataFrame(out),'mouse_upstream_integrated.tsv')
manifest=[{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in [hf,mf,V/'04_chains/all_sender_receiver_associations.tsv.gz',V/'02_sender/AMlike_DE_all_genes.tsv.gz',V/'02_sender/AllMac_DE_all_genes.tsv.gz']];(O/'provenance/input_hashes.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print('HUMAN_FOCAL');print(stats[(stats.scope.isin(['Core','MemCD4','NaiveCD4','Treg']))&(stats.gene.isin(['STAT1','GSDMD','CXCR6','IL27RA']))].to_string(index=False));print('MOUSE');print(pd.DataFrame(out).query("model=='Pooled' and scope=='ActCD4'")[['sender','x','rho','P','BH','ligand_logFC','ligand_P','ligand_BH']].to_string(index=False));print('DECOMP');print(dec[dec.gene.isin(['STAT1','GSDMD'])].to_string(index=False))