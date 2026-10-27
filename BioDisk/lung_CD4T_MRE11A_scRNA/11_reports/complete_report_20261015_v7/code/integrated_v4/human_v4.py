import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
from pathlib import Path
import sys,json,hashlib,itertools,re
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from scipy.stats import permutation_test,false_discovery_control
O=Path(sys.argv[1]);P=O.parents[1];H=P/'11_reports/GSE174725_annotation_validation_20261005_v1/results';CAND=['IKZF2','STAT1','IL1R2','PENK','TNFRSF9','LAG3','TTC39C','FOXP1','AHNAK','MDFIC','ST6GALNAC3'];MODS=['PyroExec','Inflam','ApoExec','OXPHOS','MitoISR','MitoROS','Mitophagy'];manifest=[]
def rd(f,**kw):
 f=Path(f);manifest.append(dict(path=str(f),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()));return pd.read_csv(f,sep='\t',**kw)
def wr(a,n):a.to_csv(O/'05_human'/n,sep='\t',index=False,encoding='utf-8',na_rep='NA')
def bh(p):
 a=np.array(p,float);q=np.full(a.shape,np.nan);ok=np.isfinite(a)
 if ok.any():q[ok]=false_discovery_control(a[ok])
 return q
f=H/'broad_annotation_review_SERVER_ONLY.h5ad';manifest.append(dict(path=str(f),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()));a=ad.read_h5ad(f);md=rd(H/'frozen_annotation_cells_SERVER_ONLY.tsv.gz',index_col=0).loc[a.obs_names];x=sparse.csr_matrix(a.layers['counts']);lib=md.total_counts.to_numpy();gi=pd.Index(a.var_names.str.upper());assert gi.is_unique and np.allclose(x.data,np.rint(x.data));samples=sorted(md['sample'].unique())
lr=rd(O/'inputs/LR_prior.tsv');query=sorted(set(lr.ligand)|set(lr.receptor)|set(CAND)|{'IL15RA','IL2RB','IL2RG','IL18RAP'});ix=gi.get_indexer(query);ok=ix>=0;xp=x[:,ix[ok]];y=xp.multiply((1e4/lib)[:,None]).tocsr();rows=[];roster=[]
scopes={'BALF_AM':md.review_label.eq('AM_like'),'FCN1_Mono_comparator':md.review_label.eq('FCN1_Mono'),'FCGR3A_Mono_comparator':md.review_label.eq('FCGR3A_Mono'),'Core':md.old_CD4.ne('NotInCore')}
for s in ['MemCD4','NaiveCD4','Treg']:scopes[s]=md.old_CD4.eq(s)
for scope,keep in scopes.items():
 for sample in samples:
  mask=(keep&md['sample'].eq(sample)).to_numpy();n=int(mask.sum());group=md.loc[md['sample']==sample,'group'].iloc[0];roster.append(dict(scope=scope,sample=sample,group=group,cells=n))
  mu=np.full(len(query),np.nan);pct=mu.copy();umi=mu.copy()
  if n:mu[ok]=np.asarray(y[mask].mean(axis=0)).ravel();pct[ok]=np.asarray((xp[mask]>0).mean(axis=0)).ravel()*100;umi[ok]=np.asarray(xp[mask].sum(axis=0)).ravel()
  for j,g in enumerate(query):rows.append(dict(scope=scope,sample=sample,group=group,cells=n,gene=g,feature_present=bool(ok[j]),CP10k=mu[j],pct=pct[j],UMI=umi[j]))
expr=pd.DataFrame(rows);wr(expr,'reviewed_expression_all_prior_features.tsv.gz');roster=pd.DataFrame(roster);wr(roster,'reviewed_sample_scope_roster.tsv')
# Phenotype panels re-evaluated in the reviewed CD4 roster, with controls fixed without labels.
modules=rd(P/'02_analysis/mechanism_step1_20260927_152218_s4V5eF/tables/module_definition.tsv');allmembers=set(modules.gene);cd4=scopes['Core'].to_numpy();cdmd=md.loc[cd4];raw=x[cd4];z=raw.multiply((1e4/lib[cd4])[:,None]).tocsr();z.data=np.log1p(z.data)
ref=np.mean([np.asarray(z[cdmd['sample'].eq(s).to_numpy()].mean(axis=0)).ravel() for s in samples],axis=0);pool=np.array([g not in allmembers and g not in CAND and not re.match(r'^(MT-|RPL|RPS|MRPL|MRPS|IGH|IGK|IGL)',g) for g in gi])&(ref>0)
bg='SFTPA1 SFTPA2 SFTPB SFTPC SFTPD SCGB1A1 SCGB3A1 SCGB3A2 ABCA3 SLC34A2 NAPSA LYZ LYZ1 LYZ2 S100A8 S100A9 CSF1R LST1 TYROBP FCER1G'.split();pool[np.isin(gi,bg)]=False;bins=pd.qcut(pd.Series(ref).rank(method='first'),25,labels=False).to_numpy();cellscore=pd.DataFrame(index=cdmd.index);coverage=[];controls=[]
for mod in MODS:
 gs=modules[modules.module==mod].gene.to_list();mi=gi.get_indexer(gs);mi=mi[mi>=0];ctrl=set()
 for j in mi:
  cand=np.flatnonzero(pool&(bins==bins[j]));order=np.lexsort((cand,abs(ref[cand]-ref[j])));ctrl.update(cand[order[:20]])
 ctrl=np.array(sorted(ctrl),int);expressed=int((np.asarray(raw[:,mi].sum(axis=0)).ravel()>0).sum());eligible=len(mi)>=3 and expressed>=3 and len(mi)/len(gs)>=.6 and len(ctrl)>=20
 cellscore[mod]=(np.asarray(z[:,mi].mean(axis=1)).ravel()-np.asarray(z[:,ctrl].mean(axis=1)).ravel()) if eligible else np.nan
 coverage.append(dict(module=mod,expected=len(gs),measured=len(mi),expressed=expressed,controls=len(ctrl),eligible=eligible));controls.extend(dict(module=mod,gene=gi[j],role='control') for j in ctrl)
wr(pd.DataFrame(coverage),'program_coverage.tsv');wr(pd.DataFrame(controls),'program_control_members.tsv')
score=[]
for scope in ['Core','MemCD4','NaiveCD4','Treg']:
 mask=scopes[scope].loc[cdmd.index]
 for sample in samples:
  take=mask&cdmd['sample'].eq(sample);n=int(take.sum());r=dict(scope=scope,sample=sample,group=cdmd.loc[cdmd['sample']==sample,'group'].iloc[0],cells=n)
  for mod in MODS:r[mod]=cellscore.loc[take,mod].mean() if n else np.nan
  score.append(r)
score=pd.DataFrame(score);wr(score,'reviewed_sample_program_scores.tsv')
def effect(c,s):
 c=np.asarray(c,float);s=np.asarray(s,float);p=np.nan
 if min(len(c),len(s))>=2 and np.ptp(np.r_[c,s])>0:p=permutation_test((c,s),lambda a,b:abs(b.mean()-a.mean()),vectorized=False,n_resamples=np.inf,alternative='greater').pvalue
 return dict(n_control=len(c),n_case=len(s),control_mean=c.mean() if len(c) else np.nan,case_mean=s.mean() if len(s) else np.nan,difference=s.mean()-c.mean() if len(c) and len(s) else np.nan,P=p)
r=[]
for scope,a0 in score.groupby('scope'):
 for mod in MODS:
  c=a0[a0.group=='Exposed'][mod].dropna();s=a0[a0.group=='Silicosis'][mod].dropna();b=effect(c,s);b.update(scope=scope,module=mod,family='Human_program_'+('Core7' if scope=='Core' else 'Subtypes21'));r.append(b)
r=pd.DataFrame(r);r['BH']=r.groupby('family').P.transform(bh);wr(r,'human_program_group_tests.tsv')
r=[]
for scope in ['Core','MemCD4','NaiveCD4','Treg']:
 for g in CAND:
  v=expr[(expr.scope==scope)&(expr.gene==g)];c=v[v.group=='Exposed'].CP10k.dropna();s=v[v.group=='Silicosis'].CP10k.dropna();b=effect(c,s);b.update(scope=scope,gene=g,family='Human_RNA_'+('Core11' if scope=='Core' else 'Subtypes33'));r.append(b)
r=pd.DataFrame(r);r['BH']=r.groupby('family').P.transform(bh);wr(r,'human_candidate_RNA_group_tests.tsv')
# Paired AM-to-receiver roster. Absent macrophages remain NA, not monocyte replacements.
pair=[]
for scope in ['Core','MemCD4','NaiveCD4','Treg']:
 r=roster[roster.scope==scope].set_index('sample');s=roster[roster.scope=='BALF_AM'].set_index('sample')
 for sample in samples:pair.append(dict(receiver_scope=scope,sample=sample,group=r.loc[sample,'group'],receiver_cells=int(r.loc[sample,'cells']),sender_cells=int(s.loc[sample,'cells']),paired_nonempty=bool(r.loc[sample,'cells']>0 and s.loc[sample,'cells']>0)))
wr(pd.DataFrame(pair),'human_AM_receiver_pairing.tsv')
# Preserve existing patient-only localization; use detection proportion as supplied.
f=P/'11_reports/GSE296682_patient_localization_20261005_v1/results/all_fixed_and_marker_localization.tsv.gz';loc=rd(f);target=loc[loc.gene.isin(CAND+['IL15','IL18','TNF'])].copy();valid=(target.N>0)&target.pct.notna()&target.n_positive.notna();assert np.allclose(target.loc[valid,'pct'],target.loc[valid,'n_positive']/target.loc[valid,'N'],atol=1e-12);wr(target,'patient_localization_existing.tsv')
wr(pd.DataFrame([dict(gene=g,extracted=g in set(loc.gene),interpretation='patient-only localization; no case-control inference' if g in set(loc.gene) else 'not extracted in legacy panel; unknown, not negative') for g in CAND+['IL15','IL18','TNF']]),'patient_localization_coverage.tsv')
wr(pd.DataFrame(manifest).drop_duplicates('path'),'input_manifest.tsv')
summary={'status':'PASS','reviewed_CD4_cells':int(cd4.sum()),'human_donors':5,'AM_case_donors_with_cells':int(((roster.scope=='BALF_AM')&(roster.group=='Silicosis')&(roster.cells>0)).sum()),'AM_group_or_pair_inferential_tests':0,'cell_number_gate':'none','new_CWP_models':0}
(O/'provenance/human_complete.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(json.dumps(summary));rr=pd.read_csv(O/'05_human/human_candidate_RNA_group_tests.tsv',sep='\t');print(rr[rr.gene.isin(['IKZF2','STAT1'])][['scope','gene','difference','P','BH','n_control','n_case']].to_string(index=False))