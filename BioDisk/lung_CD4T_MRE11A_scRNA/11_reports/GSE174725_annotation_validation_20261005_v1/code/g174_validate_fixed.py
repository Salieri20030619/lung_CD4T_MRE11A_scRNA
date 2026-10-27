from pathlib import Path
import os,json,itertools,math
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np,pandas as pd,anndata as ad
from scipy.stats import t
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/GSE174725_annotation_validation_20261005_v1';R=O/'results'
def rd(n):return pd.read_csv(R/n,sep='\t')
def wr(d,n):pd.DataFrame(d).to_csv(R/n,sep='\t',index=False,na_rep='NA')
assert (O/'provenance/FROZEN_VALIDATION_SPEC.txt').exists()
a=ad.read_h5ad(R/'broad_annotation_review_SERVER_ONLY.h5ad');obs=pd.read_csv(R/'frozen_annotation_cells_SERVER_ONLY.tsv.gz',sep='\t',index_col=0).loc[a.obs_names];X=a.X;C=a.layers['counts'];samples=sorted(obs['sample'].unique())
members=rd('program_members_and_common_coverage.tsv');usable=members[(members.donors_covered==5)&~members.symbols.str.contains(';',na=False)]
def gs(s):
 q=usable[usable['set']==s];assert len(q)>0,s
 assert set(q.symbols)<=set(a.var_names),s
 return list(q.symbols)
def mean_genes(genes):return np.asarray(X[:,a.var_names.get_indexer(genes)].mean(axis=1)).ravel()
programs={}
for label,rec in [('Core','Core'),('CCR7','CD4-CCR7'),('GZMA','CD4-GZMA'),('IL7R','CD4-IL7R'),('Treg','Treg')]:programs[label+'_NomUp']=[('Disease_'+rec+'_Nominal05',1)]
programs['GZMA_FDRUp']=[('Disease_CD4-GZMA_FDR05',1)]
for label in ['TAct','Inhib','Naive']:programs[label+'_adjusted']=[('Module_'+label,1),('Module_'+label+'_background',-1)]
programs['TReg_minus_IKZF2']=[('Prior_TReg_minus_IKZF2_target',1),('Prior_TReg_minus_IKZF2_frozen_control',-1)]
programs['IL7R_excess']=[('Disease_CD4-IL7R_Nominal05',1),('MatchedControl1',-1/3),('MatchedControl2',-1/3),('MatchedControl3',-1/3)]
print('AVAILABLE_SETS',sorted(members['set'].unique()),flush=True)
scores={};rows=[]
for label,terms in programs.items():
 scores[label]=sum(w*mean_genes(gs(s)) for s,w in terms)
 for s,w in terms:
  for _,v in members[members['set']==s].iterrows():rows.append(dict(program=label,component=s,weight=w,gene=v.gene,symbol=v.symbols,donors_covered=v.donors_covered,included=v.donors_covered==5 and v.symbols in a.var_names))
wr(rows,'fixed_validation_full_program_members.tsv.gz')
scopes={'BALF_AM':obs.review_label.eq('AM_like'),'FCN1_Mono_comparator':obs.review_label.eq('FCN1_Mono'),'FCGR3A_Mono_comparator':obs.review_label.eq('FCGR3A_Mono'),'Core':obs.old_CD4.ne('NotInCore')}
for s in ['MemCD4','NaiveCD4','Treg']:scopes[s]=obs.old_CD4.eq(s)
L='IL15 TNF CLSTN1 ADAM15 PLAU DSC2 CXCL2 ICAM2 VTN VEGFB GSTP1 RTN4'.split();recgenes='IKZF2 STAT1 IL15RA IL2RB IL2RG TNFRSF1A TNFRSF1B'.split();lr=rd('fixed_candidate_receptor_pairs.tsv');allgenes=sorted(set(L+recgenes+list(lr.to)))
genevalues={g:np.asarray(X[:,a.var_names.get_loc(g)].toarray()).ravel() for g in allgenes if g in a.var_names};det={g:np.asarray(C[:,a.var_names.get_loc(g)].toarray()).ravel()>0 for g in genevalues}
donors=[];coverage=[]
for clean in [False,True]:
 for scope,base in scopes.items():
  for gsm in samples:
   mask=np.asarray(base&obs['sample'].eq(gsm)&(~obs.sensitivity_flag if clean else True));n=int(mask.sum());group='Exposed' if gsm in samples[:2] else 'Silicosis';meta=dict(scope=scope,sample=gsm,group=group,cell_filter='NoFlags' if clean else 'All',cells=n)
   genes=L+['IL15RA'] if scope not in ['Core','MemCD4','NaiveCD4','Treg'] else recgenes
   for g in genes:donors.append(dict(**meta,family='sender_RNA' if scope=='BALF_AM' else 'comparator_RNA' if 'comparator' in scope else 'receiver_RNA',endpoint=g,value=float(genevalues[g][mask].mean()) if n else np.nan))
   if scope in ['Core','MemCD4','NaiveCD4','Treg']:
    for name,score in scores.items():donors.append(dict(**meta,family='receiver_program',endpoint=name,value=float(score[mask].mean()) if n else np.nan))
   for g in allgenes:
    known=g in genevalues;coverage.append(dict(**meta,gene=g,mapped_common=known,positive_cells=int(det[g][mask].sum()) if known and n else np.nan,pct=float(det[g][mask].mean()) if known and n else np.nan,mean_log=float(genevalues[g][mask].mean()) if known and n else np.nan,status='NO_CELLS' if not n else 'UNMAPPED_COMMON' if not known else 'MEASURED'))
d=pd.DataFrame(donors);wr(d,'fixed_validation_donor_scores_SERVER_ONLY.tsv.gz');cv=pd.DataFrame(coverage);wr(cv,'fixed_LR_donor_coverage_SERVER_ONLY.tsv.gz')
tests=[]
for (cellfilter,scope,family,end),q in d.groupby(['cell_filter','scope','family','endpoint']):
 for minimum in [1,20]:
  u=q[(q.cells>=minimum)&q.value.notna()];x=u[u.group=='Exposed'].value.to_numpy();y=u[u.group=='Silicosis'].value.to_numpy();effect=y.mean()-x.mean() if len(x) and len(y) else np.nan;p=lo=hi=df=np.nan;status='INSUFFICIENT_DONORS';nperm=0
  if family=='comparator_RNA':status='DESCRIPTIVE_COMPARATOR'
  elif min(len(x),len(y))>=2:
   z=np.r_[x,y];vals=[]
   for choose in itertools.combinations(range(len(z)),len(y)):
    ix=np.zeros(len(z),bool);ix[list(choose)]=True;vals.append(z[ix].mean()-z[~ix].mean())
   p=float(np.mean(np.abs(vals)>=abs(effect)-1e-12));nperm=len(vals);v1=x.var(ddof=1)/len(x);v2=y.var(ddof=1)/len(y);se=np.sqrt(v1+v2)
   if se>0:df=(v1+v2)**2/(v1*v1/(len(x)-1)+v2*v2/(len(y)-1));margin=t.ppf(.975,df)*se;lo=effect-margin;hi=effect+margin
   status='COMPUTED_EXACT'
  tests.append(dict(cell_filter=cellfilter,min_cells=minimum,scope=scope,family=family,endpoint=end,n_exposed=len(x),n_silicosis=len(y),effect=effect,CI95_low=lo,CI95_high=hi,Welch_df=df,P_exact=p,permutations=nperm,status=status,donors=';'.join(u['sample'])))
out=pd.DataFrame(tests);out['FDR_BH']=np.nan;out['BH_family_size']=0
for _,ix in out.groupby(['cell_filter','min_cells','family']).groups.items():
 valid=out.loc[ix].dropna(subset=['P_exact']).sort_values('P_exact');n=len(valid)
 if n:out.loc[valid.index,'FDR_BH']=np.minimum.accumulate((valid.P_exact.to_numpy()*n/np.arange(1,n+1))[::-1])[::-1].clip(0,1);out.loc[ix,'BH_family_size']=n
wr(out,'all_fixed_validation_tests.tsv')
gates=[]
for cf in ['All','NoFlags']:
 for scope in ['Core','MemCD4','NaiveCD4','Treg']:
  for gsm in samples:
   sender=cv[(cv.cell_filter==cf)&(cv.scope=='BALF_AM')&(cv['sample']==gsm)].set_index('gene');recv=cv[(cv.cell_filter==cf)&(cv.scope==scope)&(cv['sample']==gsm)].set_index('gene');sn=int(sender.cells.iloc[0]);rn=int(recv.cells.iloc[0])
   vals=[sender.loc[g,'pct'] for g in ['IL15','IL15RA']]+[recv.loc[g,'pct'] for g in ['IL2RB','IL2RG']]
   for threshold in [0,.05,.10]:gates.append(dict(cell_filter=cf,scope=scope,sample=gsm,sender_cells=sn,receiver_cells=rn,threshold=threshold,IL15_pct=vals[0],IL15RA_pct=vals[1],IL2RB_pct=vals[2],IL2RG_pct=vals[3],all_marginal_detected=bool(all(np.isfinite(v) and (v>0 if threshold==0 else v>=threshold) for v in vals)),both_ends20=sn>=20 and rn>=20,interpretation='RNA_COMPATIBILITY_ONLY_NOT_TRANSPRESENTATION'))
wr(gates,'IL15_complex_coverage_gates.tsv')
# Secondary annotation evidence, no new clustering or post-hoc endpoint changes.
markers='C1QA C1QB C1QC MERTK CD163 CSF1R MRC1 MSR1 APOE APOC1 TREM2 GPNMB SPP1 LGALS3 FOLR2 LYVE1 TIMD4 CD68 CTSB CTSL LGMN FCN1 S100A8 VCAN FCGR3A FABP4 PPARG MARCO CD1C FCER1A CLEC10A CSF3R FCGR3B CD3D TRAC SCGB1A1 EPCAM KRT8 KRT18'.split();profiles=[]
for cl in sorted(obs.broad06.unique()):
 mask=np.asarray(obs.broad06.eq(cl))
 for g in markers:
  known=g in a.var_names;x=C[mask,a.var_names.get_loc(g)].toarray().ravel() if known else []
  profiles.append(dict(cluster=cl,review_label=obs.loc[mask,'review_label'].iloc[0],gene=g,mapped=known,cells=int(mask.sum()),pct=float((x>0).mean()) if known else np.nan,mean_count=float(x.mean()) if known else np.nan))
wr(profiles,'mature_nonAM_macrophage_marker_audit.tsv')
summary=dict(status='COMPLETE',cells=a.n_obs,common_features=a.n_vars,programs=len(programs),tests=len(out),computed=int(out.P_exact.notna().sum()),nominal=int((out.P_exact<.05).sum()),BH=int((out.FDR_BH<.05).sum()),spec=str(O/'provenance/FROZEN_VALIDATION_SPEC.txt'))
(O/'provenance/fixed_validation_complete.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True);print(out[(out.cell_filter=='All')&(out.min_cells==1)&(out.scope=='Core')].to_string(index=False),flush=True)
