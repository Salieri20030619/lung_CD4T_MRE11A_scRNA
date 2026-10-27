from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import pandas as pd,numpy as np,json,hashlib,shutil,fcntl
from scipy import stats
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA'); A=P/'02_analysis';O=P/'11_reports/human_first_20261004_v1';T=P/'11_reports/technical_delivery_20261002_ikzf2_stat1_v1/revision_v2/final_release'
lock=open(O/'provenance/stage1.lock','w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
assert not (O/'provenance/stage1_complete.json').exists()
manifest=[]
def rd(p):
 p=Path(p);manifest.append(dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()));return pd.read_csv(p,sep='\t')
def wr(d,n):pd.DataFrame(d).to_csv(O/'results'/n,sep='\t',index=False,na_rep='NA')
def bh(p):
 p=np.asarray(p);order=np.argsort(p);v=p[order]*len(p)/np.arange(1,len(p)+1);out=np.empty(len(p));out[order]=np.minimum(1,np.minimum.accumulate(v[::-1])[::-1]);return out
def fit(d,y,model):
 d=d[np.isfinite(d[y])&np.isfinite(d.library)].copy()
 if model=='WithinG227':d=d[d.series.eq('G227')]
 g=d.group.eq('CWP').astype(float).to_numpy();X=np.column_stack([np.ones(len(d)),g]);names=['Intercept','CWP']
 if model in ['Source','SourceDep']:
  z=pd.get_dummies(d.series,drop_first=True,dtype=float);X=np.column_stack([X,z]);names+=list(z.columns)
 if model=='SourceDep':X=np.column_stack([X,np.log1p(d.library)]);names+=['log1p_depth']
 n=len(d);rank=np.linalg.matrix_rank(X);df=n-rank
 out=dict(n=n,n_case=int(g.sum()),n_control=int(n-g.sum()),samples=';'.join(d['sample']),design=';'.join(names),df=df,effect=np.nan,P=np.nan,CI_low=np.nan,CI_high=np.nan,status='INSUFFICIENT')
 if min(out['n_case'],out['n_control'])<2:return out
 if rank<X.shape[1] or df<2:out['status']='NON_IDENTIFIABLE';return out
 yy=d[y].to_numpy();b=np.linalg.lstsq(X,yy,rcond=None)[0];res=yy-X@b;se=np.sqrt((res@res/df)*np.linalg.inv(X.T@X)[1,1]);p=2*stats.t.sf(abs(b[1]/se),df) if se>0 else np.nan
 out.update(effect=b[1],P=p,CI_low=b[1]-stats.t.ppf(.975,df)*se,CI_high=b[1]+stats.t.ppf(.975,df)*se,status='FIT' if se>0 else 'CONSTANT');return out
cells=rd(A/'mechanism_step2_20260927_154919_EsN4gk/input/human_cells.tsv');assert cells.cell.is_unique
design=cells[['sample','group','series']].drop_duplicates();assert design['sample'].is_unique and len(design)==19
design['donor']=design['sample'];design['clinical_label']=design.group.map({'CWP':'coal workers pneumoconiosis','Control':'control'});design['clinical_detail_status']='disease group available; severity/stage/exposure not verified'
wr(design,'stage1_donor_design.tsv');wr(cells.assign(donor=cells['sample'],source=cells.series),'stage1_cell_dictionary.tsv.gz')
scores=rd(T/'source_tables/M2_human_sample_scores.tsv');wr(scores,'stage1_reused_all_program_inputs.tsv')
legacy=rd(T/'data/M2_all_module_group_models.tsv');fixed=rd(A/'mechanism_step3_20260928_103057_lIRL38/tables/human_source_models_BH.tsv');fixed=fixed[fixed.metric.str.startswith('score_')]
assert len(legacy)==len(fixed)==504
keys=['scope','filter','model','metric'];ch=legacy.merge(fixed,on=keys,suffixes=('_old','_fixed'),validate='one_to_one');assert np.allclose(ch.p_nominal_old,ch.p_nominal_fixed,equal_nan=True)
wr(legacy,'stage1_legacy_504_original_FDR.tsv');wr(fixed,'stage1_legacy_504_repaired_BH.tsv')
mods=['score_TAct','score_TReg','score_Inhib','score_Naive'];rows=[]
for scope,z in scores[scores['filter'].eq('All')].groupby('scope'):
 for roster,zz in [('PositiveCells',z[z.cells>0]),('Cells20',z[z.cells>=20])]:
  for model in ['Base','Source','SourceDep','WithinG227']:
   for metric in mods:
    r=dict(kind='program',scope=scope,denominator='within_scope',roster=roster,model=model,metric=metric,filter='All')
    if roster=='Cells20' and model!='Base':
     q=fixed[(fixed.scope==scope)&(fixed['filter']=='All')&(fixed.model==model)&(fixed.metric==metric)].iloc[0]
     rr=fit(zz,metric,model);assert np.isclose(rr['P'],q.p_nominal,equal_nan=True,rtol=1e-7)
     r.update(rr,origin='legacy_statistic_verified',legacy_FDR=q.FDR)
    else:r.update(fit(zz,metric,model),origin='new_finite_model',legacy_FDR=np.nan)
    rows.append(r)
counts=cells.groupby(['sample','subtype']).size().unstack(fill_value=0).reindex(design['sample']).fillna(0)
fra=[]
for denom,ss in [('Primary',list(counts)),('Core',[c for c in counts if c!='CD4state'])]:
 for sub in ss:
  q=design.copy();q['scope']=sub;q['denominator']=denom;q['cells']=q['sample'].map(counts[sub]);q['denominator_cells']=q['sample'].map(counts[ss].sum(axis=1));q['fraction']=q.cells/q.denominator_cells
  q=q.merge(scores[(scores.scope==denom)&scores['filter'].eq('All')][['sample','library']],on='sample',validate='one_to_one');fra.append(q)
  for model in ['Base','Source','SourceDep','WithinG227']:rows.append(dict(kind='fraction',scope=sub,denominator=denom,roster='All19',model=model,metric='fraction',filter='All',origin='new_finite_model',**fit(q,'fraction',model)))
wr(pd.concat(fra),'stage1_proportions_inputs.tsv');r=pd.DataFrame(rows);r['BH_family']=r.kind+'|'+r.denominator+'|'+r.roster+'|'+r.model;r['FDR']=np.nan
for _,ix in r.groupby('BH_family').groups.items():
 ix=[i for i in ix if np.isfinite(r.loc[i,'P'])];r.loc[ix,'FDR']=bh(r.loc[ix,'P'])
wr(r,'stage1_all_finite_models.tsv')
# Counts, paired availability and annotations, preserving every reviewed source.
rev=rd(T/'prior_review_repaired/scenic3_review/results/reviewed_cells.tsv.gz');raw=rd(A/'scenic3_20261001_123452_Kvp9h5/input/human_mac_cells.tsv');assert set(rev.cell)==set(raw.cell)
merged=raw.merge(rev[['cell','review_scope','label','r4']],on='cell',validate='one_to_one');mapped=merged[merged['sample'].isin(design['sample'])]
wr(merged.assign(donor=merged['sample'],source=merged.series),'stage1_sender_cell_dictionary.tsv.gz')
for name in ['annotation_frozen.tsv','cluster_donor_dominance.tsv','sender_composition.tsv']:
 d=rd(P/'11_reports/front_half_20261004_v1/results'/name);wr(d,'stage1_'+name)
rc=cells.groupby(['sample','subtype']).size().unstack(fill_value=0);rc['Primary']=rc.sum(axis=1);rc['Core']=rc['Primary']-rc['CD4state']
sc=mapped.groupby(['sample','review_scope']).size().unstack(fill_value=0);sc['MacWide']=sc.MacAM+sc.MonoMac
pair=[]
for rec in rc:
 for sen in ['MacAM','MonoMac','MacWide']:
  q=design.copy();q['receiver']=rec;q['sender']=sen;q['receiver_cells']=q['sample'].map(rc[rec]).fillna(0);q['sender_cells']=q['sample'].map(sc[sen]).fillna(0);q['paired20']=(q.receiver_cells>=20)&(q.sender_cells>=20);pair.append(q)
wr(pd.concat(pair),'stage1_pairing_all_sources.tsv')
front=rd(P/'11_reports/front_half_20261004_v1/results/group_context_models.tsv');wr(front,'stage1_reused_front_half_125.tsv');q=front[front.kind.eq('receiver_state')];assert len(q)==30 and not (q.P<.05).any()
wr(manifest,'stage1_input_manifest.tsv')
summary={'donors':19,'cells':len(cells),'subtypes':cells.subtype.value_counts().to_dict(),'legacy_models':len(fixed),'finite_models':len(r),'new_models':int(r.origin.eq('new_finite_model').sum()),'nominal':int((r.P<.05).sum()),'BH05':int((r.FDR<.05).sum()),'clinical_limit':'severity and patient metadata beyond disease/source not verified','front_half_receiver_30_nominal':0}
(O/'provenance/stage1_complete.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2));print(r[(r.P<.05)&r.roster.eq('PositiveCells')][['scope','model','metric','n','effect','P','FDR']].to_string(index=False));print('FRACTIONS',r[(r.kind=='fraction')&(r.P<.05)][['scope','denominator','model','effect','P','FDR']].to_string(index=False))
