from pathlib import Path
import numpy as np,pandas as pd,json,hashlib
from scipy import stats
import itertools,shutil
ROOT=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/stepwise_diagnostic_20261002_v1')
I=ROOT/'inputs';O=ROOT/'results';O.mkdir(exist_ok=True)
d=pd.read_csv(I/'receiver_all_positive_cell_samples.tsv',sep='\t')
s=pd.read_csv(I/'sender_IL18_MacAM_all_samples.tsv',sep='\t');s=s[s.sender_filter=='All']
MODULES=['OXPHOS','MitoISR','MitoROS','Mitophagy','PyroExec','ApoExec','TAct','TReg','Inhib','Inflam','Naive','Cycle']
results=[];members=[]
PERMUTATIONS={}

def bh(p):
 p=np.asarray(p,float);order=np.argsort(p);q=np.empty(len(p));q[order]=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];return np.minimum(q,1)

def fit_group(y,case):
 y=np.asarray(y,float);X=np.column_stack([np.ones(len(y)),np.asarray(case,float)]);rank=np.linalg.matrix_rank(X);df=len(y)-rank
 r=dict(estimate=np.nan,p=np.nan,df=df,design_rank=rank,status='FIT',effect_type='case_minus_control_OLS',p_method='OLS_t_exploratory')
 if rank<X.shape[1] or df<=0:r['status']='UNIDENTIFIABLE';return r
 if np.ptp(y)<1e-14:r['status']='CONSTANT_OUTCOME';return r
 b=np.linalg.lstsq(X,y,rcond=None)[0];res=y-X@b;se=float(np.sqrt(np.sum(res*res)/df*np.linalg.inv(X.T@X)[1,1]))
 r.update(estimate=float(b[1]),se=se)
 if se<=1e-14:r['status']='DEGENERATE_RESIDUAL_VARIANCE';return r
 t=b[1]/se;r.update(p=float(2*stats.t.sf(abs(t),df)),CI_low=float(b[1]-stats.t.ppf(.975,df)*se),CI_high=float(b[1]+stats.t.ppf(.975,df)*se))
 return r

def fit_spearman(x,y):
 x=np.asarray(x,float);y=np.asarray(y,float);n=len(x)
 r=dict(estimate=np.nan,p=np.nan,df=n-2,design_rank=1,status='FIT',effect_type='Spearman_rho',p_method='rank_correlation_t_approximation_exploratory')
 if n<3:r['status']='UNIDENTIFIABLE';return r
 if np.ptp(x)<1e-14 or np.ptp(y)<1e-14:r['status']='CONSTANT_VARIABLE';return r
 rho=float(np.corrcoef(stats.rankdata(x),stats.rankdata(y))[0,1]);rho=float(np.clip(rho,-1,1))
 approx=float(2*stats.t.sf(abs(rho)*np.sqrt((n-2)/max(1-rho*rho,1e-15)),n-2))
 r.update(estimate=rho,p=approx,p_t_approx=approx)
 if n<=9:
  if n not in PERMUTATIONS:PERMUTATIONS[n]=np.asarray(list(itertools.permutations(range(n))),dtype=np.int16)
  rx=stats.rankdata(x);ry=stats.rankdata(y);rx-=rx.mean();ry-=ry.mean()
  permrho=ry[PERMUTATIONS[n]]@rx/(np.linalg.norm(rx)*np.linalg.norm(ry))
  r['p']=float(np.mean(np.abs(permrho)>=abs(rho)-1e-12));r['p_method']='exact_two_sided_rank_permutation'
 return r

def add(z,evidence,feature,tf,seed,fn,xcol=None,ycol=None,ytransform='identity'):
 z=z.copy()
 need=[ycol]+([xcol] if xcol else [])
 z=z.replace([np.inf,-np.inf],np.nan).dropna(subset=need)
 z=z[z.receiver_cells>=1]
 if 'sender_cells' in z:z=z[z.sender_cells>=1]
 z=z.sort_values('sample')
 assert z['sample'].is_unique
 y=z[ycol].to_numpy(float)
 if ytransform=='log1p_CP10k':y=np.log1p(y)
 case=z.group.isin(['CWP','Silica']).astype(int)
 res=fn(y,case) if evidence.endswith('_GROUP') else fn(z[xcol],y)
 first=z.iloc[0] if len(z) else None
 row=dict(evidence=evidence,species=first.species if first is not None else '',filter=first['filter'] if first is not None else '',scope=first.scope if first is not None else '',tf=tf,feature=feature,seed=str(seed),stage='B0_unadjusted_min1',model='Unadjusted',threshold_receiver=1,threshold_sender=1 if evidence=='SENDER_TF' else np.nan,n=len(z),n_control=int((case==0).sum()),n_case=int((case==1).sum()),sample_ids=';'.join(z['sample']),group_distribution=json.dumps(z.group.value_counts().to_dict(),sort_keys=True),source_distribution=json.dumps(z.series.value_counts().to_dict(),sort_keys=True),receiver_cells_total=int(z.receiver_cells.sum()),receiver_cells_min=int(z.receiver_cells.min()) if len(z) else np.nan,receiver_cells_max=int(z.receiver_cells.max()) if len(z) else np.nan,sender_cells_total=float(z.sender_cells.sum()) if 'sender_cells' in z else np.nan,tf_detection_median_pct=float(z.tf_pct.median()) if tf else np.nan,tf_detection_min_pct=float(z.tf_pct.min()) if tf else np.nan,tf_detection_max_pct=float(z.tf_pct.max()) if tf else np.nan,sender_detection_median_pct=float(z.sender_pct.median()) if 'sender_pct' in z else np.nan,x_variable=xcol or 'case',y_variable=ycol,y_transform=ytransform,**res)
 row['small_sample_or_low_df']=row['n']<10 or row['df']<=3
 row['nominal_positive']=bool(np.isfinite(row['p']) and row['p']<.05)
 row['model_id']=f'B{len(results)+1:04d}'
 z['model_id']=row['model_id'];z['model_y']=y;z['model_x']=z[xcol].to_numpy() if xcol else case.to_numpy()
 results.append(row);members.append(z)

for (sp,arm,scope,tf),z in d.groupby(['species','filter','scope','tf'],sort=True):
 rr=z[z.seed==17].copy()
 add(rr,'RNA_GROUP',tf,tf,'NA',fit_group,ycol='tf_mean_CP10k',ytransform='log1p_CP10k')
for (sp,arm,scope,seed,tf),z in d.groupby(['species','filter','scope','seed','tf'],sort=True):
 add(z,'AUC_GROUP',tf,tf,seed,fit_group,ycol='mean_auc')
 for mod in MODULES:add(z,'TF_MODULE',mod,tf,seed,fit_spearman,xcol='mean_auc',ycol='score_'+mod)
 if sp=='Human':
  zz=z.merge(s,on='sample',suffixes=('','_sender'),validate='one_to_one')
  assert (zz.group==zz.group_sender).all() and (zz.series==zz.series_sender).all()
  add(zz,'SENDER_TF','IL18_MacAM',tf,seed,fit_spearman,xcol='sender_mean_CP10k',ycol='mean_auc')
for (sp,arm,scope),z in d[(d.seed==17)&(d.tf=='IKZF2')].groupby(['species','filter','scope'],sort=True):
 for mod in MODULES:add(z,'MODULE_GROUP',mod,'','NA',fit_group,ycol='score_'+mod)
r=pd.DataFrame(results);r['FDR_exploratory']=np.nan
r['FDR_family']='Baseline/'+r.species+'/'+r.evidence
for _,idx in r.groupby('FDR_family').groups.items():
 valid=r.loc[idx].p.dropna().index
 r.loc[valid,'FDR_exploratory']=bh(r.loc[valid].p)
r['FDR_status']=np.where(r.p.notna(),'NEW_BASELINE_EXPLORATORY_BH','NOT_ESTIMABLE')
if (O/'baseline_models.tsv').exists() and not (O/'baseline_tapprox_firstpass.tsv').exists():shutil.copy2(O/'baseline_models.tsv',O/'baseline_tapprox_firstpass.tsv')
r.to_csv(O/'baseline_models.tsv',sep='\t',index=False)
pd.concat(members,ignore_index=True).to_csv(O/'baseline_model_samples.tsv.gz',sep='\t',index=False)
r.groupby(['species','scope','filter','evidence']).agg(models=('model_id','size'),valid_p=('p','count'),nominal=('nominal_positive','sum'),n_min=('n','min'),n_max=('n','max')).reset_index().to_csv(O/'baseline_inventory.tsv',sep='\t',index=False)
notes={'stage':'BASELINE_PASS','models':len(r),'new_modeling_authorized':True,'biological_unit':'individual subject or independent mouse library, never cells or seeds','RNA_endpoint':'log1p of frozen original-library CP10k sample mean; exploratory OLS, not a replacement for original edgeR count model','AUC_endpoint':'frozen DropOverlap sample mean','association_method':'unadjusted Spearman, t approximation; small-sample P retained as explicitly unstable/exploratory','FDR_family':'new baseline exploratory BH within species and evidence class across all predeclared scopes, All/BG1, both TFs, all three seeds and all 12 modules; original FDR unchanged','available_scope':'widest frozen SCENIC scoring-cell populations, not reinstatement of raw upstream excluded cells','fixed_sender':'Human MacAM IL18 only; no ligand search','scope_separation':'Human Core and CD4-IL7R kept separate; Mouse Core and ActCD4 kept separate; no cross-species chain','no_causal_claim':True}
notes['association_method']='Unadjusted Spearman: exact two-sided enumeration for n<=9 (including ties), t approximation otherwise. Both P_t_approx and the calibrated primary P retained; no n<10 withholding. Small-sample instability still flagged.'
(O/'baseline_methods.json').write_text(json.dumps(notes,indent=2))
print('BASELINE_PASS',len(r),flush=True)
sel=r[(r['filter']=='All')&(((r.species=='Human')&(r.scope=='CD4-IL7R'))|((r.species=='Mouse')&(r.scope=='ActCD4')))]
print('BASELINE_GROUP_AND_SENDER')
print(sel[sel.evidence.isin(['RNA_GROUP','AUC_GROUP','SENDER_TF'])][['species','scope','evidence','tf','seed','n','n_control','n_case','estimate','p','FDR_exploratory','df','status']].to_string(index=False))
print('BASELINE_MODULE_GROUP')
print(sel[sel.evidence=='MODULE_GROUP'][['species','scope','feature','n','estimate','p','FDR_exploratory']].to_string(index=False))
print('BASELINE_TF_MODULE_SEED_RANGES')
print(sel[sel.evidence=='TF_MODULE'].groupby(['species','scope','tf','feature']).agg(rho_min=('estimate','min'),rho_max=('estimate','max'),p_min=('p','min'),p_max=('p','max'),nominal_seeds=('nominal_positive','sum')).to_string())
