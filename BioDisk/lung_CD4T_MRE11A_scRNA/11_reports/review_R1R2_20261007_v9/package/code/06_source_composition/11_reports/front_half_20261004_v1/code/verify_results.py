from pathlib import Path
import os,json,hashlib,datetime
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np,pandas as pd
from scipy import stats
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/front_half_20261004_v1');R=O/'results'
def rd(n):return pd.read_csv(R/n,sep='\t',low_memory=False)
m=rd('new_link_models.tsv');inp=rd('new_link_model_inputs.tsv.gz');gm=rd('group_context_models.tsv');loo=rd('new_link_LOO.tsv.gz');checks=[]
assert len(m)==1675 and len(rd('finite_test_matrix.tsv'))==185 and m.model_id.is_unique
ix={k:v.sort_values('sample') for k,v in inp.groupby('model_id')}
for _,r in m.iterrows():
 z=ix.get(r.model_id,pd.DataFrame(columns=['sample']));assert len(z)==r.n
 assert ';'.join(z['sample'])==(r.sample_ids if isinstance(r.sample_ids,str) else '')
 if r.status!='FIT':continue
 cs=[] if pd.isna(r.covariates) or r.covariates=='None' else r.covariates.split('+');columns=[np.ones(len(z))]
 for c in cs:
  if c=='G':columns.append(z.group.eq('CWP').to_numpy(float))
  elif c=='S':columns.extend([z.series.eq(v).to_numpy(float) for v in sorted(z.series.unique())[1:]])
  else:columns.append(stats.rankdata(z[c]))
 D=np.column_stack(columns);xx=stats.rankdata(z.x_value);yy=stats.rankdata(z.y_value)
 a=xx-D@np.linalg.lstsq(D,xx,rcond=None)[0];b=yy-D@np.linalg.lstsq(D,yy,rcond=None)[0];rho=np.corrcoef(a,b)[0,1]
 X=np.column_stack([D,xx]);beta=np.linalg.lstsq(X,yy,rcond=None)[0];df=len(z)-np.linalg.matrix_rank(X);variance=np.sum((yy-X@beta)**2)/df*np.linalg.inv(X.T@X)[-1,-1]
 p=2*stats.t.sf(abs(beta[-1])/np.sqrt(variance),df)
 checks.append(dict(model_id=r.model_id,rho_difference=abs(rho-r.rho),P_difference=abs(p-r.P),origin=r.statistic_origin))
 assert abs(rho-r.rho)<1e-9 and abs(p-r.P)<1e-8,(r.model_id,rho,r.rho,p,r.P)
for _,q in m[m.sample_policy.eq('fixed_complete')].groupby(['kind','scope','sender','x','y']):assert q.sample_ids.fillna('').nunique()==1
for fam,q in m[m.P.notna()].groupby('FDR_family'):
 u=q.drop_duplicates('test_hash');new=stats.false_discovery_control(u.P.to_numpy(),method='bh');assert np.allclose(new,u.FDR,atol=1e-12)
assert np.allclose(stats.false_discovery_control(gm.loc[gm.P.notna(),'P'],method='bh'),gm.loc[gm.P.notna(),'FDR'],atol=1e-12)
ginput=rd('group_context_inputs.tsv');gchecks=[]
for _,v in gm.iterrows():
 z=ginput[ginput.model_id.eq(v.model_id)].sort_values('sample');assert len(z)==v.n
 if v.status!='FIT':continue
 cols=[np.ones(len(z))]
 if '+S' in v.covariates:cols.extend([z.series.eq(s).to_numpy(float) for s in sorted(z.series.unique())[1:]])
 if 'library' in v.covariates:cols.append(stats.rankdata(z.depth_value))
 C=np.column_stack(cols);g=z.group.eq('CWP').to_numpy(float);y=z.y_value.to_numpy();gr=g-C@np.linalg.lstsq(C,g,rcond=None)[0];yr=y-C@np.linalg.lstsq(C,y,rcond=None)[0]
 effect=gr@yr/(gr@gr);rho=np.corrcoef(gr,yr)[0,1];df=len(z)-np.linalg.matrix_rank(C)-1;p=2*stats.t.sf(abs(rho)*np.sqrt(df/(1-rho*rho)),df)
 assert abs(effect-v.effect)<1e-9 and abs(p-v.P)<1e-8
 gchecks.append(dict(model_id=v.model_id,effect_difference=abs(effect-v.effect),P_difference=abs(p-v.P)))
pd.DataFrame(gchecks).to_csv(R/'independent_group_checks.tsv',sep='\t',index=False)
for _,r in m[(m.status=='FIT')&m.statistic_origin.eq('new_finite_analysis')].iterrows():
 v=loo[(loo.model_id==r.model_id)&(loo.status=='FIT')].rho.dropna();assert len(v)==r.LOO_valid
 if len(v):assert np.allclose([v.min(),v.max()],[r.LOO_min,r.LOO_max])
pan=rd('panel_definition_overlap.tsv');assert not pan.loc[pan.used,'in_regulon_union'].any() and pan.loc[pan.used,'gene'].is_unique
manifest=rd('input_manifest.tsv');recheck=[]
for _,v in manifest.iterrows():
 p=Path(v.path);h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 recheck.append(dict(path=str(p),unchanged=h.hexdigest()==v.sha256))
assert all(v['unchanged'] for v in recheck)
pd.DataFrame(checks).to_csv(R/'independent_statistical_checks.tsv',sep='\t',index=False)
pd.DataFrame(recheck).to_csv(O/'provenance/input_recheck.tsv',sep='\t',index=False)
final=m[m.sample_policy.eq('fixed_complete')].copy();final['steps']=final.covariates.fillna('None').str.count('\\+')
final=final.sort_values('steps').groupby(['kind','scope','sender','x','y'],as_index=False).tail(1).drop(columns='steps');final.to_csv(R/'fixed_final_models.tsv',sep='\t',index=False)
baseline=m[m.sample_policy.eq('fixed_complete')&m.covariates.fillna('None').eq('None')]
transition=baseline.merge(final,on=['kind','scope','sender','x','y'],suffixes=('_baseline','_final'),validate='one_to_one');transition.to_csv(R/'fixed_adjustment_transitions.tsv',sep='\t',index=False)
out=dict(status='PASS',time=datetime.datetime.now().isoformat(),specifications=len(m),hypotheses=len(final),new_FIT=int(((m.status=='FIT')&m.statistic_origin.eq('new_finite_analysis')).sum()),reused_statistics=int(m.statistic_origin.eq('reused_missing_links').sum()),not_estimable=int(m.status.ne('FIT').sum()),independent_rechecks=len(checks),independent_group_rechecks=len(gchecks),max_rho_difference=max(v['rho_difference'] for v in checks),max_P_difference=max(v['P_difference'] for v in checks),BH_verified=True,fixed_rosters_verified=True,input_files_unchanged=len(recheck),all_new_panels_disjoint_and_regulon_free=True)
(O/'provenance/verification.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
print('\nFINAL BY KIND');print(final.groupby(['kind','sender']).agg(rows=('model_id','size'),nominal=('nominal_P_lt_05','sum'),FDR05=('FDR_lt_05','sum')).to_string())
cols=['model_id','kind','scope','sender','x','y','n','rho','P','FDR','LOO_min','LOO_max','LOO_same_direction']
print('\nFINAL NOMINAL');print(final[final.P.lt(.05)][cols].sort_values(['kind','P']).to_string(index=False))
print('\nBASELINE VS FINAL SOURCE');print(transition[transition.kind.isin(['composition','sender_state'])][['kind','sender','x','y','n_baseline','rho_baseline','P_baseline','rho_final','P_final','FDR_final']].to_string(index=False))
print('\nGROUP SOURCE+DEPTH');print(gm[gm.covariates.str.contains('S')&gm.P.lt(.05)][['kind','scope','y','covariates','n','effect','P','FDR','LOO_min','LOO_max']].to_string(index=False))
