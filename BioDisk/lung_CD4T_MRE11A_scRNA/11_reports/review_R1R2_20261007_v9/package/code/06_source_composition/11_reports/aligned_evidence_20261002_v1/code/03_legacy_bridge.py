from pathlib import Path
import pandas as pd,numpy as np,json,importlib.util,hashlib
from stat_core import fit,bh
O=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('oldstats',O/'provenance/original_S2_stats.py');oldstats=importlib.util.module_from_spec(spec);spec.loader.exec_module(oldstats)
old=pd.read_csv(O/'inputs/existing_S2_TF_to_six_candidates.tsv',sep='\t');act=pd.read_csv(O/'inputs/S2_TF_sample_activity.tsv',sep='\t');act=act[act.variant=='DropOverlap']
hs=pd.read_csv(O/'inputs/M2_human_sample_scores.tsv',sep='\t');ms=pd.read_csv(O/'inputs/M1_sample_scores.tsv',sep='\t');ms=ms[ms.dataset=='GSE314943']
rows=[];members=[];cache={};cellchecks=[]
for (sp,arm,scope),a in act.groupby(['species','filter','scope']):
 cells=pd.read_csv(O/'inputs'/f'{sp}_{arm}_cells.tsv',sep='\t');z=cells if scope=='Core' else cells[cells.subtype==scope]
 counts=z.groupby('sample').size();ids=sorted(counts[counts>=20].index)
 dd=z[['sample','group','series']].drop_duplicates().set_index('sample').loc[ids].reset_index();dd['cells']=dd['sample'].map(counts);dd['library']=dd['sample'].map(z.groupby('sample').library.mean());dd['case']=dd.group.isin(['Silica','CWP']).astype(int)
 prop=pd.crosstab(z['sample'],z.subtype,normalize='index').reindex(ids)
 for j,c in enumerate(prop.columns):dd[f'prop_{j}']=prop[c].to_numpy(float)
 q=ms if sp=='Mouse' else hs;qscope=('All' if scope=='Core' else scope) if sp=='Mouse' else ('Core' if scope=='Core' else 'Subtype:'+scope)
 q=q[(q['filter']==arm)&(q.scope==qscope)].set_index('sample').reindex(ids)
 assert np.array_equal(q.cells.to_numpy(),dd.cells.to_numpy()) and np.allclose(q.library,dd.library)
 for c in q:
  if c.startswith('expr_'):dd[c]=q[c].to_numpy()
 dd['receiver_library']=dd.library
 cache[(sp,arm,scope)]=dd
 cellchecks.append(dict(species=sp,filter=arm,scope=scope,n=len(dd),cell_counts_equal=True,depth_equal=True,composition_columns=';'.join(prop.columns),sample_ids=';'.join(ids)))
for ix,row in old.iterrows():
 dd=cache[(row.species,row['filter'],row.scope)].copy();aa=act[(act.species==row.species)&(act['filter']==row['filter'])&(act.scope==row.scope)&(act.seed==row.seed)&(act.regulon==row.regulon)]
 dd=dd.merge(aa[['sample','mean_auc']],on='sample',how='left',validate='one_to_one');assert dd.mean_auc.notna().all()
 dd=dd[np.isfinite(dd[row.feature])&np.isfinite(dd.mean_auc)].copy()
 res=oldstats.association(dd[row.feature].to_numpy(),dd.mean_auc.to_numpy(),dd,row.model)
 covs=[] if row.model=='Spearman' else ['G']
 if row.model in ['SourceDep','SourceComp']:covs+=['S','R']
 if row.model=='Comp':covs+=['R']
 if row.model in ['Comp','SourceComp']:covs += [c for c in dd if c.startswith('prop_')][:-1]
 new=fit(dd,row.feature,'mean_auc',covs)
 diff=abs(res['rho']-row.rho) if np.isfinite(res['rho']) and np.isfinite(row.rho) else np.nan
 pdiff=abs(res['p']-row.p) if np.isfinite(res['p']) and np.isfinite(row.p) else np.nan
 assert np.isfinite(res['rho'])==np.isfinite(row.rho) and (not np.isfinite(diff) or diff<1e-8)
 assert np.isfinite(res['p'])==np.isfinite(row.p) and (not np.isfinite(pdiff) or pdiff<1e-8)
 assert np.isfinite(new['rho'])==np.isfinite(row.rho)
 assert not np.isfinite(row.rho) or abs(new['rho']-row.rho)<1e-8
 out=row.to_dict();out.update(legacy_reproduced_rho=res['rho'],legacy_reproduced_p=res['p'],rho_delta=diff,p_delta=pdiff,new_rho=new['rho'],new_p=new['p'],new_status=new['status'],new_df=new['df'],new_p_method=new['p_method'],new_near_perfect_low_df=new.get('near_perfect_low_df',False),sample_ids=';'.join(dd['sample']),new_covariates='+'.join(covs) or 'None',new_p_fills_old_missing=bool(pd.isna(row.p) and np.isfinite(new['p'])))
 rows.append(out)
 mm=dd.copy();mm.insert(0,'legacy_original_row',row.original_source_row);mm['species']=row.species;mm['filter']=row['filter'];mm['scope']=row.scope;mm['seed']=row.seed;mm['regulon']=row.regulon;mm['feature']=row.feature;mm['model']=row.model;members.append(mm)
r=pd.DataFrame(rows);r['new_FDR']=np.nan;r['new_FDR_family']='Legacy_bridge_v1/'+r.species
for _,idx in r.groupby('new_FDR_family').groups.items():
 ok=r.loc[idx].new_p.dropna().index;r.loc[ok,'new_FDR']=bh(r.loc[ok].new_p)
r.to_csv(O/'results/S2_all_1440_reproduced_and_calibrated.tsv',sep='\t',index=False)
pd.concat(members,ignore_index=True,sort=False).to_csv(O/'results/S2_model_inputs.tsv.gz',sep='\t',index=False)
pd.DataFrame(cellchecks).to_csv(O/'provenance/S2_scope_cell_alignment.tsv',sep='\t',index=False)
summary=dict(status='PASS',rows=len(r),max_rho_difference=float(r.rho_delta.max()),max_p_difference=float(r.p_delta.max()),new_P_filled=int(r.new_p_fills_old_missing.sum()),originals_unchanged=True,new_FDR='BH within species across fixed two TF x six genes x original scopes/filters/seeds/models; separate from old FDR',near_perfect_low_df=int(r.new_near_perfect_low_df.sum()),new_status_counts=r.new_status.value_counts().to_dict())
(O/'provenance/S2_reproduction.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
