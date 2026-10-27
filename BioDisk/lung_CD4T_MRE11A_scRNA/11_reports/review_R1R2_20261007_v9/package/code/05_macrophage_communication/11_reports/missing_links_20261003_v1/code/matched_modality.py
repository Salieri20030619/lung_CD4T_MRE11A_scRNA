"""Same frozen-AUC donor sensitivity for otherwise wider Human TF RNA models.
Declared during coverage review, before inspecting any association output.
Reuse the exact finite_models statistical functions, without running its matrix.
"""
from pathlib import Path
import ast
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/missing_links_20261003_v1')
src=(O/'code/finite_models.py').read_text()
# Compile only imports and statistical helper functions, never the model runner.
tree=ast.parse(src)
safe=[n for n in tree.body if isinstance(n,(ast.Import,ast.ImportFrom)) or isinstance(n,ast.FunctionDef) and n.name in ['bh','design','fit','covcols']]
exec(compile(ast.Module(body=safe,type_ignores=[]),'finite_models_statistical_helpers','exec'))
perms={}
lock=open(O/'analysis.lock','a+');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
R=O/'results';rec=pd.read_csv(R/'receiver_paired_modalities.tsv',sep='\t');pair=pd.read_csv(R/'sender_receiver_paired_inputs.tsv',sep='\t');sen=pd.read_csv(R/'sender_all_samples.tsv',sep='\t')
pair=pair.merge(sen[sen.sender.eq('TNK_other')][['species','sample','L_IFNG']].rename(columns={'L_IFNG':'TNK_IFNG'}),on=['species','sample'],how='left',validate='many_to_one')
m=pd.read_csv(R/'all_new_models.tsv',sep='\t',keep_default_na=False,na_values=['NA',''])
rows=[];members=[]
for _,r in m[m.species.eq('Human')&m.threshold.eq(1)&m.sample_policy.eq('available')].iterrows():
 tf=next((t for t in ['STAT1','IKZF2'] if r.x==t+'_RNA' or r.y==t+'_RNA'),None)
 if tf is None:continue
 source=pair[pair.sender.eq(r.sender)] if r.kind=='sender' else rec
 z=source[source.species.eq('Human')&source.scope.eq(r.scope)&source['sample'].isin(str(r.sample_ids).split(';'))].sort_values('sample')
 z=z[z[tf+'_AUC_17'].notna()]
 if len(z)==r.n:continue # same-sample result already exists; never recompute it
 covs=[] if r.covariates=='None' else r.covariates.split('+')
 result=fit(z,r.x,r.y,covs)
 q=r.to_dict();q.update(result);q.update(model_id='MATCH'+str(len(rows)+1).zfill(4),sample_policy='fixed_frozen_AUC_donors',FDR_family='sensitivity',sample_ids=';'.join(z['sample']),excluded_ids=';'.join(sorted(set(str(r.sample_ids).split(';'))-set(z['sample']))),missing_n=int(r.eligible_n)-len(z),parent_model_id=r.model_id)
 loo=[fit(z[z['sample']!=sample],r.x,r.y,covs,False)['rho'] for sample in z['sample']] if result['status']=='FIT' else [];loo=[v for v in loo if np.isfinite(v)]
 q.update(LOO_valid=len(loo),LOO_rho_min=min(loo) if loo else np.nan,LOO_rho_max=max(loo) if loo else np.nan,LOO_same_direction=bool(loo and all(np.sign(v)==np.sign(result['rho']) for v in loo)))
 rows.append(q)
 for _,v in z.iterrows():
  member=dict(model_id=q['model_id'],species='Human',sample=v['sample'],group=v.group,series=v.series,x_value=v[r.x],y_value=v[r.y],receiver_cells=v.receiver_cells)
  for c in covcols(covs):member[c]=v[c]
  members.append(member)
new=pd.DataFrame(rows)
allm=pd.concat([m,new],ignore_index=True);allm['FDR']=np.nan
for family,inds in allm.groupby('FDR_family').groups.items():
 u=allm.loc[inds];u=u[u.P.notna()].drop_duplicates('test_hash');mp=dict(zip(u.test_hash,bh(u.P)))
 allm.loc[inds,'FDR']=allm.loc[inds,'test_hash'].map(mp)
allm['nominal_positive']=allm.P.lt(.05);allm['FDR_positive']=allm.FDR.lt(.05)
allm.to_csv(R/'all_new_models_with_matched.tsv',sep='\t',index=False,na_rep='NA')
allm[allm.model_id.str.startswith('MATCH')].to_csv(R/'matched_modality_models.tsv',sep='\t',index=False,na_rep='NA')
pd.DataFrame(members).to_csv(R/'matched_modality_sample_inputs.tsv.gz',sep='\t',index=False,na_rep='NA',compression='gzip')
(O/'provenance/matched_modality_complete.json').write_text(json.dumps(dict(rows=len(new),reason='two Human raw-RNA donors lack frozen AUC; fix RNA roster to AUC donors without imputing AUC',original_models_recomputed=False),indent=2))
print((O/'provenance/matched_modality_complete.json').read_text())
