from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
from stat_core import fit,bh,subsets,COVCOL
O=Path(__file__).resolve().parents[1];I=O/'inputs';R=O/'results';R.mkdir(exist_ok=True)
audit=pd.read_csv(I/'tables/AMlike_ligand_input_audit.tsv',sep='\t')
lr=pd.read_csv(I/'step3_ligand_receptor_sample_expression.tsv',sep='\t');lr=lr[(lr.lineage=='Mac')&(lr.scope=='AMlike')&(lr['filter']=='All')&lr.canonical.isin(['IL15','TNF'])]
assert len(lr)==12
chk=audit.merge(lr,on=['sample','canonical'],suffixes=('_M6','_M3'),validate='one_to_one')
assert (chk.cells_M6==chk.cells_M3).all() and np.allclose(chk.mean_CP10k_M6,chk.mean_CP10k_M3,atol=1e-14)
md=pd.read_csv(I/'step3_scope_sample_counts.tsv',sep='\t');md=md[(md.lineage=='Mac')&(md.scope=='AMlike')&(md['filter']=='All')]
s=audit.rename(columns={'canonical':'ligand','cells':'sender_cells','mean_CP10k':'sender_mean_CP10k'}).merge(md[['sample','library','cells']],on='sample',validate='many_to_one')
assert s.sender_cells.eq(s.cells).all();s['sender_library']=s.library/s.cells
s=s[['sample','group','ligand','sender_cells','sender_library','sender_mean_CP10k','pct','raw_sum']]
s.to_csv(R/'frozen_AMlike_two_ligand_samples.tsv',sep='\t',index=False)
# Reproduce only the already existing focused Step6 sender associations.
act=pd.read_csv(I/'tables/TF_sample_activity.tsv',sep='\t');prior=pd.read_csv(I/'tables/TF_associations.tsv',sep='\t')
prior=prior[(prior.kind=='AMlike')&prior.tf.isin(['STAT1','IKZF2'])];legacy=[];legacy_members=[]
for i,row in prior.iterrows():
 z=act[(act.key==row.key)&(act.level==row.level)&(act['mode']=='Drop7')&(act.tf==row.tf)].merge(s[s.ligand==row.feature].drop(columns='group'),on='sample',validate='one_to_one').sort_values('sample')
 g=fit(z,'sender_mean_CP10k','activity',['G']);base=fit(z,'sender_mean_CP10k','activity',[])
 loo=[fit(z[z['sample']!=who],'sender_mean_CP10k','activity',['G'],False)['rho'] for who in z['sample']]
 assert len(z)==6 and abs(g['rho']-row.rho)<1e-8
 assert abs(min(loo)-row.LOO_min)<1e-8 and abs(max(loo)-row.LOO_max)<1e-8
 legacy.append(dict(**row.to_dict(),source_row=int(i+2),n=6,raw_rho=base['rho'],raw_exact_P=base['p'],reproduced_group_rho=g['rho'],group_exact_P=g['p'],group_t_approx=g['p_t_approx'],group_exact_method=g['p_method'],rho_difference=abs(g['rho']-row.rho),sample_ids=';'.join(z['sample'])))
 z=z.assign(key=row.key,level=row.level,tf=row.tf);legacy_members.append(z)
leg=pd.DataFrame(legacy);leg['new_FDR_raw']=bh(leg.raw_exact_P);leg['new_FDR_group']=bh(leg.group_exact_P);leg.to_csv(R/'Step6_existing_sender_reproduced.tsv',sep='\t',index=False);pd.concat(legacy_members,ignore_index=True).to_csv(R/'Step6_reproduction_sample_inputs.tsv',sep='\t',index=False)
# Pair the unchanged frozen sender values to the later DropOverlap AUC.
d=pd.read_csv(I/'frozen_receiver_with_candidates.tsv',sep='\t');d=d[d.species=='Mouse'].copy()
drop=[c for c in d if c.startswith('sender_') or c in ['ligand','pct','raw_sum']];d=d.drop(columns=drop)
z=d.merge(s,on=['sample','group'],validate='many_to_many')
assert len(z)==len(d)*2 and not z.duplicated(['filter','scope','seed','tf','ligand','sample']).any()
z.to_csv(R/'mouse_sender_receiver_paired_inputs.tsv',sep='\t',index=False)
rosters=[];models=[];loo_rows=[]
for key,q in z.groupby(['filter','scope','tf','seed','ligand']):
 arm,scope,tf,seed,ligand=key;seed=int(seed);q=q.sort_values('sample')
 assert len(q)==6 and q.group.eq('Saline').sum()==3 and q.group.eq('Silica').sum()==3
 for threshold in [1,5,10,20]:
  for row in q.itertuples():rosters.append(dict(filter=arm,scope=scope,tf=tf,seed=seed,ligand=ligand,threshold=threshold,sample=row.sample,receiver_cells=row.receiver_cells,sender_cells=row.sender_cells,included=row.receiver_cells>=threshold and row.sender_cells>=threshold,reason='' if row.receiver_cells>=threshold and row.sender_cells>=threshold else 'cell_threshold'))
 candidates=['G','R','D']+(['T'] if scope=='ActCD4' else [])
 for covs in subsets(candidates):
  qq=q.dropna(subset=['sender_mean_CP10k','mean_auc']+[COVCOL[c] for c in covs]);res=fit(qq,'sender_mean_CP10k','mean_auc',covs)
  row=dict(species='Mouse',filter=arm,scope=scope,tf=tf,seed=seed,ligand=ligand,sender_scope='AMlike',sender_filter='All',activity='frozen_SCENIC_DropOverlap',covariates='+'.join(covs) or 'None',sample_ids=';'.join(qq['sample']),n_control=int(qq.group.eq('Saline').sum()),n_case=int(qq.group.eq('Silica').sum()),excluded=';'.join(sorted(set(q['sample'])-set(qq['sample']))),**res)
  row['model_id']='MS'+str(len(models)+1).zfill(4);models.append(row)
  for who in qq['sample']:
   f=fit(qq[qq['sample']!=who],'sender_mean_CP10k','mean_auc',covs,False);loo_rows.append(dict(model_id=row['model_id'],omitted_sample=who,loo_rho=f['rho'],loo_df=f['df'],loo_status=f['status']))
m=pd.DataFrame(models);m['FDR']=np.nan;ok=m.p.notna();m.loc[ok,'FDR']=bh(m.loc[ok,'p']);m['FDR_family']='Mouse_sender_fixed_IL15_TNF_2TF_Core_ActCD4_All_BG1_3seeds_all_predeclared_models';m['nominal_positive']=m.p.lt(.05)
m.to_csv(R/'mouse_sender_all_models.tsv',sep='\t',index=False);pd.DataFrame(rosters).to_csv(R/'mouse_pairing_threshold_audit.tsv',sep='\t',index=False);pd.DataFrame(loo_rows).to_csv(R/'mouse_sender_LOO.tsv',sep='\t',index=False)
mm=m.groupby(['filter','scope','tf','ligand','covariates']).agg(n=('n','min'),rho_min=('rho','min'),rho_max=('rho','max'),P_min=('p','min'),P_max=('p','max'),FDR_min=('FDR','min'),FDR_max=('FDR','max'),nominal_seeds=('nominal_positive','sum'),df_min=('df','min'),statuses=('status',lambda x:';'.join(sorted(set(x))))).reset_index();mm.to_csv(R/'mouse_sender_seed_summary.tsv',sep='\t',index=False)
old=pd.read_csv(I/'prior_aligned_all_models.tsv.gz',sep='\t',low_memory=False);old.covariates=old.covariates.fillna('None');old=old[(old.species=='Mouse')&(old.branch=='Receiver')&(old.threshold==1)&(old.sample_policy=='available_cases')]
assert old.sample_ids.nunique()==1 and old.sample_ids.iloc[0]==m.sample_ids.iloc[0]
old.to_csv(R/'same_six_library_candidate_module_models.tsv',sep='\t',index=False)
side=[]
for key,qq in m.groupby(['filter','scope','tf','ligand','seed']):
 arm,scope,tf,ligand,seed=key
 for cov in ['None','G','G+R']:
  sender=qq[qq.covariates==cov].iloc[0]
  for target in old[(old['filter']==arm)&(old.scope==scope)&(old.tf==tf)&(old.seed==seed)&(old.covariates==cov)].itertuples():
   side.append(dict(filter=arm,scope=scope,tf=tf,seed=seed,ligand=ligand,covariates=cov,sample_ids=sender.sample_ids,sender_rho=sender.rho,sender_P=sender.p,sender_FDR=sender.FDR,receiver_evidence=target.evidence,receiver_feature=target.feature,receiver_rho=target.rho,receiver_P=target.p,receiver_FDR=target.FDR,interpretation='same six libraries; two associations are not mediation or causality'))
pd.DataFrame(side).to_csv(R/'same_library_two_segment_evidence.tsv',sep='\t',index=False)
summary=dict(status='PASS',Step6_sender_rows_reproduced=len(leg),Step6_max_rho_difference=float(leg.rho_difference.max()),Step6_IKZF2_activity_rows=int(act.tf.eq('IKZF2').sum()),frozen_ligands=['IL15','TNF'],sender_filter='All',receiver_filters=['All','BG1'],receiver_scopes=['Core','ActCD4'],n_libraries=6,n_control=3,n_case=3,all_thresholds_same_six=all(a['included'] for a in rosters),sender_cells_min=int(s.sender_cells.min()),sender_cells_max=int(s.sender_cells.max()),receiver_cells_min=int(z.receiver_cells.min()),models=len(m),model_status=m.status.value_counts().to_dict(),same_roster_side_by_side_rows=len(side),scenic_recomputed=False,regulons_unchanged=True,original_edgeR_unchanged=True)
(O/'provenance/analysis_verification.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2));print(mm[(mm['filter']=='All')&mm.covariates.isin(['None','G','G+R','G+R+D','G+R+D+T'])].to_string(index=False))
