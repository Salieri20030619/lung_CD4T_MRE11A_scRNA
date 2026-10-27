from pathlib import Path
import json, shutil, hashlib
import pandas as pd
import numpy as np
from scipy.stats import rankdata
from scipy.stats import spearmanr
D=Path('stepwise_diagnostic'); R=D/'results'; Q=D/'provenance'; Q.mkdir(exist_ok=True)
OLD=D/'reference_tables'; checks=[]; diffs=[]
def read(p):
 d=pd.read_csv(p,sep='\t');
 for c in ['tf','covariates']:
  if c in d:d[c]=d[c].fillna('None' if c=='covariates' else '')
 return d
def compare(label,a,b,keys,fields):
 m=a.merge(b,on=keys,suffixes=('_old','_new'),validate='one_to_one');assert len(m)==len(a),(label,len(a),len(m))
 summary={'check':label,'rows':len(m)}
 for x in fields:
  av=pd.to_numeric(m[x+'_old'],errors='coerce');bv=pd.to_numeric(m[x+'_new'],errors='coerce');valid=av.notna();assert bv[valid].notna().all(),label
  delta=(av-bv).abs();maxd=float(delta.max()) if valid.any() else 0
  summary[x+'_compared']=int(valid.sum());summary[x+'_max_abs_delta']=maxd
  m[x+'_abs_delta']=delta;assert maxd<1e-9,(label,x,maxd)
 summary['status']='PASS';checks.append(summary);m.to_csv(Q/(label+'_comparison.tsv'),sep='\t',index=False)
b=read(R/'baseline_models.tsv');r=read(R/'receiver_stepwise_models.tsv');s=read(R/'sender_stepwise_models.tsv')
keys=['species','scope','filter','tf','seed','evidence','feature']
rr=r[(r.threshold==1)&(r.covariates=='None')&(r.sample_policy=='fixed_complete_cases')]
compare('baseline_receiver_736',b[b.evidence!='SENDER_TF'],rr,keys,['estimate','p','n','df'])
ss=s[(s.threshold==1)&(s.covariates=='None')&(s.sample_policy=='fixed_complete_cases')&(s.variant=='Matched')&(s.analysis=='SENDER_TF_STEPS')].rename(columns={'rho':'estimate'})
compare('baseline_sender_24',b[b.evidence=='SENDER_TF'],ss,['species','scope','filter','tf','seed'],['estimate','p','n','df'])
old=read(OLD/'reviewed_sender_TF_models.tsv')
old=old[(old.sender=='MacAM')&(old.ligand=='IL18')&old.scope.isin(['Core','CD4-IL7R'])&old.tf.isin(['IKZF2(+)','STAT1(+)'])&~((old.sender_filter=='NoParent16')&(old.sample_set=='NoDom'))].copy()
old['tf']=old.tf.str.replace('(+)', '',regex=False);old['variant']=np.where(old.sender_filter=='NoParent16','NoParent16',old.sample_set)
old['covariates']=old.model.map({'Spearman':'None','SourceDep':'G+S+R','BothEpi':'G+S+R+D+E'})
ss=s[(s.threshold==20)&(s.sample_policy=='fixed_complete_cases')&(s.analysis=='SENDER_TF_STEPS')]
for lab,oo in [('original_sender_adjusted_144',old[old.model!='Spearman']),('original_sender_unadjusted_72',old[old.model=='Spearman'])]:
 compare(lab,oo,ss,['scope','filter','seed','tf','variant','covariates'],['rho','p','n','df'])
rr=r[(r.threshold==20)&(r.sample_policy=='fixed_complete_cases')]
old=read(OLD/'S3_mouse_TF_identity_adjustment.tsv');old['tf']=old.tf.str.upper();old['covariates']=old.model.map({'Group':'None','Depth':'R','TregAdj':'T','DepthTreg':'R+T'});old=old.rename(columns={'coef':'estimate'})
compare('mouse_identity_48',old,rr[(rr.species=='Mouse')&(rr.evidence=='AUC_GROUP')],['scope','filter','seed','tf','covariates'],['estimate','p','n','df'])
old=read(OLD/'S2_TF_group_effects.tsv');old=old[(old.variant=='DropOverlap')&old.scope.isin(['Core','CD4-IL7R','ActCD4'])&old.tf.isin(['IKZF2','STAT1'])].copy();old['covariates']=old.model.map({'Group':'None','Depth':'R','Source':'S','SourceDep':'S+R'});old=old[old.covariates.notna()].rename(columns={'coef':'estimate'});old['n']=old.n_control+old.n_case
compare('S2_AUC_groups',old,rr[rr.evidence=='AUC_GROUP'],['species','scope','filter','seed','tf','covariates'],['estimate','p','n','df'])
old=read(OLD/'S2_TF_associations.tsv');old=old[(old.variant=='DropOverlap')&old.scope.isin(['Core','CD4-IL7R','ActCD4'])&old.tf.isin(['IKZF2','STAT1'])&old.feature.str.startswith('score_')&old.model.isin(['Spearman','Group','SourceDep'])].copy();old=old[~((old.species=='Mouse')&(old.model=='SourceDep'))];old['covariates']=old.model.map({'Spearman':'None','Group':'G','SourceDep':'G+S+R'});old['feature']=old.feature.str.removeprefix('score_');old=old.rename(columns={'rho':'estimate'})
compare('S2_TF_module_associations',old,rr[rr.evidence=='TF_MODULE'],['species','scope','filter','seed','tf','covariates','feature'],['estimate','p','n','df'])
old=read(OLD/'M2_all_module_group_models.tsv');old=old[old.scope.isin(['Core','Subtype:CD4-IL7R'])&old.model.isin(['Source','SourceDep'])&old.metric.str.startswith('score_')].copy();old['scope']=old.scope.str.removeprefix('Subtype:');old['feature']=old.metric.str.removeprefix('score_');old['covariates']=old.model.map({'Source':'S','SourceDep':'S+R'});old=old.rename(columns={'coefficient':'estimate','p_nominal':'p','residual_df':'df'});old['n']=old.n_control+old.n_case
compare('M2_module_groups_96',old,rr[(rr.species=='Human')&(rr.evidence=='MODULE_GROUP')],['scope','filter','covariates','feature'],['estimate','p','n','df'])
def bh(p):
 p=np.array(p,float);order=np.argsort(p);q=np.empty(len(p));q[order]=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];return np.minimum(q,1)
for name,d in [('baseline',b),('receiver',r),('sender',s)]:
 assert d.model_id.is_unique
 assert (d.n==d.n_control+d.n_case).all()
 assert all(len(t.split(';'))==n and len(set(t.split(';')))==n for t,n in zip(d.sample_ids,d.n))
 assert d.p.dropna().between(0,1).all()
 assert (d.loc[d.p.notna(),'df']>0).all()
 assoc=d.evidence.isin(['TF_MODULE','SENDER_TF']) if 'evidence' in d else pd.Series(True,index=d.index)
 assert (d.df==d.n-d.design_rank-assoc.astype(int)).all()
 if 'threshold' in d:assert (d.receiver_cells_min>=d.threshold).all()
 if name=='sender':assert (d.sender_cells_min>=d.threshold).all()
 if 'complete_case_loss' in d:assert (d.complete_case_loss==0).all()
 u=d.drop_duplicates(['FDR_family','unique_test_id']) if 'unique_test_id' in d else d
 for _,f in u.groupby('FDR_family'):
  f=f[f.p.notna()];assert np.allclose(bh(f.p),f.FDR_exploratory,atol=1e-12)
 checks.append({'check':name+'_structural_audit','rows':len(d),'unique_tests':len(u),'status':'PASS','max_complete_case_loss':int(d.complete_case_loss.max()) if 'complete_case_loss' in d else None})
 # Membership tables independently verify every saved model roster, not merely counts in the model table.
 file={'baseline':'baseline_model_samples.tsv.gz','receiver':'receiver_stepwise_model_samples.tsv.gz','sender':'sender_stepwise_model_samples.tsv.gz'}[name]
 m=read(R/file);agg=m.groupby('model_id')['sample'].agg(lambda x:';'.join(sorted(x)))
 assert d.set_index('model_id').sample_ids.sort_index().equals(agg.sort_index()),name
 checks.append({'check':name+'_full_membership_roster','models':len(agg),'member_rows':len(m),'status':'PASS'})
# Explicit source structure and rank covariance decomposition at both fixed rosters.
inp=read(D/'inputs/receiver_all_positive_cell_samples.tsv');send=read(D/'inputs/sender_IL18_MacAM_all_samples.tsv');send=send[send.sender_filter=='All']
trans=read(R/'filter_transition_subset_removals.tsv');max_rho=max_p=0.
for row in trans.itertuples():
 z=inp[(inp.species=='Human')&(inp.scope==row.scope)&(inp['filter']==row.filter)&(inp.tf==row.tf)&(inp.seed==row.seed)].merge(send[['sample','sender_cells','sender_mean_CP10k']],on='sample',validate='one_to_one')
 before=z[(z.receiver_cells>=row.threshold_before)&(z.sender_cells>=row.threshold_before)]
 after=z[(z.receiver_cells>=row.threshold_after)&(z.sender_cells>=row.threshold_after)]
 removed=set() if row.removed_subset=='NONE' else set(row.removed_subset.split(';'))
 assert removed.issubset(set(before['sample'])-set(after['sample']))
 ids=set(row.sample_ids.split(';'));assert ids==set(before['sample'])-removed
 z=z[z['sample'].isin(ids)];rho,p=spearmanr(z.sender_mean_CP10k,z.mean_auc)
 max_rho=max(max_rho,abs(rho-row.rho));max_p=max(max_p,abs(p-row.p));assert len(z)==row.n and row.df==row.n-2
assert max_rho<1e-12 and max_p<1e-12
checks.append({'check':'all_filter_transition_subsets','rows':len(trans),'rho_max_abs_delta':max_rho,'p_max_abs_delta':max_p,'status':'PASS'})
out=[];composition=[]
for threshold in [1,20]:
 for seed in [17,29,43]:
  z=inp[(inp.species=='Human')&(inp.scope=='CD4-IL7R')&(inp['filter']=='All')&(inp.tf=='IKZF2')&(inp.seed==seed)].merge(send[['sample','sender_cells','sender_mean_CP10k','sender_pct']],on='sample',validate='one_to_one')
  z=z[(z.receiver_cells>=threshold)&(z.sender_cells>=threshold)].copy();z['xrank']=rankdata(z.sender_mean_CP10k);z['yrank']=rankdata(z.mean_auc)
  x=z.xrank-z.xrank.mean();y=z.yrank-z.yrank.mean();wx=z.xrank-z.groupby('series').xrank.transform('mean');wy=z.yrank-z.groupby('series').yrank.transform('mean')
  total=float(x@y);within=float(wx@wy);between=total-within
  out.append(dict(threshold=threshold,seed=seed,n=len(z),total_rank_crossproduct=total,within_source_rank_crossproduct=within,between_source_rank_crossproduct=between,unadjusted_rho=np.corrcoef(x,y)[0,1],source_residual_rho=np.corrcoef(wx,wy)[0,1],sender_rank_variance_removed=1-float(wx@wx/(x@x)),auc_rank_variance_removed=1-float(wy@wy/(y@y))))
  if seed==17:
   comp=z.groupby(['series','group']).agg(n=('sample','size'),sample_ids=('sample',lambda x:';'.join(x)),IL18_meanCP10k=('sender_mean_CP10k','mean'),IKZF2_AUC_mean=('mean_auc','mean')).reset_index();comp.insert(0,'threshold',threshold);composition.append(comp)
pd.DataFrame(out).to_csv(R/'source_rank_covariance_decomposition.tsv',sep='\t',index=False);pd.concat(composition).to_csv(R/'source_group_design.tsv',sep='\t',index=False)
old=read(OLD/'existing_RNA_DE_focused.tsv');old.to_csv(R/'original_edgeR_focused_unchanged.tsv',sep='\t',index=False)
pd.concat([b[b.evidence=='RNA_GROUP'].assign(method_origin='new_baseline_meanCP10k'),r[(r.evidence=='RNA_GROUP')&(r.threshold==1)&(r.sample_policy=='fixed_complete_cases')].assign(method_origin='new_adjusted_meanCP10k')],ignore_index=True).to_csv(R/'new_meanCP10k_RNA_comparison.tsv',sep='\t',index=False)
(Q/'full_reproduction_audit.json').write_text(json.dumps({'status':'PASS','checks':checks,'thresholds':[1,5,10,20],'original_FDR_unchanged':True,'old_small_n_NA_not_compared_to_new_authorized_P':True},indent=2))
pd.DataFrame(checks).to_csv(Q/'full_reproduction_audit_summary.tsv',sep='\t',index=False)
print(json.dumps(checks,indent=2));print(pd.DataFrame(out).to_string(index=False));print(pd.concat(composition).to_string(index=False))
