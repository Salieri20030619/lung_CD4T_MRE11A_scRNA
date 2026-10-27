from pathlib import Path
import os,json,hashlib,fcntl,datetime
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import pandas as pd,numpy as np
from scipy import stats
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/phenotype_candidate_discovery_20261005_v1';R=O/'results';H=P/'11_reports/human_first_20261004_v1/results'
lock=open(O/'provenance/compute.lock','w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
assert not (O/'provenance/compute_complete.json').exists()
(O/'provenance/DISCOVERY_FROZEN_SPEC.txt').write_bytes((O/'code/DISCOVERY_FROZEN_SPEC.txt').read_bytes())
manifest=[]
def rd(p):
 p=Path(p);manifest.append(dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()));return pd.read_csv(p,sep='\t',low_memory=False)
def wr(d,n):pd.DataFrame(d).to_csv(R/n,sep='\t',index=False,na_rep='NA')
def bh(d,keys,p='P',out='FDR'):
 d[out]=np.nan;d[out+'_family']=''
 for k,ii in d.groupby(keys,dropna=False).groups.items():
  ix=[i for i in ii if np.isfinite(d.loc[i,p])]
  if ix:d.loc[ix,out]=stats.false_discovery_control(d.loc[ix,p]);d.loc[ix,out+'_family']=str(k)+';m='+str(len(ix))
 return d
def ols(d,y,model):
 d=d[np.isfinite(d[y])&np.isfinite(d.library)].copy()
 if model=='WithinG227':d=d[d.series=='G227']
 g=d.group.eq('CWP').to_numpy(dtype=float);X=np.column_stack([np.ones(len(d)),g])
 if model in ['Source','SourceDep']:X=np.column_stack([X,pd.get_dummies(d.series,drop_first=True,dtype=float)])
 if model=='SourceDep':X=np.column_stack([X,np.log1p(d.library)])
 n=len(d);rank=np.linalg.matrix_rank(X);df=n-rank;r=dict(n=n,n_case=int(g.sum()),n_control=int(n-g.sum()),samples=';'.join(d['sample']),df=df,effect=np.nan,P=np.nan,CI_low=np.nan,CI_high=np.nan,status='INSUFFICIENT')
 if min(r['n_case'],r['n_control'])<2:return r
 if rank<X.shape[1] or df<2:r['status']='NON_IDENTIFIABLE';return r
 b=np.linalg.lstsq(X,d[y],rcond=None)[0];res=d[y]-X@b;se=np.sqrt((res@res/df)*np.linalg.inv(X.T@X)[1,1])
 if se<=1e-12:r['status']='CONSTANT';return r
 ci=stats.t.ppf(.975,df)*se;r.update(effect=b[1],P=2*stats.t.sf(abs(b[1]/se),df),CI_low=b[1]-ci,CI_high=b[1]+ci,status='FIT');return r
def corr(d,step,x='x',y='y'):
 n=len(d);r=dict(n=n,n_case=int(d.group.eq('CWP').sum()),n_control=int(d.group.eq('Control').sum()),samples=';'.join(d['sample']),rho=np.nan,P=np.nan,CI_low=np.nan,CI_high=np.nan,status='INSUFFICIENT',df=np.nan)
 if min(r['n_case'],r['n_control'])<2:return r
 cols=[np.ones(n)]
 if step!='Unadjusted':cols.append(d.group.eq('CWP').astype(float))
 if step in ['Source','SourceBothDepth']:cols.append(pd.get_dummies(d.series,drop_first=True,dtype=float))
 if step=='SourceBothDepth':cols.extend([stats.rankdata(d.sender_depth),stats.rankdata(d.receiver_depth)])
 X=np.column_stack(cols);rank=np.linalg.matrix_rank(X);df=n-rank-1;r['df']=df
 if rank<X.shape[1] or df<3:r['status']='NON_IDENTIFIABLE_OR_DF_LT3';return r
 a=stats.rankdata(d[x]);b=stats.rankdata(d[y]);a-=X@np.linalg.lstsq(X,a,rcond=None)[0];b-=X@np.linalg.lstsq(X,b,rcond=None)[0]
 if min(np.std(a),np.std(b))<1e-10:r['status']='CONSTANT_OR_NO_RESIDUAL_VARIATION';return r
 rho=float(np.clip(np.corrcoef(a,b)[0,1],-1,1));pv=2*stats.t.sf(abs(rho)*np.sqrt(df/max(1e-15,1-rho*rho)),df);ci=np.tanh(np.arctanh(np.clip(rho,-1+1e-12,1-1e-12))+np.array([-1,1])*stats.norm.ppf(.975)/np.sqrt(df-1));r.update(rho=rho,P=pv,CI_low=ci[0],CI_high=ci[1],status='SAME_DATA_DISCOVERY_APPROX_CI');return r
scores=rd(H/'stage1_reused_all_program_inputs.tsv');scores['receiver']=scores.scope.str.replace('Subtype:','',regex=False);mods=[c for c in scores if c.startswith('score_')];assert len(mods)==12
des=rd(H/'stage1_donor_design.tsv');meta=rd(H/'stage2_scope_sample_counts.tsv');members=rd(H/'stage1_complete_program_members.tsv');wr(members,'original_program_members.tsv');wr(des,'donor_design_SERVER_ONLY.tsv')
rows=[]
for (rec,filt),v in scores.groupby(['receiver','filter']):
 for roster,z in [('Cells20',v[v.cells>=20]),('PositiveCells',v[v.cells>0])]:
  for model in ['Base','Source','SourceDep','WithinG227']:
   for metric in mods:rows.append(dict(receiver=rec,filter=filt,roster=roster,model=model,metric=metric,**ols(z,metric,model)))
phen=bh(pd.DataFrame(rows),['filter','roster','model']);phen['FDR_global']=np.nan;ix=phen.index[phen.P.notna()];phen.loc[ix,'FDR_global']=stats.false_discovery_control(phen.loc[ix,'P']);wr(phen,'all_12_program_phenotype_models.tsv')
old=rd(H/'stage1_legacy_504_repaired_BH.tsv');old['receiver']=old.scope.str.replace('Subtype:','',regex=False);check=phen[phen.roster=='Cells20'].merge(old,on=['receiver','filter','model','metric'],suffixes=('_new','_old'));assert len(check)==504 and np.allclose(check.P,check.p_nominal,equal_nan=True,rtol=1e-6);wr(check,'legacy_504_reproduction.tsv.gz')
frac=rd(H/'stage1_all_finite_models.tsv');frac=frac[frac.kind=='fraction'];wr(frac,'all_composition_phenotypes.tsv')
# Freeze eligible phenotype endpoints before reading ligand association outcomes.
eligible=phen[(phen.P<.05)&(~phen.receiver.isin(['Core','Primary']))][['receiver','metric']].drop_duplicates();eligible['kind']='program';eligible['denominator']='within_scope'
ef=frac[frac.P<.05][['scope','metric','denominator']].drop_duplicates().rename(columns={'scope':'receiver'});ef['kind']='fraction';endpoints=pd.concat([eligible,ef],ignore_index=True).sort_values(['kind','receiver','metric','denominator']);wr(endpoints,'frozen_nominal_phenotype_endpoints.tsv');(O/'provenance/endpoint_freeze.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'endpoints':endpoints.to_dict('records'),'rule':'any predeclared phenotype model nominal; not ligand selected'},indent=2))
print('PHENOTYPES',len(phen),'nominal',int(phen.P.lt(.05).sum()),'endpoints',len(endpoints),flush=True);print(endpoints.to_string(index=False),flush=True)
expr=rd(H/'stage2_LR_expression_all_donors.tsv.gz');existing=rd(H/'stage2_evidence_all_candidates.tsv.gz');assoc=rd(H/'stage2_same_donor_associations.tsv.gz');assert len(existing)==231168 and len(assoc)==17596
old15=set('IKZF2 STAT1 IL15 TNF IL18 ICAM2 VEGFB PLAU RTN4 GSTP1 CXCL2 DSC2 CLSTN1 ADAM15 VTN'.split());existing['previously_studied_15']=existing.ligand.isin(old15);existing['source_identity']=existing.sender.map({'MacAM':'reviewed_macrophage','MonoMac':'mixed_monocyte_macrophage','MacWide':'union_not_independent'});existing['FDR_family']=existing[['receiver','filter','model','test','direction']].astype(str).agg('|'.join,axis=1)
family=assoc.groupby(['receiver','filter','model','test','direction']).agg(rows=('P','size'),finite=('P','count'),nominal=('P',lambda x:int((x<.05).sum())),BH05=('FDR',lambda x:int((x<.05).sum()))).reset_index();wr(family,'original_association_BH_families.tsv')
wr(existing,'complete_231168_candidate_audit.tsv.gz');wr(assoc,'original_17596_associations.tsv.gz')
targets=rd(H/'stage2_targetset_members.tsv.gz');bg=rd(H/'stage2_background_members.tsv.gz');ov=[]
for key,z in targets.groupby(['receiver','filter','model','test','direction']):
 for prog,gs in members.groupby('program'):
  overlap=sorted(set(z.gene)&set(gs.gene));ov.append(dict(zip(['receiver','filter','model','test','direction'],key),program=prog,program_n=len(gs),target_n=len(z),overlap_n=len(overlap),members=';'.join(overlap)))
wr(ov,'all_12_program_actual_disease_target_overlap.tsv')
fracinputs=rd(H/'stage1_proportions_inputs.tsv');rows=[];datasets=[]
alllig=sorted(existing.ligand.unique());wr(pd.DataFrame({'ligand':alllig}),'full_prior_688_ligands.tsv')
for e in endpoints.to_dict('records'):
 rec,metric,kind,den=e['receiver'],e['metric'],e['kind'],e['denominator']
 if kind=='program':base=scores[(scores.receiver==rec)&(scores['filter']=='All')][['sample','group','series','cells','library',metric]].rename(columns={'cells':'receiver_cells','library':'receiver_depth',metric:'y'})
 else:
  base=fracinputs[(fracinputs.scope==rec)&(fracinputs.denominator==den)][['sample','group','series','denominator_cells','library','fraction']].rename(columns={'denominator_cells':'receiver_cells','library':'receiver_depth','fraction':'y'})
 for sender in ['MacAM','MonoMac']:
  sm=meta[(meta.lineage=='Mac')&(meta.scope==sender)][['sample','cells','mean_depth']].rename(columns={'cells':'sender_cells','mean_depth':'sender_depth'});xx=expr[(expr.key==sender)&expr.gene.isin(alllig)][['sample','gene','mean_CP10k','pct']]
  d=base.merge(sm,on='sample').merge(xx,on='sample').rename(columns={'gene':'ligand','mean_CP10k':'x'});d=d[np.isfinite(d[['x','y','sender_depth','receiver_depth']]).all(axis=1)&d.sender_cells.ge(20)&d.receiver_cells.ge(20)]
  for lig in alllig:
   v=d[d.ligand==lig];tag=dict(e,sender=sender,ligand=lig,filter='All');
   for step in ['Unadjusted','Group','Source','SourceBothDepth']:rows.append(dict(tag,step=step,**corr(v,step)))
  datasets.append(d.assign(**e,sender=sender))
 models=bh(pd.DataFrame(rows),['kind','receiver','metric','denominator','step']);models['FDR_global']=np.nan;ii=models.index[models.P.notna()];models.loc[ii,'FDR_global']=stats.false_discovery_control(models.loc[ii,'P']);wr(models,'complete_direct_phenotype_discovery_models.tsv.gz')
 print('ENDPOINT_DONE',rec,metric,den,len(rows),flush=True)
wr(pd.concat(datasets,ignore_index=True),'direct_phenotype_donor_inputs_SERVER_ONLY.tsv.gz')
# All gene-endpoint decision records, not a P-sorted winner list.
keys=['kind','receiver','metric','denominator','sender','ligand'];trace=models.pivot(index=keys,columns='step',values=['rho','P','FDR','n','df']).reset_index();trace.columns=['_'.join(x).rstrip('_') for x in trace.columns];trace['previously_studied_15']=trace.ligand.isin(old15);trace['any_nominal']=trace[[c for c in trace if c.startswith('P_')]].lt(.05).any(axis=1);trace['source_nominal']=trace.P_Source.lt(.05);trace['source_depth_nominal']=trace.P_SourceBothDepth.lt(.05);trace['same_sign_4']=trace[[c for c in trace if c.startswith('rho_')]].apply(lambda r:len(set(np.sign(r.dropna())))==1 and r.notna().sum()==4,axis=1)
lr=rd(H/'stage2_all_LR_supported_and_unsupported.tsv.gz');lr=lr[(lr['filter']=='All')&(lr.gate=='Pct10')];lrmap=lr.groupby(['receiver','sender','from']).pair_supported.any().to_dict();trace['LR_10pct_compatible']=[lrmap.get((x.receiver,x.sender,x.ligand),False) for x in trace.itertuples()]
ee=existing[(existing['filter']=='All')&(existing.test=='Nominal05')];positive=ee.groupby(['receiver','sender','ligand']).aupr_corrected.max().to_dict();trace['any_nominal_target_positive_activity']=[positive.get((x.receiver,x.sender,x.ligand),np.nan)>0 for x in trace.itertuples()]
trace['decision_reason']=np.select([trace.previously_studied_15,trace.sender.ne('MacAM'),~trace.any_nominal,~trace.source_nominal & ~trace.source_depth_nominal],['previously_studied_kept_in_audit','mixed_source_not_macrophage_specific','no_nominal_direct_phenotype_link','nominal_only_before_source_adjustment'],default='new_macrophage_direct_phenotype_lead_review_all_axes')
wr(trace,'all_gene_endpoint_decision_audit.tsv.gz');leads=trace[trace.decision_reason=='new_macrophage_direct_phenotype_lead_review_all_axes'];wr(leads,'new_macrophage_leads_all_axes.tsv')
wr(manifest,'input_manifest.tsv');summary=dict(phenotype_models=len(phen),phenotype_finite=int(phen.P.notna().sum()),phenotype_nominal=int(phen.P.lt(.05).sum()),endpoints=len(endpoints),original_candidates=len(existing),original_associations=len(assoc),new_models=len(models),new_finite=int(models.P.notna().sum()),new_nominal=int(models.P.lt(.05).sum()),new_BH05=int(models.FDR.lt(.05).sum()),gene_endpoint_audit=len(trace),new_macrophage_lead_records=len(leads));(O/'provenance/compute_complete.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2),flush=True);print(leads[['receiver','metric','denominator','ligand','rho_Group','P_Group','rho_Source','P_Source','FDR_Source','rho_SourceBothDepth','P_SourceBothDepth','FDR_SourceBothDepth','n_Source','same_sign_4','LR_10pct_compatible','any_nominal_target_positive_activity']].to_string(index=False),flush=True)
