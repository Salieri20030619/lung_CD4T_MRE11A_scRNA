from pathlib import Path
import os,json,ast
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import pandas as pd,numpy as np
from scipy import stats
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/human_gap_diagnostics_20261005_v1';R=O/'results';V=P/'11_reports/human_alignment_20261005_v1/results';H=P/'11_reports/human_first_20261004_v1/results';A=P/'02_analysis'
L='IL15 TNF CLSTN1 ADAM15 PLAU DSC2 CXCL2 ICAM2 VTN VEGFB GSTP1 RTN4'.split();SC=['Core','CD4-CCR7','CD4-GZMA','CD4-IL7R','Treg'];steps=['Unadjusted','Group','Source','SourceSenderDepth','SourceBothDepth']
def rd(p):return pd.read_csv(p,sep='\t',low_memory=False)
def wr(d,n):pd.DataFrame(d).to_csv(R/n,sep='\t',index=False,na_rep='NA')
def bh(d,keys):
 d['FDR']=np.nan;d['BH_family']=''
 for key,ix in d.groupby(keys,dropna=False).groups.items():
  ix=[i for i in ix if np.isfinite(d.loc[i,'P'])]
  if ix:d.loc[ix,'FDR']=stats.false_discovery_control(d.loc[ix,'P']);d.loc[ix,'BH_family']=str(key)+';m='+str(len(ix))
 return d
def cor(d,step,x='x',y='y',extra=None):
 n=len(d);out=dict(n=n,n_case=int(d.group.eq('CWP').sum()),n_control=int(d.group.eq('Control').sum()),samples=';'.join(d['sample']),df=np.nan,rho=np.nan,CI_low=np.nan,CI_high=np.nan,P=np.nan,status='INSUFFICIENT')
 if min(out['n_case'],out['n_control'])<2:return out
 cols=[np.ones(n)]
 if step!='Unadjusted':cols.append(d.group.eq('CWP'))
 if step not in ['Unadjusted','Group']:cols.append(pd.get_dummies(d.series,drop_first=True,dtype=float))
 if step in ['SourceSenderDepth','SourceBothDepth']:cols.append(stats.rankdata(d.sender_depth))
 if step in ['SourceReceiverDepth','SourceBothDepth']:cols.append(stats.rankdata(d.receiver_depth))
 if extra:cols.append(stats.rankdata(d[extra]))
 X=np.column_stack(cols);rank=np.linalg.matrix_rank(X);df=n-rank-1;out['df']=df
 if rank<X.shape[1] or df<3:out['status']='NON_IDENTIFIABLE_OR_DF_LT3';return out
 a=stats.rankdata(d[x]);b=stats.rankdata(d[y]);a-=X@np.linalg.lstsq(X,a,rcond=None)[0];b-=X@np.linalg.lstsq(X,b,rcond=None)[0]
 if min(np.std(a),np.std(b))<1e-10:out['status']='NO_RESIDUAL_VARIATION';return out
 rho=float(np.clip(np.corrcoef(a,b)[0,1],-1,1));p=2*stats.t.sf(abs(rho)*np.sqrt(df/max(1e-15,1-rho*rho)),df);z=np.arctanh(np.clip(rho,-1+1e-12,1-1e-12));ci=np.tanh(z+np.array([-1,1])*1.95996398454/np.sqrt(df-1));out.update(rho=rho,P=p,CI_low=ci[0],CI_high=ci[1],status='EXPLORATORY_APPROX_FISHER_CI');return out
def complete(d,two=False,extra=None):
 cols=['x','y','receiver_depth']+(['sender_depth'] if two else [])+([extra] if extra else []);ok=np.isfinite(d[cols]).all(axis=1)&d.receiver_cells.ge(20)
 if two:ok&=d.sender_cells.ge(20)
 return d[ok].copy()
assert (O/'provenance/raw_complete.json').exists();assert 'GAP_EXPORT_COMPLETE' in (O/'provenance/gap_export.log').read_text();(O/'provenance/GAP_SPEC.txt').write_bytes((O/'code/GAP_SPEC.txt').read_bytes())
des=rd(V/'stage1_donor_design.tsv');gene=rd(V/'alignment_gene_donors.tsv.gz');tf=rd(V/'TF_readout_donors.tsv.gz');tf=tf[tf['filter']=='All'];prog=rd(V/'TF_excluded_program_donors.tsv.gz');meas=rd(R/'raw_donor_gene_measurements.tsv.gz');meta=rd(H/'stage2_scope_sample_counts.tsv');rows=[];datasets={}
for (scope,t,readout,variant,seed),dd in tf[tf.scope.isin(SC)].groupby(['scope','tf','readout','variant','seed']):
 base=dd[['sample','group','series','cells','mean_depth','value']].rename(columns={'cells':'receiver_cells','mean_depth':'receiver_depth','value':'y'})
 for sender in ['MacAM','MonoMac','MacWide']:
  for lig in L:
   sx=gene[(gene.key==sender)&(gene.gene==lig)][['sample','cells','mean_depth','mean_CP10k']].rename(columns={'cells':'sender_cells','mean_depth':'sender_depth','mean_CP10k':'x'});d=complete(base.merge(sx,on='sample'),True);tag=dict(receiver=scope,tf=t,readout=readout,variant=variant,seed=seed,sender=sender,ligand=lig);datasets[tuple(tag.values())]=d
   for step in steps:rows.append(dict(tag,step=step,**cor(d,step)))
up=bh(pd.DataFrame(rows),['readout','variant','seed','step']);wr(up,'all_candidate_TF_stepwise_models.tsv.gz');wr(up[up.P<.05],'all_nominal_candidate_TF_models.tsv')
# Validate final nested model against preceding round exactly, same n and rho/P.
prev=rd(V/'fixed_TF_bridge_associations.tsv.gz');prev=prev[prev.link=='ligand_to_TF'];keys=['receiver','tf','readout','variant','seed','sender','ligand'];audit=up[up.step=='SourceBothDepth'].merge(prev,on=keys,suffixes=('_new','_old'));assert len(audit)==len(prev);assert np.allclose(audit.rho_new,audit.rho_old,equal_nan=True,atol=1e-10);assert np.allclose(audit.P_new,audit.P_old,equal_nan=True,rtol=1e-8);wr(audit[keys+['n_new','n_old','rho_new','rho_old','P_new','P_old']],'previous_final_link_reproduction.tsv.gz')
trace=up.pivot(index=keys,columns='step',values=['rho','P','FDR','n','df']);trace.columns=['_'.join(x) for x in trace.columns];trace=trace.reset_index();trace['any_nominal']=trace[[c for c in trace if c.startswith('P_')]].lt(.05).any(axis=1);wr(trace,'candidate_TF_stepwise_traces.tsv.gz');wr(trace[trace.any_nominal],'candidate_TF_any_nominal_traces.tsv')
# TF->program and matching triplets for every prior final nominal upstream lead, no best seed selection.
pr=[];chains=[]
for scope,test in [(s,'Nominal05') for s in SC]+[('CD4-GZMA','FDR05')]:
 for (t,readout,variant,seed),dd in tf[tf.scope==scope].groupby(['tf','readout','variant','seed']):
  pg=prog[(prog.receiver==scope)&(prog['filter']=='All')&(prog.model=='Source')&(prog.test==test)&(prog.direction=='Up')&(prog.excluded_tf==t)][['sample','target_score']]
  base=dd[['sample','group','series','cells','mean_depth','value']].rename(columns={'cells':'receiver_cells','mean_depth':'receiver_depth','value':'x'}).merge(pg,on='sample').rename(columns={'target_score':'y'});d=complete(base)
  tag=dict(receiver=scope,test=test,tf=t,readout=readout,variant=variant,seed=seed)
  for step in ['Unadjusted','Group','Source','SourceReceiverDepth']:pr.append(dict(tag,step=step,**cor(d,step)))
  if test!='Nominal05':continue
  leads=up[(up.receiver==scope)&(up.tf==t)&(up.readout==readout)&(up.variant==variant)&(up.seed==seed)&(up.step=='SourceBothDepth')&(up.P<.05)]
  for _,v in leads.iterrows():
   sx=gene[(gene.key==v.sender)&(gene.gene==v.ligand)][['sample','cells','mean_depth','mean_CP10k']].rename(columns={'cells':'sender_cells','mean_depth':'sender_depth','mean_CP10k':'ligand_x'});dd=base.merge(sx,on='sample');dd=dd[np.isfinite(dd.ligand_x)];dd=complete(dd,True)
   for link,x,y in [('TF_program','x','y'),('ligand_program','ligand_x','y')]:chains.append(dict(tag,sender=v.sender,ligand=v.ligand,upstream_rho=v.rho,upstream_P=v.P,upstream_family_FDR=v.FDR,link=link,**cor(dd,'SourceBothDepth',x,y)))
pr=bh(pd.DataFrame(pr),['readout','variant','seed','step']);wr(pr,'TF_program_stepwise_models.tsv');wr(chains,'nominal_upstream_same_roster_program_checks.tsv')
# IL15 observed attrition, with >0 observations separated from existing percentage thresholds.
joint=rd(V/'IL15_matched_donor_coverage.tsv.gz');joint=joint[(joint['filter']=='All')&(joint.gate==.1)].copy();flow=[];allpairs=[]
for sender in ['MacAM','MonoMac','MacWide']:
 sm=gene[(gene.key==sender)&gene.gene.isin(['IL15','IL15RA'])].pivot(index=['sample','group','series','cells','mean_depth'],columns='gene',values='pct').reset_index().rename(columns={'cells':'sender_cells','mean_depth':'sender_depth'})
 for scope in ['Primary']+SC+['CD4state']:
  rm=gene[(gene.key==scope+'_All')&gene.gene.isin(['IL2RB','IL2RG'])].pivot(index=['sample','cells','mean_depth'],columns='gene',values='pct').reset_index().rename(columns={'cells':'receiver_cells','mean_depth':'receiver_depth'});d=sm.merge(rm,on='sample');d['sender']=sender;d['receiver']=scope
  for group,v in d.groupby('group'):
   present=(v.sender_cells>0)&(v.receiver_cells>0);n20=(v.sender_cells>=20)&(v.receiver_cells>=20);observed=(v[['IL15','IL15RA','IL2RB','IL2RG']]>0).all(axis=1)
   for stage,ok in [('PairedCellsPositive',present),('BothEnds20',n20),('FourGenesObserved',n20&observed),('FourMarginal5pct',n20&(v[['IL15','IL15RA','IL2RB','IL2RG']]>=.05).all(axis=1)),('FourMarginal10pct',n20&(v[['IL15','IL15RA','IL2RB','IL2RG']]>=.1).all(axis=1))]:flow.append(dict(sender=sender,receiver=scope,group=group,stage=stage,n=int(ok.sum()),samples=';'.join(v.loc[ok,'sample']),at_least_two=int(ok.sum())>=2))
  allpairs.append(d)
wr(flow,'IL15_fixed_condition_attrition.tsv');wr(pd.concat(allpairs),'IL15_all_continuous_paired_donors.tsv.gz')
measure=meas[meas.gene.isin(['IL15','IL15RA','IL2RB','IL2RG'])].copy();n=measure.cells.to_numpy();p=measure.pct.to_numpy();zz=1.95996398454;cen=(p+zz*zz/(2*n))/(1+zz*zz/n);half=zz*np.sqrt(p*(1-p)/n+zz*zz/(4*n*n))/(1+zz*zz/n);measure['wilson_low']=cen-half;measure['wilson_high']=cen+half;wr(measure,'IL15_measurements_with_descriptive_sampling_intervals.tsv.gz')
ms=measure[measure.cells>=20].groupby(['scope','gene','group','series']).agg(donors=('sample','size'),cells=('cells','sum'),counts=('counts','sum'),positive_cells=('detected','sum'),zero_donors=('detected',lambda x:sum(x==0)),median_pct=('pct','median'),min_pct=('pct','min'),max_pct=('pct','max'),median_depth=('mean_depth','median')).reset_index();wr(ms,'IL15_source_measurement_summary.tsv')
quart=rd(R/'within_donor_depth_quartiles.tsv.gz');qt=quart.pivot(index=['lineage','scope','sample','group','series','gene','cells'],columns='quartile',values='pct').reset_index();qt['high_minus_low']=qt[4]-qt[1];wr(qt,'within_donor_depth_detection_contrasts.tsv.gz')
# Independent program recomputations and frozen-control validation.
ind=rd(R/'independence_program_donors.tsv.gz');m2=rd(P/'11_reports/technical_delivery_20261002_ikzf2_stat1_v1/revision_v2/final_release/source_tables/M2_human_sample_scores.tsv');m2=m2[m2['filter']=='All'].copy();m2['scope']=m2.scope.str.replace('Subtype:','',regex=False)
check=ind[ind.program=='TReg_original'].merge(m2[m2.scope=='CD4-GZMA'][['sample','score_TReg']],on='sample');ok=check.score.notna()&check.score_TReg.notna();err=np.max(abs(check.loc[ok,'score']-check.loc[ok,'score_TReg']));assert err<1e-6,err
tree=ast.parse((P/'11_reports/human_alignment_20261005_v1/code/align_compute.py').read_text());exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='ols'],type_ignores=[]),'<frozen OLS>','exec'),globals())
disease=[]
for program,g in ind[ind.scope=='CD4-GZMA'].groupby('program'):
 for roster,cut in [('PositiveCells',1),('Cells20',20)]:
  d=g[g.cells>=cut].copy();d['cells']=np.maximum(d.cells,20);d['value']=d.score
  for model in ['Base','Source','SourceDep','WithinG227']:
   q=d[d.series=='G227'] if model=='WithinG227' else d;disease.append(dict(program=program,roster=roster,model=model,**ols(q,'Base' if model=='WithinG227' else model)))
disease=bh(pd.DataFrame(disease),['roster','model']);wr(disease,'GZMA_TReg_independence_disease_models.tsv')
# IL7R program composition and single-covariate controls, finite, no network rerun.
context=meas[(meas.scope=='CD4-IL7R')&(meas.gene=='IKZF2')][['sample','n_genes_mean','mt_mean']].merge(m2[m2.scope=='CD4-IL7R'][['sample','score_Naive','score_TAct']],on='sample');ir=[]
for (variant,seed),dd in tf[(tf.scope=='CD4-IL7R')&(tf.tf=='IKZF2')&(tf.readout=='AUC')].groupby(['variant','seed']):
 base=dd[['sample','group','series','cells','mean_depth','value']].rename(columns={'cells':'receiver_cells','mean_depth':'receiver_depth','value':'x'}).merge(context,on='sample')
 for program,g in ind[ind.scope=='CD4-IL7R'].groupby('program'):
  d=base.merge(g[['sample','score','members_n']],on='sample').rename(columns={'score':'y'});d=complete(d)
  for adjust in ['Baseline','n_genes_mean','mt_mean','score_Naive','score_TAct','WithinG227']:
   q=d[d.series=='G227'] if adjust=='WithinG227' else d;extra=None if adjust in ['Baseline','WithinG227'] else adjust;ir.append(dict(program=program,variant=variant,seed=seed,adjust=adjust,members_n=g.members_n.iloc[0],**cor(q,'SourceReceiverDepth',extra=extra)))
ir=bh(pd.DataFrame(ir),['variant','seed','adjust']);wr(ir,'IL7R_independence_state_diagnostics.tsv')
summary=dict(upstream_models=len(up),upstream_nominal=int((up.P<.05).sum()),upstream_final_nominal=int(((up.P<.05)&(up.step=='SourceBothDepth')).sum()),upstream_final_BH=int(((up.FDR<.05)&(up.step=='SourceBothDepth')).sum()),program_models=len(pr),TReg_reconstruction_max_error=float(err),same_roster_checks=len(chains));(O/'provenance/compute_complete.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2));print('FLOW_CORE_IL7R',pd.DataFrame(flow)[lambda x:x.receiver.isin(['Core','CD4-IL7R'])&x.group.eq('CWP')].drop(columns=['samples']).to_string(index=False));print('TREG_MINUS_SELF',disease.to_string(index=False));print('IL7R_DIAGNOSTICS',ir[['program','variant','seed','adjust','n','df','rho','P','FDR','status']].to_string(index=False));print('FINAL_NOMINAL_RNA',up[(up.step=='SourceBothDepth')&(up.readout=='RNA')&(up.P<.05)].to_string(index=False))
