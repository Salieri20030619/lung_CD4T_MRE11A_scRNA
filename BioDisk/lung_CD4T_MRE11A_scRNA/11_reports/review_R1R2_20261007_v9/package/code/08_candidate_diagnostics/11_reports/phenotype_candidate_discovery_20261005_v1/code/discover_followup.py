from pathlib import Path
import os,json,datetime,hashlib
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import pandas as pd,numpy as np
from scipy import stats
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/phenotype_candidate_discovery_20261005_v1';R=O/'results';H=P/'11_reports/human_first_20261004_v1/results'
# Import the frozen pure statistical helpers without executing the discovery pipeline.
src=(O/'code/discover_compute.py').read_text();helper=src[src.index('def wr('):src.index("scores=rd(")];exec(helper)
rd=lambda p:pd.read_csv(p,sep='\t',low_memory=False)
spec=O/'code/DISCOVERY_SHORTLIST_SPEC.txt';(O/'provenance/DISCOVERY_SHORTLIST_SPEC.txt').write_bytes(spec.read_bytes());(O/'provenance/shortlist_freeze.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'spec_sha256':hashlib.sha256(spec.read_bytes()).hexdigest(),'genes':['ALCAM','HBEGF','HMGB1']}))
sel=pd.DataFrame([dict(ligand='ALCAM',receiver='CD4-IL7R',metric='score_Mitophagy'),dict(ligand='HBEGF',receiver='CD4-IL7R',metric='score_Mitophagy'),dict(ligand='HMGB1',receiver='CD4-GZMA',metric='score_TReg')]);sel['sender']='MacAM';wr(sel,'shortlist_frozen.tsv')
inputs=rd(R/'direct_phenotype_donor_inputs_SERVER_ONLY.tsv.gz');scores=rd(H/'stage1_reused_all_program_inputs.tsv');scores['receiver']=scores.scope.str.replace('Subtype:','',regex=False);expr=rd(H/'stage2_LR_expression_all_donors.tsv.gz');meta=rd(H/'stage2_scope_sample_counts.tsv');rows=[];loo=[];disease=[];donors=[];strat=[]
for tag in sel.to_dict('records'):
 v=inputs.copy()
 for k,val in tag.items():v=v[v[k]==val]
 v=v[v.kind=='program'].copy();assert v['sample'].is_unique
 bg=scores[(scores.receiver==tag['receiver'])&(scores['filter']=='BG1')&scores.cells.ge(20)][['sample','cells','library',tag['metric']]].rename(columns={'cells':'receiver_cells','library':'receiver_depth',tag['metric']:'y'})
 bg=v.drop(columns=['receiver_cells','receiver_depth','y']).merge(bg,on='sample');common=set(bg['sample']);allcom=v[v['sample'].isin(common)]
 for label,z in [('All_eligible',v),('All_fixed_BG1_roster',allcom),('BG1_fixed_same_roster',bg),('BG1_eligible',bg)]:
  donors.append(z.assign(followup_roster=label))
  for step in ['Unadjusted','Group','Source','SourceBothDepth']:rows.append(dict(tag,roster=label,step=step,**corr(z,step)))
  zz=z.rename(columns={'receiver_depth':'library'})
  for model in ['Base','Source','SourceDep','WithinG227']:disease.append(dict(tag,roster=label,model=model,**ols(zz,'y',model)))
 for source,z in v.groupby('series'):
  for step in ['Unadjusted','Group']:strat.append(dict(tag,source=source,step=step,**corr(z,step)))
 for step in ['Unadjusted','Group','Source','SourceBothDepth']:
  full=corr(v,step)
  for sample in v['sample']:
   z=v[v['sample']!=sample];r=corr(z,step);loo.append(dict(tag,step=step,omitted=sample,baseline_rho=full['rho'],same_sign=np.sign(full['rho'])==np.sign(r['rho']) if np.isfinite(r['rho']) else np.nan,**r))
models=bh(pd.DataFrame(rows),['roster','step']);# additionally one entire followup family, duplicate BG1 aliases not counted twice
models['FDR_followup_unique']=np.nan;ix=models.index[(models.roster!='BG1_eligible')&models.P.notna()];models.loc[ix,'FDR_followup_unique']=stats.false_discovery_control(models.loc[ix,'P']);wr(models,'selected_fixed_roster_associations.tsv');dy=bh(pd.DataFrame(disease),['roster','model']);wr(dy,'selected_paired_roster_disease_models.tsv');wr(loo,'selected_donor_leave_one_out_SERVER_ONLY.tsv');wr(strat,'selected_source_stratified_models.tsv');wr(pd.concat(donors),'selected_fixed_donor_inputs_SERVER_ONLY.tsv.gz')
ll=pd.DataFrame(loo);summary=ll.groupby(['ligand','receiver','metric','step']).agg(n_omit=('omitted','size'),n_finite=('rho','count'),same_sign=('same_sign',lambda x:int(pd.to_numeric(x,errors='coerce').sum())),rho_min=('rho','min'),rho_max=('rho','max'),P_min=('P','min'),P_max=('P','max')).reset_index();wr(summary,'selected_influence_summary.tsv')
# Preserve complete old context for the three selected directions and all old genes.
e=rd(H/'stage2_evidence_all_candidates.tsv.gz');se=e[e.ligand.isin(sel.ligand)];wr(se,'selected_all_existing_specifications.tsv.gz')
lr=rd(H/'stage2_all_LR_supported_and_unsupported.tsv.gz');wr(lr[lr['from'].isin(sel.ligand)],'selected_LR_all_supported_and_unsupported.tsv.gz')
genes=set(sel.ligand)|set(lr.loc[lr['from'].isin(sel.ligand),'to']);cv=expr[expr.gene.isin(genes)];wr(cv,'selected_LR_donor_expression_SERVER_ONLY.tsv.gz');cov=cv.groupby(['key','gene','group']).apply(lambda z:pd.Series({'donors20':int(z.cells.ge(20).sum()),'observed_donors20':int((z.cells.ge(20)&z.pct.gt(0)).sum()),'pct10_donors20':int((z.cells.ge(20)&z.pct.ge(.1)).sum()),'median_detection20':z.loc[z.cells.ge(20),'pct'].median(),'total_cells20':z.loc[z.cells.ge(20),'cells'].sum()}),include_groups=False).reset_index();wr(cov,'selected_expression_coverage_summary.tsv')
senderde=rd(H/'stage2_sender_disease_DE.tsv.gz');wr(senderde[senderde.gene.isin(sel.ligand)],'selected_sender_disease_DE.tsv')
members=rd(H/'stage1_complete_program_members.tsv');targets=rd(H/'stage2_targetset_members.tsv.gz');background=rd(H/'stage2_background_members.tsv.gz');targetrows=[]
for tag in sel.to_dict('records'):
 gs=members[members.program==tag['metric'].replace('score_','')]
 for filt in ['All','BG1']:
  for model in ['Base','Source']:
   de=rd(H/f"stage2_DE_{tag['receiver']}_{filt}_{model}.tsv.gz");d=gs.merge(de[['gene','logFC','PValue','FDR']],on='gene',how='left');d['receiver']=tag['receiver'];d['ligand']=tag['ligand'];d['filter']=filt;d['model']=model
   bb=set(background[(background.receiver==tag['receiver'])&(background['filter']==filt)&(background.model==model)].gene);d['in_original_background']=d.gene.isin(bb)
   for test in ['Nominal05','FDR05']:
    for direction in ['Up','Down']:
     st=set(targets[(targets.receiver==tag['receiver'])&(targets['filter']==filt)&(targets.model==model)&(targets.test==test)&(targets.direction==direction)].gene);d[test+'_'+direction]=d.gene.isin(st)
   targetrows.append(d)
wr(pd.concat(targetrows),'selected_complete_program_DE_members.tsv')
lt=[]
for f in H.glob('stage2_LT_*.tsv.gz'):
 d=rd(f);d=d[d.ligand.isin(sel.ligand)]
 if len(d):lt.append(d)
wr(pd.concat(lt),'selected_complete_existing_ligand_target_members.tsv.gz');wr(background[background.receiver.isin(sel.receiver)],'selected_complete_background_members.tsv.gz')
mouse=rd(H/'stage3_reused_mouse_ligand_activities.tsv');wr(mouse[mouse.human_ligand.isin(sel.ligand)],'selected_mouse_existing_activities.tsv');mapping=rd(H/'stage3_mouse_human_orthologs.tsv');wr(mapping[mapping.human.isin(genes|set(members.gene))],'selected_mouse_human_orthologs.tsv')
mp=rd(P/'02_analysis/mechanism_step1_20260927_152218_s4V5eF/tables/group_effects.tsv');wr(mp[mp.metric.isin(sel.metric)&mp.dataset.isin(['GSE314943','GSE183682'])],'selected_mouse_program_directions.tsv')
ann=rd(H/'stage1_annotation_frozen.tsv');dom=rd(H/'stage1_cluster_donor_dominance.tsv');wr(ann,'reviewed_sender_annotation.tsv');wr(dom,'reviewed_cluster_donor_dominance.tsv')
trace=rd(R/'all_gene_endpoint_decision_audit.tsv.gz');trace['shortlist_status']='not_selected_after_full_pool_review';trace.loc[trace.previously_studied_15,'shortlist_status']='previous_direction_preserved_not_new_discovery'
for tag in sel.to_dict('records'):
 mask=np.ones(len(trace),bool)
 for k,val in tag.items():mask&=trace[k].eq(val)
 trace.loc[mask,'shortlist_status']='selected_independent_hypothesis'
trace.loc[trace.ligand.eq('NUCB2')&trace.receiver.eq('CD4-IL7R')&trace.metric.eq('score_Mitophagy'),'shortlist_status']='reserve_nominal_state_correlate_ERAP1_prior_not_established_signaling_receptor'
wr(trace,'all_gene_endpoint_final_decisions.tsv.gz')
print('FIXED_MODELS',models[['ligand','roster','step','n','rho','P','FDR_followup_unique']].to_string(index=False));print('INFLUENCE',summary.to_string(index=False));print('PAIRED_PHENOTYPE',dy[(dy.roster=='All_eligible')][['ligand','model','n','effect','P','FDR']].to_string(index=False));(O/'provenance/followup_complete.json').write_text(json.dumps({'directions':3,'fixed_models':len(models),'loo_models':len(ll),'no_new_network':True},indent=2))
