from pathlib import Path
import os,json,pickle,hashlib,datetime
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np,pandas as pd
from scipy import stats
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/breakpoint_TNF_IL18_20261004_v1';R=O/'results';F=P/'11_reports/front_half_20261004_v1';M=P/'11_reports/missing_links_20261003_v1'
def rd(p):return pd.read_csv(p,sep='\t',low_memory=False)
def wr(d,n):pd.DataFrame(d).to_csv(R/n,sep='\t',index=False,na_rep='NA')
m=rd(R/'aligned_chain_models.tsv');inputs=rd(R/'aligned_model_inputs.tsv.gz');pair=rd(F/'results/paired_inputs.tsv')
with (O/'provenance/frozen_plan.pkl').open('rb') as f:cache,specs,data=pickle.load(f)
assert len(m)==1368 and int(m.calculation.eq('new_missing_signature').sum())==850
assert m.spec_id.is_unique and m.n.max()==19
index={k:v.sort_values('sample') for k,v in inputs.groupby('spec_id')};checks=[]
for _,v in m.iterrows():
 z=index[v.spec_id];assert len(z)==v.n and ';'.join(z['sample'])==v.sample_ids
 if v.calculation=='reused_historical':
  old=cache[v.signature];assert abs(old['rho']-v.rho)<1e-12 and abs(old['P']-v.P)<1e-12;continue
 if v.calculation!='new_missing_signature':continue
 C=[np.ones(len(z))];cs=[] if pd.isna(v.covariates) else v.covariates.split('+')
 for c in cs:
  if c=='G':C.append(z.group.eq('CWP').to_numpy(float))
  elif c=='S':C.extend([z.series.eq(s).to_numpy(float) for s in sorted(z.series.unique())[1:]])
  else:C.append(stats.rankdata(z[c]))
 C=np.column_stack(C);x=stats.rankdata(z.x_value);y=stats.rankdata(z.y_value);xx=x-C@np.linalg.lstsq(C,x,rcond=None)[0];yy=y-C@np.linalg.lstsq(C,y,rcond=None)[0]
 rho=np.corrcoef(xx,yy)[0,1];X=np.column_stack([C,x]);b=np.linalg.lstsq(X,y,rcond=None)[0];df=len(z)-np.linalg.matrix_rank(X);se=np.sqrt(np.sum((y-X@b)**2)/df*np.linalg.inv(X.T@X)[-1,-1]);p=2*stats.t.sf(abs(b[-1])/se,df)
 assert abs(rho-v.rho)<1e-9 and abs(p-v.P)<1e-8,(v.spec_id,rho,v.rho,p,v.P)
 checks.append(dict(spec_id=v.spec_id,rho_difference=abs(rho-v.rho),P_difference=abs(p-v.P)))
for _,q in m[m.policy.eq('aligned_steps')].groupby(['chain','modality']):assert q.sample_ids.nunique()==1
for _,q in m.groupby('FDR_family'):
 u=q[q.P.notna()].drop_duplicates('signature');assert np.allclose(stats.false_discovery_control(u.P,method='bh'),u.FDR)
wr(checks,'independent_new_fit_checks.tsv')
# Shared-roster sign compatibility, no mediation estimate.
signs=[]
for key,q in m[m.policy.isin(['aligned_steps','cells20','Treg_composition','TNK_competition'])].groupby(['chain','modality','policy','covariates'],dropna=False):
 lt=q[q.edge.eq('ligand_TF')];tf=q[q.edge.eq('TF_phenotype')&q.y.eq('State_Effector')];lp=q[q.edge.eq('ligand_phenotype')&q.y.eq('State_Effector')]
 if any(len(v)!=1 for v in [lt,tf,lp]):continue
 lt,tf,lp=[v.iloc[0] for v in [lt,tf,lp]];assert lt.sample_ids==tf.sample_ids==lp.sample_ids
 row=dict(zip(['chain','modality','policy','covariates'],key),n=lt.n,sample_ids=lt.sample_ids,LT_spec=lt.spec_id,TFP_spec=tf.spec_id,LP_spec=lp.spec_id,LT_rho=lt.rho,LT_P=lt.P,LT_FDR=lt.FDR,TFP_rho=tf.rho,TFP_P=tf.P,TFP_FDR=tf.FDR,LP_rho=lp.rho,LP_P=lp.P,LP_FDR=lp.FDR,product_sign=int(np.sign(lt.rho)*np.sign(tf.rho)),observed_LP_sign=int(np.sign(lp.rho)),sign_compatible=bool(np.sign(lt.rho)*np.sign(tf.rho)==np.sign(lp.rho)),all_three_nominal=bool(max(lt.P,tf.P,lp.P)<.05),all_three_FDR05=bool(max(lt.FDR,tf.FDR,lp.FDR)<.05),interpretation='sign consistency only; no indirect effect or causal mediation estimate')
 signs.append(row)
signs=pd.DataFrame(signs);wr(signs,'same_donor_chain_sign_audit.tsv')
near=m[m.edge.eq('ligand_proximal')];wr(near,'proximal_response_all_specs.tsv')
focal=m[(m.policy=='aligned_steps')&(m.covariates=='G+S+receiver_library+sender_library')];wr(focal,'common_final_chain_edges.tsv')
# Detection and zeros on every fixed roster; don't gate on expression.
gene=rd(M/'results/gene_coverage_by_sample.tsv');gene=gene[(gene.species=='Human')&(gene.role=='receiver')&gene.gene.isin(['TNFRSF1A','TNFRSF1B','IL18R1','IL18RAP'])];wr(gene,'receiver_component_detection_reused.tsv')
senderold=rd(P/'11_reports/technical_delivery_20261002_ikzf2_stat1_v1/revision_v2/final_release/prior_review_repaired/scenic3_review/results/reviewed_sender_expression.tsv');senderold=senderold[(senderold.sender_filter=='All')&senderold.ligand.isin(['TNF','IL18'])];wr(senderold,'sender_detection_reused.tsv')
coverage=[];perdonor=[]
rosters=rd(R/'chain_roster_summary.tsv')
for _,v in rosters.iterrows():
 sender,scope,axis_tf=v.chain.split('/');axis=axis_tf.split('-')[0];tf=axis_tf.split('-')[1];ids=v.sample_ids.split(';');z=pair[(pair.sender==sender)&(pair.scope==scope)&pair['sample'].isin(ids)]
 recs=['R_TNFR1','R_TNFR2'] if axis=='TNF' else ['R_IL18']
 for r in recs:
  coverage.append(dict(chain=v.chain,roster=v.roster,receptor=r,n=len(z),ligand_positive_subjects=int(z['L_'+axis].gt(0).sum()),receptor_codetected_subjects=int(z[r].gt(0).sum()),receptor_zero_subjects=int(z[r].eq(0).sum()),receptor_pct_median=z[r].median(),receptor_positive_cells=int(np.rint((z[r]*z.receiver_cells/100).sum())),total_receiver_cells=int(z.receiver_cells.sum()),ligand_unique_values=z['L_'+axis].nunique(),TF_RNA_zero_subjects=int(z[tf+'_RNA'].eq(0).sum()),active_ligand_assayed=False,functional_receptor_assayed=False))
  for _,s in z.iterrows():perdonor.append(dict(chain=v.chain,roster=v.roster,sample=s['sample'],group=s.group,series=s.series,sender_cells=s.sender_cells,receiver_cells=s.receiver_cells,ligand_RNA=s['L_'+axis],receptor=r,receptor_pct=s[r],receptor_positive_cells=int(round(s[r]*s.receiver_cells/100))))
wr(coverage,'chain_detection_summary.tsv');wr(perdonor,'chain_detection_by_donor.tsv')
for name in ['module_coverage.tsv','regulon_overlap_audit.tsv']:wr(rd(M/'results'/name)[lambda d:d.species.eq('Human')],'reused_'+name)
wr(rd(F/'results/panel_definition_overlap.tsv'),'phenotype_panel_overlap.tsv')
wr(rd(F/'results/treg_original_fraction_residual_side_by_side.tsv'),'treg_definitions_reused.tsv')
legacy=rd(F/'results/reused_reviewed_sender_TF.tsv');legacy=legacy[(legacy['filter']=='All')&legacy.ligand.isin(['TNF','IL18'])&legacy.scope.isin(['Core','CD4-IL7R'])];wr(legacy,'legacy_reviewed_sender_TF_context.tsv')
competition=rd(M/'results/all_new_models_final.tsv');competition=competition[(competition.species=='Human')&((competition.axis=='IFNG')|((competition.axis=='STAT1')&competition.y.str.startswith('IFN_state')))];wr(competition,'existing_IFNG_STAT1_IFN_context.tsv')
manifest=rd(R/'input_manifest.tsv');recheck=[]
for _,r in manifest.iterrows():
 h=hashlib.sha256()
 with Path(r.path).open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 recheck.append(dict(path=r.path,unchanged=h.hexdigest()==r.sha256))
assert all(r['unchanged'] for r in recheck);wr(recheck,'input_recheck.tsv')
verification=dict(status='PASS',new_unique_fits_verified=len(checks),historical_statistic_copies_verified=int(m.calculation.eq('reused_historical').sum()),max_rho_difference=max(v['rho_difference'] for v in checks),max_P_difference=max(v['P_difference'] for v in checks),fixed_chain_rosters_verified=True,BH_verified=True,input_files_unchanged=len(recheck),same_donor_sign_rows=len(signs))
(O/'provenance/verification.json').write_text(json.dumps(verification,indent=2));print(json.dumps(verification,indent=2))
print('\nMACAM_SIGN_BASE_FINAL');print(signs[signs.chain.str.startswith('MacAM/')&signs.policy.eq('aligned_steps')&(signs.covariates.isna()|signs.covariates.eq('G+S+receiver_library+sender_library'))].drop(columns=['sample_ids','interpretation','LT_spec','TFP_spec','LP_spec']).to_string(index=False))
print('\nMACAM_KEY_TRACE');print(m[m.sender.eq('MacAM')&m.modality.eq('AUC17')&m.edge.eq('ligand_TF')&~m.policy.eq('seed_sensitivity')][['chain','policy','covariates','n','rho','P','FDR','df','endpoint_VIF','max_covariate_VIF','LOO_min','LOO_max']].to_string(index=False))
print('\nMACAM_PROXIMAL');print(focal[focal.sender.eq('MacAM')&focal.edge.eq('ligand_proximal')][['scope','axis','modality','y','n','rho','P','FDR']].to_string(index=False))
print('\nCOVERAGE');print(pd.DataFrame(coverage)[lambda d:d.chain.str.startswith('MacAM/')].to_string(index=False))
print('\nSOURCE_FOCAL');within=rd(R/'within_source_descriptions.tsv');print(within[within.chain.str.startswith('MacAM/')&within.modality.eq('AUC17')&within.edge.isin(['ligand_TF','TF_phenotype'])&within.source.eq('G227')][['chain','edge','covariates','n','n_case','rho','P','status','LOO_min','LOO_max']].to_string(index=False))
