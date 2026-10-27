from pathlib import Path
import json,datetime,hashlib
import pandas as pd,numpy as np
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/missing_links_20261003_v1';R=O/'results'
def rd(p):return pd.read_csv(p,sep='\t')
r=rd(R/'receiver_paired_modalities.tsv');s=rd(R/'sender_all_samples.tsv');a=rd(R/'pairing_audit.tsv');g=rd(R/'gene_coverage_by_sample.tsv');c=rd(R/'module_coverage.tsv')
old=rd(P/'11_reports/stepwise_diagnostic_20261002_v1/inputs/receiver_all_positive_cell_samples.tsv')
old=old[old['filter'].eq('All')&old.seed.eq(17)&old.tf.eq('IKZF2')]
z=r.merge(old,on=['species','scope','sample'],suffixes=('','_old'),validate='one_to_one')
assert z.receiver_cells.eq(z.receiver_cells_old).all()
assert np.allclose(z.receiver_library,z.receiver_library_old,atol=1e-8)
assert np.allclose(z.IKZF2_AUC_17,z.mean_auc,atol=1e-12)
verified_receiver_n=len(z)
mouse=rd(P/'11_reports/mouse_sender_20261002_v1/results/frozen_AMlike_two_ligand_samples.tsv')
checks=[]
for ligand in ['IL15','TNF']:
 q=s[s.species.eq('Mouse')&s.sender.eq('AM')].merge(mouse[mouse.ligand.eq(ligand)],on='sample',suffixes=('','_old'),validate='one_to_one')
 assert q.sender_cells.eq(q.sender_cells_old).all() and np.allclose(q['L_'+ligand],q.sender_mean_CP10k,atol=1e-10)
 checks.append(dict(check='Mouse_AM_'+ligand,rows=len(q),max_abs_diff=float(abs(q['L_'+ligand]-q.sender_mean_CP10k).max())))
human=rd(P/'11_reports/technical_delivery_20261002_ikzf2_stat1_v1/revision_v2/final_release/prior_review_repaired/scenic3_review/results/reviewed_sender_expression.tsv')
for ligand in ['IL15','TNF','IL18','CXCL16']:
 q=s[s.species.eq('Human')&s.sender.eq('AM')].merge(human[human.sender.eq('MacAM')&human.sender_filter.eq('All')&human.ligand.eq(ligand)],on='sample',suffixes=('','_old'),validate='one_to_one')
 assert q.sender_cells.eq(q.sender_cells_old).all() and np.allclose(q['L_'+ligand],q.sender_mean_CP10k,atol=1e-10)
 checks.append(dict(check='Human_AM_'+ligand,rows=len(q),max_abs_diff=float(abs(q['L_'+ligand]-q.sender_mean_CP10k).max())))
pd.DataFrame(checks).to_csv(R/'frozen_expression_alignment_checks.tsv',sep='\t',index=False)
joint=[]
for (sp,scope),q in r.groupby(['species','scope']):
 for sender in ['AM','TNK_other']:
  z=q.merge(s[s.species.eq(sp)&s.sender.eq(sender)],on='sample',suffixes=('','_sender'),validate='one_to_one')
  for ax in ['IL15','TNFR1','TNFR2','IL18','CXCL16','IL27','IFNB1','IFNG']:
   lg='TNF' if ax.startswith('TNFR') else ax
   col='L_'+lg+('_complex' if lg in ['IL15','IL27'] else '')
   joint.append(dict(species=sp,scope=scope,sender=sender,axis=ax,paired_donors=len(z),ligand_detected_donors=int(z[col].gt(0).sum()),receptor_codetected_donors=int(z['R_'+ax].gt(0).sum()),both_sample_detected=int((z[col].gt(0)&z['R_'+ax].gt(0)).sum()),same_sender_cell_complex_donors=int(z['joint_'+lg].gt(0).sum()) if lg in ['IL15','IL27'] else np.nan,interpretation='RNA only; no expression gate applied to broad models'))
pd.DataFrame(joint).to_csv(R/'topology_pairing_summary.tsv',sep='\t',index=False,na_rep='NA')
summary=dict(status='PASS',time=datetime.datetime.now().isoformat(),receiver_rows=len(r),sender_rows=len(s),frozen_receiver_matched_rows=verified_receiver_n,frozen_cell_depth_auc_alignment=True,expression_alignment_checks=checks,IFNB1_zero_all_senders=bool(s.L_IFNB1.eq(0).all()),IL27_same_cell_codetection_zero=bool(s.joint_IL27.eq(0).all()),AUC_missing_donors=['CTL070','CTL103'],sender_without_receiver=['Human CTL105 TNK_other'],low_coverage_modules=c[~c.eligible][['species','module','variant','used']].to_dict('records'),proceed='finite broad RNA associations; IFNB1 constant models retained; IL27 incomplete sender codetection is a limitation, not a ligand-activity claim')
(O/'provenance/coverage_review.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
