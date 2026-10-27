from pathlib import Path
import os,json,hashlib,shutil
os.environ['OPENBLAS_NUM_THREADS']='1'
import pandas as pd,numpy as np
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/independent_candidates_revision_20261005_v1';R=O/'results';D=R/'source_evidence';H=P/'11_reports/human_first_20261004_v1';V=P/'11_reports/human_alignment_20261005_v1';T=P/'11_reports/technical_delivery_20261002_ikzf2_stat1_v1/revision_v2/final_release'
pool='IL15 TNF IKZF2 STAT1 CLSTN1 ADAM15 PLAU DSC2 CXCL2 ICAM2 VTN VEGFB GSTP1 RTN4 IL18'.split()
manifest=[]
for root,names in [(H/'results',['stage1_donor_design.tsv','stage1_complete_program_members.tsv','stage2_scope_sample_counts.tsv','stage2_DE_status.tsv','stage2_targetset_status.tsv','stage2_targetset_members.tsv.gz','stage2_background_members.tsv.gz','stage2_receptor_complex_audit.tsv.gz','stage2_prior_dimensions.tsv']), (V/'results',['previous_nine_influence_summary.tsv','TF_disease_target_overlap.tsv','TF_disease_overlap_complete_members.tsv.gz']), (P/'11_reports/human_gap_diagnostics_20261005_v1/results',['GZMA_IKZF2_source_composition_summary.tsv']), (T/'data',['S2_regulon_exclusion_coverage.tsv']), (T/'source_tables',['S4_external_time_matched_phenotype_deltas.tsv'])]:
 for name in names:
  f=root/name
  if not f.exists():print('OPTIONAL_MISSING',f);continue
  dest=D/name;shutil.copyfile(f,dest);manifest.append(dict(path=str(f),copy=str(dest.relative_to(O)),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
for name in ['stage2.R','stage2_followup.py','STAGE2_SPEC.txt']:
 f=H/'code'/name
 if not f.exists():f=H/'provenance'/name
 if f.exists():shutil.copyfile(f,O/'provenance'/('REUSED_'+name))
lt=[]
for f in sorted((H/'results').glob('stage2_LT_*.tsv.gz')):
 q=pd.read_csv(f,sep='\t');q=q[q.ligand.isin(pool)];q['source_file']=str(f);lt.append(q)
pd.concat(lt,ignore_index=True).to_csv(R/'fixed15_ligand_target_complete_members.tsv.gz',sep='\t',index=False)
q=pd.read_csv(H/'results/stage2_all_LR_supported_and_unsupported.tsv.gz',sep='\t');q=q[q['from'].isin(pool)];q.to_csv(R/'fixed15_LR_all_positive_negative.tsv.gz',sep='\t',index=False)
# Restore source-level statuses and exact P/FDR; no new significance family fitted to selected candidates.
def rd(n):return pd.read_csv(next(D.glob(n+'.tsv*')),sep='\t',low_memory=False)
q=rd('human_ligand_program');q=q[q.ligand.isin(pool)];keys=['receiver','filter','model','test','direction'];allq=rd('human_ligand_program');sz=allq.groupby(keys).P.count().rename('BH_m').reset_index();q=q.merge(sz,on=keys);q['BH_family_definition']='all finite ligand x sender tests within receiver/filter/DE-model/target-threshold/direction';q['association_adjustment']='group+source+rank_sender_depth+rank_receiver_depth';q.to_csv(R/'independent_ligand_program_ALL_RESULTS.tsv.gz',sep='\t',index=False)
q=rd('IL18_TNF_direct_phenotype');q=q[q.edge.isin(['ligand_phenotype','ligand_proximal'])].drop_duplicates(['signature']);q['interpretation']='INDEPENDENT_LIGAND_ENDPOINT;UPSTREAM_TF_LINK_NOT_REQUIRED';q.to_csv(R/'independent_TNF_IL18_phenotype_ALL_RESULTS.tsv.gz',sep='\t',index=False)
q=rd('human_TF_program_bridge');q=q[q.link=='TF_to_program'];q.to_csv(R/'independent_TF_disease_program_ALL_RESULTS.tsv.gz',sep='\t',index=False)
q=rd('mouse_ligand_activity');q=q[q.human_ligand.isin(pool)];q.to_csv(R/'independent_mouse_NicheNet_ALL_RESULTS.tsv',sep='\t',index=False)
q=rd('old_TF_program');q=q[q.tf.isin(['IKZF2','STAT1'])];q.to_csv(R/'independent_TF_old_program_ALL_RESULTS.tsv.gz',sep='\t',index=False)
q=rd('human_program_disease');q=q[(q.kind=='program')&(q['filter']=='All')&q.scope.isin(['CD4-CCR7','CD4-GZMA','CD4-IL7R','Treg','CD4state'])];q.to_csv(R/'all_subtype_disease_program_context.tsv',sep='\t',index=False)
print('INFLUENCE',pd.read_csv(D/'previous_nine_influence_summary.tsv',sep='\t').to_string(index=False))
print('IL18_LOCALIZATION_EXISTS',(R/'IL18_localization_all_author_scopes.tsv').exists())
print('MOUSE_SPEC',rd('mouse_ligand_activity').groupby(['receiver','filter','gate','test','direction','omit']).size().to_string())
print('EXTRA_MANIFEST',len(manifest));pd.DataFrame(manifest).to_csv(O/'provenance/additional_source_manifest.tsv',sep='\t',index=False)
# Original metadata/definitions remain a full server-side source trace.
print('ASSEMBLED',flush=True)
