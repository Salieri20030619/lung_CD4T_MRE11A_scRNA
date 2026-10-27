from pathlib import Path
import pandas as pd,json,hashlib,shutil,numpy as np
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');B=P/'11_reports';T=B/'technical_delivery_20261002_ikzf2_stat1_v1/revision_v2/final_release';H=B/'human_first_20261004_v1/results';S=B/'independent_candidates_revision_20261005_v1/results/source_evidence';O=B/'phase_delivery_IKZF2_STAT1_20261006_v2';out={};manifest=[]
def rd(f):
 f=Path(f);manifest.append(dict(path=str(f),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()));return pd.read_csv(f,sep='\t',low_memory=False)
def put(k,d):
 d=d.drop(columns=[c for c in ['sample','samples','sample_ids','donor','omitted','most_changed_donor'] if c in d],errors='ignore');out[k]=json.loads(d.to_json(orient='records',double_precision=15))
for tag in ['human_TF_disease','old_RNA_DE','old_TF_group','old_TF_program','IL7R_IKZF2_excess','mouse_TF_identity','IL7R_IKZF2_influence','GZMA_TReg_independence']:
 d=rd(next(S.glob(tag+'.tsv*')))
 if tag=='human_TF_disease':d=d[d['filter'].eq('All')]
 if tag=='old_TF_group':d=d[(d['filter'].isin(['All','BG1']))&d.scope.isin(['Core','ActCD4','CD4-GZMA','CD4-IL7R'])&d.variant.isin(['Original','DropOverlap'])&d.model.isin(['Group','Depth','Source','SourceDep'])]
 if tag=='old_TF_program':d=d[(d['filter']=='All')&d.scope.isin(['Core','CD4-GZMA','CD4-IL7R'])&(d.variant=='DropOverlap')&(d.model=='SourceDep')&d.feature.isin(['score_Inflam','score_PyroExec','score_TReg','score_Inhib'])]
 put(tag,d)
for name in ['S2_regulon_exclusion_coverage','TF_network_all_edges','network_display_edges']:
 d=rd(T/'data'/f'{name}.tsv')
 if name=='TF_network_all_edges':d=d.groupby(['species','filter','seed','tf']).agg(targets=('target','nunique')).reset_index()
 put(name,d)
for name in ['M3_geneset_status','M5_geneset_status','M5_focused_ligand_ranks','S4_transfer_coverage','S4_direction_replication','Human_preliminary_lineage_summary','M1_module_definition']:
 put(name,rd(T/'source_tables'/f'{name}.tsv'))
put('stage1_models',rd(H/'stage1_all_finite_models.tsv'))
for name in ['stage2_targetset_status.tsv','stage2_DE_status.tsv','stage2_prior_dimensions.tsv']:
 put(name.split('.')[0],rd(H/name))
d=rd(H/'stage2_scope_sample_counts.tsv');put('scope_coverage',d.groupby(['lineage','scope','filter','group']).agg(cells=('cells','sum'),positive_donors=('cells',lambda x:int(x.gt(0).sum())),donors20=('cells',lambda x:int(x.ge(20).sum()))).reset_index())
d=rd(H/'stage1_pairing_all_sources.tsv');put('pairing',d.groupby(['receiver','sender','group']).agg(paired_donors=('paired20','sum'),receiver_cells=('receiver_cells','sum'),sender_cells=('sender_cells','sum')).reset_index())
d=rd(H/'stage2_ligand_activities.tsv.gz');put('human_NN_main',d[(d['filter']=='All')&((d['rank']<=5)|d.test_ligand.isin(['IL15','TNF','IL18','ICAM2','VEGFB']))])
d=rd(H/'stage2_evidence_all_candidates.tsv.gz');put('human_NN_fixed',d[(d['filter']=='All')&d.ligand.isin(['IL15','TNF','IL18','ICAM2','VEGFB'])]);put('human_NN_inventory',d.groupby(['filter','status'],dropna=False).size().reset_index(name='rows'))
d=rd(H/'stage2_same_donor_associations.tsv.gz');put('human_NN_BH',d[d.FDR<.05])
d=rd(B/'human_alignment_20261005_v1/results/previous_nine_influence_summary.tsv');put('NN_BH_influence',d)
d=rd(H/'stage3_reused_mouse_ligand_activities.tsv');put('mouse_NN',d[(d['filter'].isin(['All','BG1']))&d.human_ligand.isin(['IL15','TNF','IL18','CD86','CD274','CLSTN1','ADAM15','VEGFB'])])
lt=[]
for rec in ['Core','CD4-GZMA','CD4-IL7R','Treg']:
 f=H/f'stage2_LT_{rec}_All_Source_Nominal05_Up.tsv.gz'
 if f.exists():
  d=rd(f);d=d[d.ligand.isin(['IL15','TNF','ICAM2','VEGFB'])];lt.append(d.sort_values('weight',ascending=False).groupby('ligand').head(5))
put('human_LT_examples',pd.concat(lt,ignore_index=True))
br=B/'bootstrap_TNF_IL18_20261004_v1/results'
for name in ['effect_summary','paired_difference_summary','chain_summary','all_seeds_chain_summary','diagnostic_summary','original_diagnostics']:
 put('bootstrap_'+name,rd(br/f'{name}.tsv'))
put('human_legacy_edgeR_audit',rd(B/'original_investigation_audit_20261002/results/human_original_TF_refit_comparison.tsv'))
d=rd(B/'original_investigation_audit_20261002/results/human_original_edgeR_roster.tsv');print('ROSTER_COLS',list(d));print(d.groupby([c for c in ['analysis','scope','Group','group'] if c in d]).size() if any(c in d for c in ['analysis','scope']) else d.shape)
for name in ['Core_Full_focal_refit_vs_original','Core_NoBG_focal_refit_vs_original']:
 put('original_'+name,rd(B/f'original_investigation_audit_20261002/results/{name}.tsv'))
summary=[]
for sp in ['Human','Mouse']:
 for arm in ['All','BG1']:
  root=P/f'02_analysis/scenic1_20260929_103513_xRoEnr/{sp}/{arm}';cells=rd(root/'cells.tsv');runs=rd(root/'run_summary.tsv')
  for _,r in runs.iterrows():summary.append(dict(species=sp,filter=arm,donors=cells['sample'].nunique(),**r.to_dict()))
put('SCENIC_runs',pd.DataFrame(summary))
for f in (B/'bootstrap_TNF_IL18_20261004_v1/deliverables').glob('*.png'):shutil.copy2(f,O/'figures'/f.name)
for name in ['F05_regulon_activity_effects.png','F06_regulon_network_display.png','F07_human_module_associations.png','F07_mouse_module_associations.png']:
 f=T/'figures'/name
 if f.exists():shutil.copy2(f,O/'figures'/name)
(O/'provenance/report_evidence_aggregate.json').write_text(json.dumps(out,ensure_ascii=False,indent=2));(O/'provenance/report_source_manifest.json').write_text(json.dumps(manifest,indent=2));print('EXTRACTED',len(out),{k:len(v) for k,v in out.items()})
