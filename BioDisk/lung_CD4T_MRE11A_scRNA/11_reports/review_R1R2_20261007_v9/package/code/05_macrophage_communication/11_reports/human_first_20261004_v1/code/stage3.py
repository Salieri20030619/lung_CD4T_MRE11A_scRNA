from pathlib import Path
import pandas as pd,numpy as np,json,shutil
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/human_first_20261004_v1';R=O/'results';A=P/'02_analysis';M3=A/'mechanism_step3_20260928_103057_lIRL38/tables';M5=A/'mechanism_step5_20260928_153825_Q4sB1x/tables'
def rd(p):return pd.read_csv(p,sep='\t')
def wr(d,n):pd.DataFrame(d).to_csv(R/n,sep='\t',index=False,na_rep='NA')
mapping=rd(R/'stage3_mouse_human_orthologs.tsv');mp=mapping.dropna(subset=['human']);mp=mp[~mp.human.duplicated(False)&~mp.mouse.duplicated(False)];mm=dict(zip(mp.mouse,mp.human))
mouse=rd(M5/'ligand_activities.tsv');mouse=mouse[(mouse.gate=='Pct10')&(mouse.omit=='Full')].copy();mouse['human_ligand']=mouse.test_ligand.map(mm);wr(mouse,'stage3_reused_mouse_ligand_activities.tsv')
sets=rd(M5/'geneset_members.tsv');sets=sets[(sets.gate=='Pct10')&(sets.omit=='Full')].copy();sets['human_gene']=sets.gene.map(mm);wr(sets,'stage3_reused_mouse_target_members.tsv')
status=rd(M5/'geneset_status.tsv');wr(status[(status.gate=='Pct10')&(status.omit=='Full')],'stage3_mouse_targetset_status.tsv')
phen=rd(A/'mechanism_step1_20260927_152218_s4V5eF/tables/group_effects.tsv');phen=phen[phen.dataset.isin(['GSE314943','GSE183682'])&phen.metric.isin(['score_TAct','score_TReg','score_Inhib','score_Naive'])];wr(phen,'stage3_reused_mouse_program_directions.tsv')
ha=rd(R/'stage2_ligand_activities.tsv.gz');ha=ha.rename(columns={'test_ligand':'human_ligand','receiver':'human_receiver','rank':'human_rank','aupr_corrected':'human_aupr_corrected'})
mo=mouse.rename(columns={'receiver':'mouse_receiver','rank':'mouse_rank','aupr_corrected':'mouse_aupr_corrected'})
cross=ha.assign(_cross=1).merge(pd.DataFrame({'mouse_receiver':['Core','ActCD4'],'_cross':[1,1]}),on='_cross').drop(columns='_cross').merge(mo[['human_ligand','mouse_receiver','filter','test','direction','mouse_rank','mouse_aupr_corrected','test_ligand']],on=['human_ligand','mouse_receiver','filter','test','direction'],how='left')
cross['alignment_class']=np.where((cross.human_receiver=='Core')&(cross.mouse_receiver=='Core'),'mixed_CD4_context_comparison','NO_ESTABLISHED_ONE_TO_ONE_SUBTYPE_ALIGNMENT')
cross['positive_in_both']=cross.human_aupr_corrected.gt(0)&cross.mouse_aupr_corrected.gt(0);cross['comparison_status']=np.where(cross.mouse_aupr_corrected.isna(),'MOUSE_MATCHED_ACTIVITY_UNAVAILABLE_SEE_MOUSE_STATUS',np.where(cross.positive_in_both,'BOTH_ACTIVITIES_POSITIVE','AT_LEAST_ONE_COMPUTED_ACTIVITY_NONPOSITIVE'));cross.loc[cross.mouse_aupr_corrected.isna(),'positive_in_both']=pd.NA;cross['prior_caveat']='human_v1_vs_mouse_v2; ranks not numerically comparable; not replication';wr(cross,'stage3_all_ligand_context_comparisons.tsv.gz')
hsets=rd(R/'stage2_targetset_members.tsv.gz');hbg=rd(R/'stage2_background_members.tsv.gz');mbg=rd(R/'stage3_mouse_background_members.tsv');mbg['human_gene']=mbg.gene.map(mm)
ov=[]
for key,v in hsets.groupby(['receiver','filter','model','test','direction']):
 rec,filt,model,test,direction=key;hg=set(v.gene);hb=set(hbg[(hbg.receiver==rec)&(hbg['filter']==filt)&(hbg.model==model)].gene)
 for mr in ['Core','ActCD4']:
  mg=set(sets[(sets.receiver==mr)&(sets['filter']==filt)&(sets.test==test)&(sets.direction==direction)].human_gene.dropna());mb=set(mbg[(mbg.receiver==mr)&(mbg['filter']==filt)].human_gene.dropna());common=hb&mb;hi=hg&common;mi=mg&common
  ov.append(dict(human_receiver=rec,mouse_receiver=mr,filter=filt,model=model,test=test,direction=direction,common_background=len(common),human_targets_common=len(hi),mouse_targets_common=len(mi),same_direction_overlap=len(hi&mi),overlap_members=';'.join(sorted(hi&mi)),interpretation='directional context; no one-to-one subtype validation'))
wr(ov,'stage3_target_program_overlap.tsv')
# All mapped tested-gene directions, not restricted to favorable candidate genes.
summ=[];joined=[]
for hp in sorted(R.glob('stage2_DE_*.tsv.gz')):
 h=rd(hp);hr=h.receiver.iloc[0];filt=h['filter'].iloc[0];model=h.model.iloc[0]
 for mr in ['Core','ActCD4','Naive']:
  mpth=M3/f'DE_CD4_{mr}_{filt}.tsv'
  if not mpth.exists():continue
  m=rd(mpth);m['human_gene']=m.gene.map(mm);m=m.dropna(subset=['human_gene']);q=h.merge(m[['gene','human_gene','logFC','PValue','FDR']],left_on='gene',right_on='human_gene',suffixes=('_human','_mouse'));q['mouse_receiver']=mr;q['same_direction']=np.sign(q.logFC_human)==np.sign(q.logFC_mouse)
  q['alignment_class']='mixed_CD4_context' if hr=='Core' and mr=='Core' else 'naive_annotation_context' if hr=='CD4-CCR7' and mr=='Naive' else 'not_established_subtype_equivalence'
  common=q[(q.PValue_human<.05)&(q.PValue_mouse<.05)&(q.logFC_human.abs()>=.5)&(q.logFC_mouse.abs()>=.5)]
  summ.append(dict(human_receiver=hr,mouse_receiver=mr,filter=filt,model=model,tested_overlap=len(q),both_nominal=len(common),both_nominal_same_direction=int(common.same_direction.sum()),both_FDR05=int(((q.FDR_human<.05)&(q.FDR_mouse<.05)).sum()),alignment_class=q.alignment_class.iloc[0] if len(q) else 'none'))
  joined.append(q)
wr(summ,'stage3_gene_direction_summary.tsv');wr(pd.concat(joined),'stage3_all_mapped_gene_directions.tsv.gz')
# No P/FDR invented for GSE183682 1+1 at each timepoint.
ext=phen[phen.dataset=='GSE183682'].copy();ext['inference']='one_library_per_group_per_time; direction_only';wr(ext,'stage3_GSE183682_direction_only.tsv')
(O/'provenance/stage3_complete.json').write_text(json.dumps({'mapping_unique_pairs':len(mp),'human_mouse_ligand_comparisons':len(cross),'gene_comparisons':sum(map(len,joined)),'no_mouse_reanalysis':True},indent=2));print('STAGE3_COMPLETE',len(cross),len(summ),flush=True)
