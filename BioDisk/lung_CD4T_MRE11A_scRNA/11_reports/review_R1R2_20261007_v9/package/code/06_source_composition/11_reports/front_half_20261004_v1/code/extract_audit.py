from pathlib import Path
import os,json,hashlib,fcntl,datetime,sys,itertools
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import numpy as np,pandas as pd
from scipy.io import mmread
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/front_half_20261004_v1';A=P/'02_analysis'
RV=P/'11_reports/technical_delivery_20261002_ikzf2_stat1_v1/revision_v2/final_release/prior_review_repaired/scenic3_review/results'
M=P/'11_reports/missing_links_20261003_v1/results';S=A/'scenic3_20261001_123452_Kvp9h5';H=A/'mechanism_step2_20260927_154919_EsN4gk/input'
lock=open(O/'analysis.lock','a+');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
assert not (O/'provenance/extraction_complete.json').exists()
manifest={}
def track(p):
 p=Path(p)
 if str(p) not in manifest:
  h=hashlib.sha256()
  with p.open('rb') as f:
   for b in iter(lambda:f.read(2**20),b''):h.update(b)
  manifest[str(p)]={'path':str(p),'bytes':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns,'sha256':h.hexdigest()}
 return p
def rd(p):return pd.read_csv(track(p),sep='\t',low_memory=False)
def wr(d,n):pd.DataFrame(d).to_csv(O/'results'/n,sep='\t',index=False,na_rep='NA')
review=rd(RV/'reviewed_cells.tsv.gz');meta=rd(S/'tables/human_refined_cells.tsv.gz')
assert review.cell.is_unique and meta.cell.is_unique and set(review.cell)==set(meta.cell)
meta=meta.merge(review[['cell','review_scope','label']],on='cell',validate='one_to_one');mapped=meta[meta['sample'].ne('Unknown')&meta['sample'].notna()].copy()
annotation=rd(RV/'reviewed_annotation.tsv');wr(annotation,'annotation_frozen.tsv')
composition=mapped.groupby(['sample','group','series','review_scope']).size().unstack(fill_value=0).reset_index()
composition['candidate_myeloid_cells']=composition[['MacAM','MonoMac','Excluded']].sum(axis=1)
composition['MacWide']=composition.MacAM+composition.MonoMac
for k in ['MacAM','MonoMac','MacWide']:composition[k+'_fraction_candidate']=composition[k]/composition.candidate_myeloid_cells
composition['MacAM_fraction_MacWide']=composition.MacAM/composition.MacWide
wr(composition,'sender_composition.tsv')
dominance=[]
for cl,q in mapped.groupby('r4'):
 n=q['sample'].value_counts();dominance.append(dict(cluster=cl,cells=len(q),subjects=len(n),dominant_sample=n.index[0],dominant_fraction=n.iloc[0]/len(q),scope=q.review_scope.iloc[0],label=q.label.iloc[0]))
wr(dominance,'cluster_donor_dominance.tsv')
states=[]
for scope in ['MacAM','MonoMac','MacWide']:
 q=mapped[mapped.review_scope.isin(['MacAM','MonoMac']) if scope=='MacWide' else mapped.review_scope.eq(scope)]
 for ss,v in q.groupby('sample'):
  row=dict(sender=scope,sample=ss,group=v.group.iloc[0],series=v.series.iloc[0],sender_cells=len(v),sender_library=v.library.mean(),sender_epi=v.epi_fraction.mean())
  row.update({c:float(v[c].mean()) for c in ['AM','C1qMac','SPP1Mac','FOLR2Mac','Mono']});states.append(row)
states=pd.DataFrame(states)
sen=rd(RV/'reviewed_sender_expression.tsv');sen=sen[sen.sender_filter.eq('All')]
sw=sen.pivot(index=['sender','sample'],columns='ligand',values='sender_mean_CP10k').add_prefix('L_').reset_index()
states=states.merge(sw,on=['sender','sample'],validate='one_to_one').merge(composition[['sample','MacAM_fraction_MacWide']],on='sample',validate='many_to_one')
assert len(states)==57;wr(states,'sender_state_and_ligands.tsv')
# Exact reviewed counts recheck, not refitting old models.
cnt=sen[['sender','sample','sender_cells']].drop_duplicates().merge(states[['sender','sample','sender_cells']],on=['sender','sample'],suffixes=('_old','_new'),validate='one_to_one');assert (cnt.sender_cells_old==cnt.sender_cells_new).all();wr(cnt,'sender_count_recheck.tsv')
print('SOURCE_AUDIT_DONE',flush=True)
cells=rd(H/'human_cells.tsv');genes=rd(H/'human_genes.tsv').symbol.str.upper().to_numpy();x=mmread(track(H/'human_counts.mtx')).T.tocsr()
assert x.shape==(len(cells),len(genes)) and cells.cell.is_unique and (cells.library>0).all()
assert np.isfinite(x.data).all() and (x.data>=0).all() and np.allclose(x.data,np.rint(x.data))
take=cells.core.astype(str).str.lower().eq('true').to_numpy();cells=cells[take].reset_index(drop=True);x=x[take];lib=cells.library.to_numpy()
panels={'Treg':'FOXP3 IL2RA CTLA4 IKZF4 LRRC32 ENTPD1 TIGIT','Activation':'CD69 ICOS TNFRSF4 BATF IRF4 NR4A1 NR4A2 EGR1','Naive':'CCR7 SELL TCF7 LEF1 KLF2 MAL','Inhibitory':'PDCD1 HAVCR2 TOX TOX2','Effector':'GZMA GZMB PRF1 IFNG NKG7'}
panels={k:v.split() for k,v in panels.items()};assert len(sum(panels.values(),[]))==len(set(sum(panels.values(),[])))
edges=rd(A/'scenic1_20260929_103513_xRoEnr/Human/All/regulon_edges.tsv');union=set(edges.loc[edges.tf.str.upper().isin(['IKZF2','STAT1']),'target'].str.upper())|{'IKZF2','STAT1'}
scores={};definition=[];coverage=[]
for name,gs in panels.items():
 post=[g for g in gs if g not in union];vals=[];used=[]
 for g in gs:
  loc=np.flatnonzero(genes==g);v=np.log1p(np.asarray(x[:,loc].sum(axis=1)).ravel()/lib*1e4) if len(loc) else np.full(len(cells),np.nan)
  good=g in post and np.isfinite(v).all() and v.std()>0
  if good:vals.append((v-v.mean())/v.std());used.append(g)
  definition.append(dict(panel=name,gene=g,in_regulon_union=g in union,present=len(loc)>0,used=good))
 eligible=len(used)>=2 and len(used)>=.6*len(post)
 scores['State_'+name]=np.mean(vals,axis=0) if eligible else np.full(len(cells),np.nan)
 coverage.append(dict(panel=name,original_n=len(gs),post_removal_n=len(post),used_n=len(used),eligible=eligible,used_genes=';'.join(used),removed_regulon=';'.join(g for g in gs if g in union)))
wr(definition,'panel_definition_overlap.tsv');wr(coverage,'panel_coverage.tsv')
receiver=[]
for scope in ['Core','CD4-IL7R']:
 q=cells if scope=='Core' else cells[cells.subtype.eq(scope)]
 for ss,v in q.groupby('sample'):
  row=dict(species='Human',scope=scope,sample=ss,group=v.group.iloc[0],series=v.series.iloc[0],receiver_cells=len(v),receiver_library=v.library.mean())
  row.update({c:float(val[v.index].mean()) for c,val in scores.items()})
  receiver.append(row)
receiver=pd.DataFrame(receiver)
old=rd(M/'receiver_paired_modalities.tsv');old=old[old.species.eq('Human')]
keys=['species','scope','sample'];over=receiver.merge(old[keys+['receiver_cells','receiver_library']],on=keys,suffixes=('_new','_old'),validate='one_to_one');assert (over.receiver_cells_new==over.receiver_cells_old).all() and np.allclose(over.receiver_library_new,over.receiver_library_old)
cols=['group','series','receiver_cells','receiver_library'];receiver=receiver.merge(old.drop(columns=cols),on=keys,validate='one_to_one')
fr=cells.groupby(['sample','subtype']).size().unstack(fill_value=0);fr=fr.div(fr.sum(axis=1),axis=0).add_prefix('Fraction_').reset_index();receiver=receiver.merge(fr,on='sample',validate='many_to_one');wr(receiver,'receiver_audit_inputs.tsv')
wr(cells.groupby(['sample','group','series','subtype']).size().reset_index(name='cells'),'receiver_subtype_counts.tsv')
wr(cells[['cell','sample','group','series','subtype']],'receiver_cell_roster.tsv.gz')
tnk=rd(M/'sender_all_samples.tsv');tnk=tnk[tnk.species.eq('Human')&tnk.sender.eq('TNK_other')][['sample','L_IFNG']].rename(columns={'L_IFNG':'TNK_IFNG'})
paired=receiver.merge(states.drop(columns=['group','series']),on='sample',validate='many_to_many').merge(tnk,on='sample',validate='many_to_one');wr(paired,'paired_inputs.tsv')
wr(receiver[['sample','group','series']].drop_duplicates().groupby(['series','group']).size().reset_index(name='subjects'),'source_group_design.tsv')
print('RECEIVER_AUDIT_DONE',len(receiver),flush=True)
# Reuse old evidence exactly; all legacy P/FDR retained.
legacy=rd(RV/'reviewed_sender_TF_models.tsv');legacy=legacy[legacy.tf.isin(['IKZF2(+)','STAT1(+)'])&legacy.scope.isin(['Core','CD4-IL7R'])];wr(legacy,'reused_reviewed_sender_TF.tsv')
aligned=rd(P/'11_reports/aligned_evidence_20261002_v1/results/aligned_all_models.tsv.gz');q=aligned[aligned.species.eq('Human')&aligned['filter'].eq('All')&aligned.scope.isin(['Core','CD4-IL7R'])&aligned.evidence.eq('TF_MODULE')&aligned.feature.isin(['TAct','TReg','Naive','Inhib'])];wr(q,'reused_aligned_TF_states.tsv')
oldmodels=rd(M/'all_new_models_final.tsv');q=oldmodels[oldmodels.axis.isin(['STAT1','IKZF2','IFNG'])&~oldmodels.y.str.contains('OXPHOS',na=False)];wr(q,'reused_TF_state_and_TNK_competition.tsv')
for stage in ['mechanism_step3_20260928_103057_lIRL38','mechanism_step5_20260928_153825_Q4sB1x']:
 d=rd(A/stage/'tables/ligand_activities.tsv');q=d[d.test_ligand.str.upper().isin(['IL15','TNF','IL18','CXCL16'])];wr(q,'reused_'+stage.split('_')[1]+'_nichenet_ranks.tsv')
lr=rd(A/'mechanism_step3_20260928_103057_lIRL38/tables/supported_LR_edges.tsv');wr(lr[lr['from'].str.upper().isin(['IL15','TNF','IL18','CXCL16'])],'reused_mouse_LR.tsv')
target=rd(A/'mechanism_step5_20260928_153825_Q4sB1x/tables/geneset_members.tsv');target=target[(target['filter']=='All')&(target.gate=='Pct10')&(target.omit=='Full')];overlaps=[]
for key,v in target.groupby(['receiver','test','direction']):
 gs=set(v.gene.str.upper())
 for name,ps in panels.items():overlaps.append(dict(receiver=key[0],test=key[1],direction=key[2],panel=name,geneset_n=len(gs),overlap=';'.join(sorted(gs&set(ps))),overlap_n=len(gs&set(ps)),interpretation='mouse disease geneset overlap only; not ligand-specific targets or human validation'))
wr(overlaps,'mouse_target_program_overlap.tsv')
wr(rd(P/'11_reports/original_investigation_audit_20261002/results/human_original_TF_refit_comparison.tsv'),'reused_human_original_TF_RNA.tsv')
wr(rd(P/'11_reports/mouse_sender_20261002_v1/results/mouse_sender_all_models.tsv'),'reused_mouse_sender_models.tsv')
wr(manifest.values(),'input_manifest.tsv')
(O/'provenance/extraction_complete.json').write_text(json.dumps(dict(time=datetime.datetime.now().isoformat(),human_subjects=receiver['sample'].nunique(),receiver_rows=len(receiver),sender_rows=len(states),paired_rows=len(paired),panels=coverage,no_old_models_rerun=True),indent=2))
print('EXTRACTION_COMPLETE',json.dumps(coverage),flush=True)
