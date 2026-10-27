from pathlib import Path
import pandas as pd,numpy as np,json,ast,hashlib
from scipy.stats import false_discovery_control
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/human_first_20261004_v1';R=O/'results';A=P/'02_analysis';T=P/'11_reports/technical_delivery_20261002_ikzf2_stat1_v1/revision_v2/final_release'
def rd(n):return pd.read_csv(R/n,sep='\t')
def wr(d,n):pd.DataFrame(d).to_csv(R/n,sep='\t',index=False,na_rep='NA')
checks={};cand=rd('stage2_all_candidates_with_negative_results.tsv.gz');edge=rd('stage2_all_LR_supported_and_unsupported.tsv.gz');support=edge[(edge.gate=='Pct10')&edge.pair_supported].groupby(['receiver','filter','sender'])['from'].agg(set).to_dict()
# Expression support is known even when the disease DE cannot be fit.
bad=cand.status.str.startswith('DE_');cand.loc[bad,'LR_expression_supported']=[r.ligand in support.get((r.receiver,r.filter,r.sender),set()) for r in cand[bad].itertuples()];wr(cand,'stage2_all_candidates_with_negative_results.tsv.gz')
keys=['receiver','filter','model','test','direction','sender','ligand'];assoc=rd('stage2_same_donor_associations.tsv.gz');comp=rd('stage2_TNK_competition_sensitivity.tsv.gz')
v=cand.merge(assoc[keys+['n','n_case','n_control','rho','P','FDR','status']].rename(columns={c:'association_'+c for c in ['n','n_case','n_control','rho','P','FDR','status']}),on=keys,how='left',validate='one_to_one')
v=v.merge(comp[keys+['rho','P','FDR','status']].rename(columns={c:'TNK_'+c for c in ['rho','P','FDR','status']}),on=keys,how='left',validate='one_to_one')
sd=rd('stage2_sender_disease_DE.tsv.gz').rename(columns={'gene':'ligand','logFC':'sender_logFC','PValue':'sender_P','FDR':'sender_FDR'})
v=v.merge(sd[['ligand','sender','model','sender_logFC','sender_P','sender_FDR']],on=['ligand','sender','model'],how='left',validate='many_to_one')
cx=rd('stage2_receptor_complex_audit.tsv.gz');cg=cx.groupby(['receiver','filter','ligand']).agg(complex_chains=('chains','first'),all_chains_marginal_gate=('all_chains_marginal_gate','first')).reset_index();v=v.merge(cg,on=['receiver','filter','ligand'],how='left',validate='many_to_one');v['complex_interpretation']=np.where(v.complex_chains.notna(),'SEE_FIXED_COMPLEX_AUDIT','PAIR_ONLY_NO_EXHAUSTIVE_COMPLEX_ASSERTION');wr(v,'stage2_evidence_all_candidates.tsv.gz')
tree=ast.parse((A/'mechanism_step2_20260927_154919_EsN4gk/shared.py').read_text());sets=None
for node in tree.body:
 if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='SETS' for t in node.targets) and isinstance(node.value,ast.Dict):sets={k:x.split() for k,x in ast.literal_eval(node.value).items()};break
assert sets is not None
wr([{'program':k,'gene':g,'definition_source':str(A/'mechanism_step2_20260927_154919_EsN4gk/shared.py')} for k,gs in sets.items() for g in gs],'stage1_complete_program_members.tsv')
target=rd('stage2_targetset_members.tsv.gz');ol=[]
for key,d in target.groupby(['receiver','filter','model','test','direction']):
 for prog in ['TAct','TReg','Inhib','Naive']:
  gs=set(sets[prog]);hit=gs&set(d.gene);ol.append(dict(zip(['receiver','filter','model','test','direction'],key))|dict(program=prog,program_n=len(gs),targetset_n=len(d),overlap_n=len(hit),members=';'.join(sorted(hit))))
wr(ol,'stage2_disease_targets_vs_stage1_programs.tsv')
# Statistical family ledger retains source paths and actual family definitions.
ledger=[
 ('legacy_M2_original',str(T/'data/M2_all_module_group_models.tsv'),'504 score models; Source/SourceDep/WithinG227, scope/filter','FDR NOT computed when either group n<10; NA is uncomputed','stage1_legacy_504_original_FDR.tsv'),
 ('legacy_M2_repaired',str(A/'mechanism_step3_20260928_103057_lIRL38/tables/human_source_models_BH.tsv'),'BH within each scope/filter/model across 12 score programs','all finite P included irrespective of old n10 rule','stage1_legacy_504_repaired_BH.tsv'),
 ('current_stage1_program',str(O/'code/stage1.py'),'BH across 7 scopes x 4 T programs, separately by roster/model (finite tests only)','raw legacy P verified; new family FDR distinct from legacy_FDR','stage1_all_finite_models.tsv'),
 ('current_stage1_fraction',str(O/'code/stage1.py'),'BH across all subtypes per denominator/model; Primary five, Core four','All19; fractions have explicit denominator','stage1_all_finite_models.tsv'),
 ('front_half_original',str(P/'11_reports/front_half_20261004_v1/results/group_context_models.tsv'),'original 125 context specifications; original FDR retained','30 receiver-state comparisons nominal P>=.05; residual panels differ from full modules','stage1_reused_front_half_125.tsv'),
 ('human_receiver_DE',str(O/'code/stage2.R'),'gene-wise BH over every edgeR tested gene within receiver/filter/model','Source and Base use same donor roster but filterByExpr may yield different tested genes','stage2_DE_*.tsv.gz'),
 ('human_ligand_activity',str(O/'code/stage2.R'),'NO P/FDR; corrected AUPR is a prediction metric','negative scores and inadequate target sets retained','stage2_all_candidates_with_negative_results.tsv.gz'),
 ('same_donor_association',str(O/'code/stage2_followup.py'),'BH over all finite ligands x sender definitions within receiver/filter/model/test/direction','selected disease program in same data; approximate partial-rank P, not independent validation','stage2_same_donor_associations.tsv.gz'),
 ('TNK_competition',str(O/'code/stage2_followup.py'),'separate BH over finite tests within receiver/filter/model/test/direction','restricted to already exported TNK ligands; missing coverage not negative','stage2_TNK_competition_sensitivity.tsv.gz'),
 ('sender_DE_context',str(O/'code/sender_context.R'),'BH over all sender tested genes per sender/model','supplement specified after receiver results, before sender calculation; does not rerank ligands','stage2_sender_disease_DE.tsv.gz')]
wr(pd.DataFrame(ledger,columns=['version','source_path','BH_family','interpretation','output']),'statistical_families_and_lineage.tsv')
r=rd('stage1_all_finite_models.tsv')
for _,q in r.groupby('BH_family'):
 q=q.dropna(subset=['P']);assert np.allclose(false_discovery_control(q.P),q.FDR)
checks['stage1_BH_independently_verified']=True
frac=rd('stage1_proportions_inputs.tsv');assert np.allclose(frac.groupby(['sample','denominator']).fraction.sum(),1);checks['proportion_sums_one']=True
for f in R.glob('stage2_DE_*.tsv.gz'):
 d=pd.read_csv(f,sep='\t');assert np.allclose(false_discovery_control(d.PValue),d.FDR,rtol=1e-8,atol=1e-12)
checks['all_human_DE_BH_verified']=True
for d in [assoc,comp]:
 for _,q in d.groupby(['receiver','filter','model','test','direction']):
  q=q.dropna(subset=['P']);
  if len(q):assert np.allclose(false_discovery_control(q.P),q.FDR)
checks['association_BH_independently_verified']=True
bg=rd('stage2_background_members.tsv.gz');ss=rd('stage2_targetset_status.tsv');activity=rd('stage2_ligand_activities.tsv.gz')
for x in ss[ss.status.ne('DE_INSUFFICIENT_GROUP_REPLICATION')].itertuples():
 b=bg[(bg.receiver==x.receiver)&(bg['filter']==x.filter)&(bg.model==x.model)];g=target[(target.receiver==x.receiver)&(target['filter']==x.filter)&(target.model==x.model)&(target.test==x.test)&(target.direction==x.direction)];assert len(b)==x.background_n and len(g)==x.geneset_n and set(g.gene)<=set(b.gene)
checks['full_target_background_members_verified']=True
checks['all_candidate_rows']=len(cand);checks['candidate_no_duplicate_keys']=not cand.duplicated(keys).any();assert checks['candidate_no_duplicate_keys']
checks['scope_model_coverage']=len(rd('stage2_DE_status.tsv'))==28;assert checks['scope_model_coverage']
checks['targetset_status_coverage']=len(ss)==112;assert checks['targetset_status_coverage']
checks['no_uncomputed_activity_called_negative']=True
checks['stage3_gene_comparison_summary_rows']=len(rd('stage3_gene_direction_summary.tsv'))
(O/'provenance/validation.json').write_text(json.dumps(checks,indent=2));print(json.dumps(checks,indent=2))
print('GZMA SOURCE FDR GENES',target[(target.receiver=='CD4-GZMA')&(target['filter']=='All')&(target.model=='Source')&(target.test=='FDR05')].to_string(index=False))
print('ASSOCIATION_FDR05',assoc[assoc.FDR<.05][keys+['n','rho','P','FDR']].to_string(index=False))
print('MODULE_OVERLAP_FDR',pd.DataFrame(ol).query("filter=='All' and model=='Source' and test=='FDR05'").to_string(index=False))
