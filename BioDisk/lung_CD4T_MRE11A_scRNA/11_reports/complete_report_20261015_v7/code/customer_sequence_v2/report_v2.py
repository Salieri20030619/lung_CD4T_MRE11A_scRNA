from pathlib import Path
import sys,json,itertools,hashlib
import numpy as np,pandas as pd
from scipy.io import mmread
from scipy.stats import permutation_test,false_discovery_control
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(sys.argv[1]);P=O.parents[1];M=P/'02_analysis/mechanism_step1_20260927_152218_s4V5eF'
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
def read(n):return pd.read_csv(O/n,sep='\t')
def wr(a,n):a.to_csv(O/n,sep='\t',index=False,encoding='utf-8',na_rep='NA')
def save(fig,n):
 for ext in ['png','svg']:fig.savefig(O/'figures'/(n+'.'+ext),dpi=170,bbox_inches='tight',facecolor='white')
 plt.close(fig)
cands=pd.concat([read(f'04_candidate_intersection/{s}_primary_intersection.tsv') for s in ['Naive','ActCD4']],ignore_index=True)
pt=read('03A_state_pseudotime/genomewide_trajectory_tests.tsv.gz');sl=read('03A_state_pseudotime/sample_level_slopes.tsv.gz');counts=read('04_candidate_intersection/stepwise_gene_counts.tsv');dirx=read('05_supporting_evidence/GSE183682_candidate_directions.tsv');human=read('05_supporting_evidence/GSE174725_candidate_group_comparisons.tsv');aut=read('05_supporting_evidence/candidate_autophagy_association_summary.tsv')
loc=read('05_supporting_evidence/GSE296682_candidate_coverage_status.tsv')
# Independent SciPy validation of all trajectory P values; biological unit is library.
validation=[]
for (state,run),v in pt.groupby(['state','run']):
 matrix=sl[(sl.state==state)&(sl.run==run)&(sl.type=='gene')].pivot(index='sample',columns='feature',values='slope').reindex(columns=v.feature).to_numpy()
 check=permutation_test((matrix,),lambda a,axis:np.abs(np.mean(a,axis=axis)),permutation_type='samples',vectorized=True,n_resamples=np.inf,alternative='greater',axis=0).pvalue
 ok=v.P.notna().to_numpy();err=float(np.max(np.abs(check[ok]-v.loc[ok,'P'])))
 assert err<1e-10,(state,run,err)
 assert np.allclose(false_discovery_control(v.loc[ok,'P']),v.loc[ok,'BH'],rtol=1e-10)
 validation.append(dict(state=state,run=run,valid_tests=int(ok.sum()),max_P_difference=err))
for name in ['Omnibus','D7_D0','D28_D0','D56_D0']:
 v=read('03B_real_time/'+name+'.tsv');assert np.allclose(false_discovery_control(v.PValue),v.FDR,rtol=1e-10)
wr(pd.DataFrame(validation),'provenance/independent_trajectory_P_BH_validation.tsv')
# Evidence cards keep negative and unmeasured evidence.
cards=[]
for r in cands.itertuples():
 v=pt[(pt.state==r.state)&(pt.feature==r.canonical)];d=dirx[(dirx.scope==r.state)&(dirx.gene==r.canonical)].set_index('day');h=human[(human['filter']=='ReviewedAll')&(human.scope=='Core')&(human.gene==r.canonical)].iloc[0];a=aut[aut.gene==r.canonical].iloc[0]
 consistent=int((np.sign(v.mean_slope)==np.sign(r.DE_logFC)).sum());nom=int((v.P<.05).sum());direct=int((np.sign(d.difference)==np.sign(r.DE_logFC)).sum())
 priority='expression/state validation first' if nom==4 and direct==2 else 'resolve discordance before mechanism work'
 cards.append(dict(gene=r.canonical,mouse_state=r.state,DE_log2FC=r.DE_logFC,DE_P=r.DE_P,DE_BH=r.DE_BH,PT_slope=r.mean_slope,PT_P=r.PT_P,PT_BH=r.PT_BH,root_direction_consistent=f'{consistent}/4',root_nominal_P05=f'{nom}/4',time_omnibus_P=r.time_omnibus_P,time_omnibus_BH=r.time_omnibus_BH,nominal_same_direction_days=r.matched_days,external_D7_difference=d.loc[7,'difference'],external_D56_difference=d.loc[56,'difference'],human_Core_difference=h.difference,human_Core_P=h.P,human_Core_BH=h.BH,autophagy_partial_rho=a.mean_partial_rho,autophagy_P=a.P,autophagy_BH=a.BH,GO_autophagy_member=a.in_GO_autophagy,patient_localization_status=loc[loc.gene==r.canonical].status.iloc[0],experimental_stage=priority,evidence_level='exploratory nominal intersection; not a confirmed autophagy regulator'))
cards=pd.DataFrame(cards);wr(cards,'04_candidate_intersection/candidate_evidence_cards.tsv')
# Every prior favoured gene receives a reproducible disposition, including not-tested records.
old=['IL1R2','PENK','TNFRSF9','LAG3','TTC39C','FOXP1','IKZF2','STAT1','MRE11'];audit=[]
for state in ['Naive','ActCD4']:
 a=read(f'04_candidate_intersection/{state}_all_gene_decisions.tsv.gz')
 for gene in old:
  v=a[a.canonical==gene]
  if not len(v):audit.append(dict(state=state,gene=gene,selected=False,reason='not in the tested disease-DE universe'));continue
  v=v.iloc[0];failed=[]
  for field in ['pass_DE','pass_PT','PT_direction_match','pass_time_omnibus','time_nominal_direction_match']:
   if not bool(v[field]):failed.append(field)
  audit.append(dict(state=state,gene=gene,selected=bool(v.selected),reason='passes fixed primary rule' if v.selected else ';'.join(failed),DE_P=v.DE_P,DE_BH=v.DE_BH,PT_P=v.PT_P,PT_BH=v.PT_BH,time_omnibus_P=v.time_omnibus_P))
wr(pd.DataFrame(audit),'04_candidate_intersection/previous_candidates_disposition.tsv')
fig,axs=plt.subplots(2,1,figsize=(15,7))
labels=['Disease tested','Disease P<0.05','+ state PT P<0.05','+ PT direction','+ time omnibus P<0.05','+ time contrast/direction']
for ax,state in zip(axs,['Naive','ActCD4']):
 a=counts[counts.state==state];ax.axis('off')
 for j,(label,n) in enumerate(zip(labels,a.genes)):
  ax.text(.02+j*.165,.55,f'{label}\n{n:,}',transform=ax.transAxes,ha='left',va='center',fontsize=10,bbox=dict(boxstyle='round,pad=.6',fc='#e9f1f8',ec='#6888a4'))
  if j<5:ax.annotate('',xy=(.01+(j+1)*.165,.55),xytext=(.14+j*.165,.55),xycoords='axes fraction',arrowprops=dict(arrowstyle='->',lw=1))
 ax.set_title(state+' -> '+', '.join(cands[cands.state==state].canonical),fontsize=14,loc='left')
fig.suptitle('Prospective intersection: identical thresholds for every gene',fontsize=17);fig.tight_layout(rect=[0,0,1,.93]);save(fig,'04_candidate_screening_flow')
# Candidate state trajectories: all six independent sample curves.
ce=pd.read_csv(M/'comm_input/CD4_cells.tsv',sep='\t');gn=pd.read_csv(M/'comm_input/CD4_genes.tsv',sep='\t');x=mmread(M/'comm_input/CD4_counts.mtx.gz').T.tocsr();y=x.multiply((1e4/ce.library.to_numpy())[:,None]).tocsr();y.data=np.log1p(y.data);gidx=pd.Index(gn.canonical);p=read('03A_state_pseudotime/pseudotime_cell_membership.tsv.gz').set_index('index');rows=[]
fig,axs=plt.subplots(1,len(cands),figsize=(16,4.4))
for ax,r in zip(axs,cands.itertuples()):
 take=ce.subtype==r.state;cc=ce.loc[take,'cell'];t=p.loc[cc,'K15_R1'].to_numpy();t=(t-t.min())/(t.max()-t.min());z=y[take.to_numpy(),gidx.get_loc(r.canonical)].toarray().ravel();samples=ce.loc[take,'sample'].to_numpy()
 for sample in sorted(set(samples)):
  mask=samples==sample;bins=pd.qcut(pd.Series(t[mask]).rank(method='first'),4,labels=False).to_numpy();xx=[];yy=[]
  for k in range(4):
   xx.append(t[mask][bins==k].mean());yy.append(z[mask][bins==k].mean());rows.append(dict(state=r.state,gene=r.canonical,sample=sample,bin=k,cells=int((bins==k).sum()),mean_pt=xx[-1],mean_logCP10k=yy[-1]))
  ax.plot(xx,yy,'o-',ms=4,lw=1,alpha=.7,color='#2c6f9f' if sample.startswith('Saline') else '#c84d3e')
 ax.set(title=f'{r.state}: {r.canonical}',xlabel='Within-state DPT',ylabel='Mean log1p(CP10k)')
fig.suptitle('Candidate trajectories: sample curves are descriptive; P uses library slopes',fontsize=14);fig.tight_layout(rect=[0,0,1,.91]);save(fig,'03A_candidate_gene_trajectories');wr(pd.DataFrame(rows),'03A_state_pseudotime/candidate_sample_curves.tsv')
# Real-time gene curves retain each original sample, with equal-sample means.
time=read('03B_real_time/sample_TMM_logCPM.tsv');time['canonical']=time.symbol.str.upper();meta=read('03B_real_time/sample_inclusion.tsv')
fig,axs=plt.subplots(1,len(cands),figsize=(16,4.5))
for ax,r in zip(axs,cands.itertuples()):
 row=time[time.canonical==r.canonical].iloc[0];means=[]
 for j,day in enumerate([0,7,28,56]):
  ids=meta[meta.day==day]['sample'];vals=row[ids].astype(float).to_numpy();ax.scatter(j+np.linspace(-.07,.07,len(vals)),vals,color='#286c93',s=30);means.append(vals.mean())
 ax.plot(range(4),means,'-',color='#c44f3d');ax.set(title=r.canonical,xticks=range(4),xticklabels=['D0','D7','D28','D56'],xlabel='Measured day',ylabel='TMM log2CPM')
fig.suptitle('GSE274372 real time: all 12 TLike libraries; not a pseudotime axis',fontsize=14);fig.tight_layout(rect=[0,0,1,.91]);save(fig,'03B_candidate_real_time')
# Supplemental real-time program means, no absolute cross-cohort comparison.
s=pd.read_csv(M/'tables/sample_scores.tsv',sep='\t');s=s[(s.dataset=='GSE274372')&(s.scope=='All')&(s['filter']=='All')];assert len(s)==12
wr(s,'03B_real_time/all12_existing_program_scores.tsv')
fig,axs=plt.subplots(2,4,figsize=(15,7))
for ax,mod in zip(axs.flat,['PyroExec','Inflam','ApoExec','OXPHOS','MitoISR','MitoROS','Mitophagy']):
 means=[]
 for j,day in enumerate([0,7,28,56]):
  vals=s[s.day==day]['score_'+mod].to_numpy();ax.scatter(j+np.linspace(-.07,.07,len(vals)),vals,color='#286c93');means.append(vals.mean())
 ax.plot(range(4),means,color='#c44f3d');ax.set(title=mod,xticks=range(4),xticklabels=['D0','D7','D28','D56'],ylabel='Sample program score')
axs.flat[-1].axis('off');axs.flat[-1].text(.1,.7,'All samples shown.\nDescriptive real-time program curves.\nNo cell-count cutoff.\nNo comparison of absolute\nscores across cohorts.',transform=axs.flat[-1].transAxes)
fig.tight_layout();save(fig,'03B_program_real_time')
# Main report follows the customer sequence; full negatives and limitations are kept.
death=read('01_CD4_death/disease_program_comparisons.tsv');mito=read('02_mito_autophagy/disease_program_comparisons.tsv');progpt=read('03A_state_pseudotime/program_trajectory_tests.tsv');auto_cov=read('02_mito_autophagy/autophagy_coverage.tsv')
methods='''<p>All nonempty, measurable sample/state units are retained. No minimum-20 or minimum-10 cell rule is applied to the new main analyses. Cell-level quality/identity annotations and valid expression filtering remain unchanged. Cells and pseudotime bins are not independent biological samples. Every main mouse scope has all 3 Saline + 3 Silica libraries; the real-time model has all 12 libraries (three per day). Exact animal/pool-to-library mapping remains unresolved; inference is conditional on the historical replicate design.</p>
<p>Program comparisons use frozen M1 matched-expression scores, plus a new broad autophagy panel from installed GO:0006914 GOALL annotations (577 symbols; 553 measured), scored as mean log1p(CP10k) minus matched-control expression. Each library contributes one mean. Scores measure transcripts, not death, mitochondrial damage or autophagic flux. Two-sided exact group-label permutations use absolute mean difference. BH families separate Core/subtypes, death/mitochondrial-autophagy programs and background sensitivity. Core is a composition-sensitive aggregate, not a fourth independent subtype.</p>
<p>Existing four DPT orderings are reused unchanged. They cover all 1,734 Naive/ActCD4 cells; no cells were lost by the old graph coverage criterion. K15_R1 is fixed primary; the other three runs are sensitivity analyses. Each state is rescaled to its observed DPT range. Within each library, log1p(CP10k) is regressed on DPT with log1p library depth; the resulting six slopes are the units for exact two-sided sign-flip tests (64 allocations), with BH over all tested genes per state/run. Program trajectories have their own BH families. This is a test of reproducible within-state trends conditional on the inferred graph and a symmetric sample-slope null, not cell-level inference or proof of lineage/time. Treg is not forced onto the Naive/ActCD4 trajectory.</p>
<p>GSE274372 all-12 TLike pseudobulk uses edgeR robust quasi-likelihood, categorical day, full original library sizes and recalculated TMM. Gene expression filtering remains standard filterByExpr. Omnibus, D7-D0, D28-D0 and D56-D0 have separate complete tested-gene BH families. The old minimum-10 model remains unchanged and is not relabelled as the new model. TLike annotation is not confirmed pure CD4, and later days lack contemporaneous saline controls.</p>'''
selection='''<p>For each discovery state, require: disease DE raw P&lt;0.05; primary state-trajectory raw P&lt;0.05; trajectory slope matching the disease logFC direction; real-time omnibus raw P&lt;0.05; and at least one D7/D28/D56 versus D0 contrast with raw P&lt;0.05 in the same direction. No newly tuned fold-change cutoff or BH significance requirement is imposed. All model P and BH values, excluded genes and stepwise member lists are preserved. Real-time support is from a broader TLike scope, not a claim of identical cross-cohort state annotation.</p>'''
parts=['<!doctype html><html lang="en"><meta charset="utf-8"><title>CD4 customer-sequence revision v2</title><style>body{font:15px Arial;max-width:1450px;margin:35px auto;line-height:1.5;color:#162c3d;padding:0 18px}h1,h2,h3{color:#17496a}table{border-collapse:collapse;font-size:12px;display:block;overflow-x:auto;margin:12px 0}td,th{padding:7px;border:1px solid #c7d3df;text-align:left}img{max-width:100%;height:auto}.note{background:#edf3f7;padding:14px;border-left:4px solid #497b9d}</style><body>',
'<h1>Silicosis CD4 project: customer-sequence revision v2</h1><p class="note">Independent project revision, 2026-10-15. No sample cell-count threshold. Original reports and v1 preserved. This is an internal English working analysis, not an authorized customer translation.</p>',
'<h2>1. Compare pyroptosis and apoptosis within CD4 subpopulations</h2><p>All six libraries are retained in Core, Naive, ActCD4 and Treg. Naive and ActCD4 PyroExec increase (+0.044881 and +0.051628), while Naive ApoExec decreases (-0.007577) and ActCD4 ApoExec is essentially unchanged (+0.000056). With all six libraries retained, Treg ApoExec decreases (-0.010557). None of these group comparisons reaches nominal P&lt;0.05. This distinction prevents conflating pyroptosis with apoptosis or selecting a state solely by its absolute score.</p><img src="figures/01_death_programs.png">',death[death['filter']=='All'].to_html(index=False,float_format=lambda x:f'{x:.5g}'),
'<h2>2. Compare mitochondrial and autophagy programs in those states</h2><p>Naive and ActCD4 mitophagy-associated scores decrease (-0.048500 and -0.038711); Treg increases (+0.037061). Broad GO autophagy changes are small (Naive -0.002203; ActCD4 -0.002235; Treg +0.000406) and do not reach nominal significance. Current data do not establish broad autophagy activation, autophagy suppression or flux changes. MitoROS is an antioxidant-response transcription panel, not a direct ROS/damage measurement.</p><img src="figures/02_mito_autophagy_programs.png">',auto_cov.to_html(index=False),
'<h2>3A. Follow genes and autophagy programs within suitable CD4 state trajectories</h2><p>Separate within-state tests identify 668 nominal gene trends in Naive and 683 in ActCD4 under the fixed primary ordering. Broad autophagy trajectories are not supported (Naive P=1; ActCD4 P=0.75); mitophagy trajectories are also not nominally supported. A reproducible gene trajectory therefore does not by itself imply an autophagy mechanism.</p><img src="figures/03A_program_trajectories.png"><img src="figures/03A_candidate_gene_trajectories.png">',progpt[progpt.run=='K15_R1'][['state','feature','mean_slope','P','BH','n_libraries']].to_html(index=False,float_format=lambda x:f'{x:.5g}'),
'<h2>3B. Analyse measured mouse time separately</h2><p>All 12 GSE274372 TLike libraries are included, including the 8-cell D0 library. The model tests 732 genes: 55 have nominal omnibus time P&lt;0.05 (2 BH&lt;0.05); D7/D28/D56 contrasts yield 93/81/67 nominal genes. These are measured-time effects, not pseudotime. Differences from the old model reflect the full roster, corresponding expression filter and explicitly retained full-library normalization.</p><img src="figures/03B_candidate_real_time.png"><img src="figures/03B_program_real_time.png">',
'<h2>4. Intersect evidence and decide the next experimental stage</h2>',selection,'<img src="figures/04_candidate_screening_flow.png">',counts.to_html(index=False),
'<p>The rule yields AHNAK in Naive and MDFIC, TTC39C, ST6GALNAC3 in ActCD4. All are exploratory nominal candidates, not FDR-confirmed discovery genes. AHNAK, TTC39C and ST6GALNAC3 retain nominal trajectory evidence in all four runs; MDFIC does so in two runs, with P=0.0625 in the other two. No original candidate was guaranteed retention; the prior-target disposition table records every failed step.</p>',cards.to_html(index=False,float_format=lambda x:f'{x:.5g}'),
'<h3>Independent and human evidence changes the follow-up priority</h3><p>GSE183682 supports positive AHNAK differences within Naive at D7/D56 and positive ST6GALNAC3 differences within ActCD4 at both times. In contrast, ActCD4 TTC39C and MDFIC have negative differences at both external times. Each external condition/time has one library: these are direction checks only. In the reviewed human BALF Core comparison (2 exposed controls, 3 silicosis donors), no candidate has nominal P&lt;0.05; AHNAK has an opposite overall direction. Human subtype results and pre-recorded flag exclusions remain in the supporting tables. No human-mouse state equivalence or clinical validation is claimed.</p><p>None of the four belongs to the installed broad GO autophagy panel, and no candidate has a nominal candidate-autophagy association after removing candidate overlap and residualizing depth and pseudotime within libraries. Legacy regulatory/communication subsets are linked only as historical supplemental evidence; new genes absent from a restricted old panel are unassessed, not negative. GSE296682 remains patient localization only; unextracted candidate features are explicitly unknown.</p>',
'<h3>Experiment decision table</h3><table><tr><th>Candidate/state</th><th>Immediate purpose</th><th>Decision before mechanism work</th></tr><tr><td>AHNAK / Naive</td><td>Prioritize independent state-specific expression and localization checks; both external mouse time directions agree.</td><td>Resolve the opposite human Core direction; do not call this a human-validated target.</td></tr><tr><td>ST6GALNAC3 / ActCD4</td><td>Prioritize expression/detection checks in independently defined ActCD4; both external mouse time directions agree.</td><td>Confirm adequate RNA/protein detection and state specificity before perturbation.</td></tr><tr><td>TTC39C / ActCD4</td><td>Keep as a nominal intersection candidate with robust root sensitivity.</td><td>Resolve opposite external mouse directions before prioritizing mechanism work.</td></tr><tr><td>MDFIC / ActCD4</td><td>Keep as a secondary exploratory candidate.</td><td>Resolve root/graph sensitivity and opposite external mouse directions.</td></tr></table>',
'<p>First validate state identity, candidate expression and disease direction in independent biological samples. If reproduced, test candidate perturbation with matched controls and rescue, while keeping pyroptosis, apoptosis, mitochondrial function and autophagy as separate outcomes. Do not assume that candidate upregulation drives cell death. A protein/functional panel should distinguish gasdermin cleavage/membrane injury from apoptosis and measure autophagic flux rather than a single static marker. These are proposed validation steps, not results of this analysis. Assay rationale: <a href="https://www.nature.com/articles/nature15514">gasdermin-D cleavage study</a>; <a href="https://pubmed.ncbi.nlm.nih.gov/33634751/">autophagy assay guidelines</a>.</p>',
'<h2>Methods, cohort boundaries and reproducibility</h2>',methods,
'<p>GSE314943 is the silica-mouse discovery cohort. GSE274372 is real-time support. GSE183682 contributes direction only. GSE174725 provides limited reviewed human support; GSE296682 provides patient localization only. GSE288887/CWP is excluded from main discovery and target selection; PRJNA1289989 was not downloaded. No cross-species expression matrix was pooled. The old CWP report cannot be relabelled a pure SiO2-selected human cohort.</p><p>Code, full tables, sample memberships, panel definitions and run exit records are included. Complete gene-trajectory P values were independently checked using SciPy sign-flip enumeration; BH was checked within every state/run and time contrast. Inferred DPT ordering is conditional on the shared graph; biological direction is not causal ancestry. Method references: <a href="https://scanpy.readthedocs.io/en/1.9.x/generated/scanpy.tl.dpt.html">Scanpy DPT</a>; <a href="https://amigo.geneontology.org/amigo/term/GO%3A0006914">GO:0006914</a>.</p></body></html>']
(O/'project_report.html').write_text('\n'.join(parts),encoding='utf-8')
wr(pd.DataFrame([dict(stage='01',directory='01_CD4_death',question='Which CD4 state changes in pyroptosis/apoptosis?',status='complete'),dict(stage='02',directory='02_mito_autophagy',question='How do mitochondrial/autophagy programs change?',status='complete'),dict(stage='03A',directory='03A_state_pseudotime',question='Which within-state gene/program trends repeat across samples?',status='complete'),dict(stage='03B',directory='03B_real_time',question='Which measured-time effects exist with all12 libraries?',status='complete'),dict(stage='04',directory='04_candidate_intersection',question='Which genes satisfy the fixed intersection and what should be validated?',status='complete'),dict(stage='05',directory='05_supporting_evidence',question='Which supportive, contradictory and unmeasured evidence accompanies each candidate?',status='complete_with_explicit_coverage_limits')]),'project_stage_index.tsv')
qa={'status':'PASS','independent_trajectory_tests_checked':sum(v['valid_tests'] for v in validation),'trajectory_families':len(validation),'time_BH_families':4,'no_sample_cell_count_gate':True,'main_samples_each_scope':6,'time_samples':12,'primary_candidates':cards.gene.to_list(),'final_customer_translation_generated':False}
(O/'provenance/report_numeric_validation.json').write_text(json.dumps(qa,indent=2),encoding='utf-8');print(json.dumps(qa));print(cards[['gene','mouse_state','experimental_stage']].to_string(index=False))