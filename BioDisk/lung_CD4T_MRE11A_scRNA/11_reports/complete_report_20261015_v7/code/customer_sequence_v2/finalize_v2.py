from pathlib import Path
import sys,json,re,hashlib,zipfile,csv
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
O=Path(sys.argv[1]);P=O.parents[1]
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
def read(n):return pd.read_csv(O/n,sep='\t')
def save(fig,n):
 for ext in ['png','svg']:fig.savefig(O/'figures'/(n+'.'+ext),dpi=170,bbox_inches='tight',facecolor='white')
 plt.close(fig)
legend=[Line2D([0],[0],color='#2c6f9f',marker='o',label='Saline: each line is one library'),Line2D([0],[0],color='#c84d3e',marker='o',label='Silica: each line is one library')]
a=read('03A_state_pseudotime/program_sample_quantile_curves.tsv');tests=read('03A_state_pseudotime/program_trajectory_tests.tsv')
fig,axs=plt.subplots(2,4,figsize=(16,7.5))
for i,state in enumerate(['Naive','ActCD4']):
 for j,mod in enumerate(['PyroExec','ApoExec','Mitophagy','Autophagy_GO']):
  ax=axs[i,j];v=a[(a.state==state)&(a.module==mod)]
  for sample,b in v.groupby('sample'):ax.plot(b.mean_pt,b.score,'o-',lw=.8,ms=3,alpha=.7,color='#2c6f9f' if sample.startswith('Saline') else '#c84d3e')
  r=tests[(tests.state==state)&(tests.run=='K15_R1')&(tests.feature==mod)].iloc[0]
  ax.set(title=f'{state}: {mod}\nslope {r.mean_slope:+.3g}; P {r.P:.3g}; BH {r.BH:.3g}',xlabel='Within-state rescaled DPT',ylabel='Program score' if j==0 else '')
fig.suptitle('Within-state program trajectories: P uses six depth-adjusted library slopes',fontsize=15);fig.legend(handles=legend,loc='lower center',ncol=2,frameon=False);fig.tight_layout(rect=[0,.06,1,.94]);save(fig,'03A_program_trajectories')
a=read('03A_state_pseudotime/candidate_sample_curves.tsv');cards=read('04_candidate_intersection/candidate_evidence_cards.tsv')
fig,axs=plt.subplots(1,4,figsize=(16,4.8))
for ax,r in zip(axs,cards.itertuples()):
 v=a[(a.state==r.mouse_state)&(a.gene==r.gene)]
 for sample,b in v.groupby('sample'):ax.plot(b.mean_pt,b.mean_logCP10k,'o-',lw=1,ms=4,alpha=.7,color='#2c6f9f' if sample.startswith('Saline') else '#c84d3e')
 ax.set(title=f'{r.mouse_state}: {r.gene}',xlabel='Within-state DPT',ylabel='Mean log1p(CP10k)')
fig.suptitle('Candidate trajectories: sample curves are descriptive; P uses library slopes',fontsize=14);fig.legend(handles=legend,loc='lower center',ncol=2,frameon=False);fig.tight_layout(rect=[0,.10,1,.90]);save(fig,'03A_candidate_gene_trajectories')
t=read('03B_real_time/sample_TMM_logCPM.tsv');t['canonical']=t.symbol.str.upper();m=read('03B_real_time/sample_inclusion.tsv');days=[0,7,28,56]
fig,axs=plt.subplots(1,4,figsize=(16,4.5))
for ax,r in zip(axs,cards.itertuples()):
 row=t[t.canonical==r.gene].iloc[0];means=[]
 for day in days:
  v=row[m[m.day==day]['sample']].astype(float).to_numpy();ax.scatter(day+np.linspace(-.65,.65,len(v)),v,c='#286c93',s=30);means.append(v.mean())
 ax.plot(days,means,c='#c44f3d');ax.set(title=r.gene,xticks=days,xticklabels=['D0','D7','D28','D56'],xlabel='Measured day (true spacing)',ylabel='TMM log2CPM')
fig.suptitle('GSE274372 real time: all 12 TLike libraries; points=samples, line=day means',fontsize=14);fig.tight_layout(rect=[0,0,1,.90]);save(fig,'03B_candidate_real_time')
s=read('03B_real_time/all12_existing_program_scores.tsv');fig,axs=plt.subplots(2,4,figsize=(15,7))
for ax,mod in zip(axs.flat,['PyroExec','Inflam','ApoExec','OXPHOS','MitoISR','MitoROS','Mitophagy']):
 means=[]
 for day in days:
  v=s[s.day==day]['score_'+mod].to_numpy();ax.scatter(day+np.linspace(-.65,.65,len(v)),v,c='#286c93');means.append(v.mean())
 ax.plot(days,means,c='#c44f3d');ax.set(title=mod,xticks=days,xticklabels=['D0','D7','D28','D56'],xlabel='Measured day',ylabel='Sample program score')
axs.flat[-1].axis('off');axs.flat[-1].text(.05,.6,'All 12 samples retained.\nPoints: one library each.\nLines: equal-library day means.\nDescriptive curves only.\nNo cross-cohort absolute comparison.',transform=axs.flat[-1].transAxes)
fig.tight_layout();save(fig,'03B_program_real_time')
# Complete gene universes and mapping exclusions.
pt=read('03A_state_pseudotime/genomewide_trajectory_tests.tsv.gz');universes=[];miss=[]
for state in ['Naive','ActCD4']:
 de=read(f'04_candidate_intersection/{state}_all_gene_decisions.tsv.gz')
 universes.append(dict(state=state,disease_tested=len(de),trajectory_tested=int(de.PT_P.notna().sum()),time_tested_in_state_universe=int(de.time_omnibus_P.notna().sum()),tested_in_all_three=int((de.PT_P.notna()&de.time_omnibus_P.notna()).sum())))
 for r in de[de.PT_P.isna()].itertuples():miss.append(dict(state=state,gene=r.canonical,reason='not mapped to canonical counted expression matrix; no trajectory P assigned'))
pd.DataFrame(universes).to_csv(O/'04_candidate_intersection/gene_universe_summary.tsv',sep='\t',index=False,encoding='utf-8');pd.DataFrame(miss,columns=['state','gene','reason']).to_csv(O/'04_candidate_intersection/trajectory_mapping_exclusions.tsv',sep='\t',index=False,encoding='utf-8')
# Complete time source provenance; verify prior input manifests before finalization.
records=[]
for n in ['provenance/input_manifest_phase123A.tsv','05_supporting_evidence/input_manifest_support.tsv']:
 for r in read(n).to_dict('records'):
  f=Path(r['path']);assert hashlib.sha256(f.read_bytes()).hexdigest()==r['sha256'];records.append(r)
base=P/'02_analysis/GSE274372/stage10_20260908_141459_a8IPaq/tables'
for name in ['sample_summary.tsv','TLike_counts.tsv.gz','features.tsv.gz']:
 f=base/name;records.append(dict(path=str(f),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest(),role='real-time model source'))
pd.DataFrame(records).drop_duplicates('path').to_csv(O/'provenance/complete_input_manifest.tsv',sep='\t',index=False,encoding='utf-8')
s=(O/'project_report.html').read_text(encoding='utf-8')
s=s.replace('<h2>Methods, cohort boundaries and reproducibility</h2>','<h2>Methods, cohort boundaries and reproducibility</h2>'+pd.DataFrame(universes).to_html(index=False)+'<p>Trajectory mapping exclusions are listed separately; unmapped features receive no trajectory P and are never treated as observed negative expression.</p>')
(O/'project_report.html').write_text(s,encoding='utf-8')
for image in re.findall(r'<img src="([^"]+)"',s):assert (O/image).exists()
readme='''CD4 CUSTOMER-SEQUENCE REVISION V2 — INTERNAL WORKING PROJECT

Start: project_report.html. project_stage_index.tsv maps stages to folders.

01_CD4_death: within-state pyroptosis/apoptosis changes, all nonempty libraries.
02_mito_autophagy: mitochondrial/mitophagy and independent broad GO-autophagy panel.
03A_state_pseudotime: six-library slopes within Naive/ActCD4; four fixed DPT runs.
03B_real_time: all12 GSE274372 TLike libraries; real day, not pseudotime.
04_candidate_intersection: prospective selection, counts/members, full decisions,
  old-candidate disposition, candidate evidence cards and experiment staging.
05_supporting_evidence: independent mouse directions, reviewed human donor data,
  patient-localization coverage and explicitly limited legacy network links.

No minimum 20-cell or 10-cell sample gate. Nonempty measurable units retained.
Expression/QC/identity rules remain valid and distinct from sample cell counts.
Current candidate list: AHNAK (Naive); MDFIC, TTC39C, ST6GALNAC3 (ActCD4).
These are nominal exploratory intersections, not FDR-confirmed mechanistic targets.
Prioritize expression/state checks for AHNAK and ST6GALNAC3. Resolve external
mouse discordance for TTC39C/MDFIC and root sensitivity for MDFIC.
No candidate establishes an autophagy mechanism; human support remains limited.

Reproduction order (server paths remain external dependencies):
Rscript --vanilla code/annotation_v2.R <this_directory>
python code/compute_v2.py <this_directory>
Rscript --vanilla code/time_v2.R <this_directory>
python code/intersect_v2.py <this_directory>
python code/support_v2.py <this_directory>
python code/report_v2.py <this_directory>
python code/finalize_v2.py <this_directory>

R: /mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript
Python: /mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
Use single-thread BLAS/OMP. No upstream launchers are executed or overwritten.
Large upstream matrices remain on the server; this is not a standalone raw-data image.

Validation: all44,456 gene trajectory P values independently checked against SciPy;
8 trajectory and4 time-model BH families checked; seven PNG figures visually reviewed.
Original inputs hash-checked unchanged. All output paths remain inside this revision.
Historical report/v1 are preserved. GSE288887 excluded from discovery/selection.
No PRJNA1289989 download. No Chinese customer version or external delivery generated.
'''
(O/'README.txt').write_text(readme,encoding='utf-8')
print('FINALIZE_PASS');print(pd.DataFrame(universes).to_string(index=False));print('MAPPING_EXCLUSIONS',miss)