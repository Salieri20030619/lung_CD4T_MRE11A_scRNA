from pathlib import Path
import json,numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
D=Path(__file__).parent/'revision_delivery';data=json.loads((D/'report_evidence_aggregate.json').read_text(encoding='utf-8'))
df=lambda k:pd.DataFrame(data[k])
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
def save(fig,name):fig.savefig(D/name,dpi=200,bbox_inches='tight');plt.close(fig)
# Frozen cohort flow, summary values from source inventories.
fig,ax=plt.subplots(figsize=(11,3.1));ax.axis('off')
boxes=[(.01,.55,.27,.36,'Human lung RNA\n19 donors: 10 Ctrl / 9 CWP\nPrimary: 6,523 CD4 cells'),(.36,.55,.28,.36,'SCENIC / AUCell\nCore >=20: 17 donors, 5,684 cells\nBG1: 15 donors, 5,298 cells'),(.72,.55,.27,.36,'Subtype / paired models\nCCR7 9 | GZMA 12 | IL7R 16\nTreg 10 | CD4state 7*')]
for x,y,w,h,s in boxes:
 ax.add_patch(plt.Rectangle((x,y),w,h,facecolor='#e9f3f5',edgecolor='#367c8c',lw=1.4));ax.text(x+w/2,y+h/2,s,ha='center',va='center',fontsize=10)
for a,b in [(.29,.35),(.65,.71)]:ax.annotate('',xy=(b,.73),xytext=(a,.73),arrowprops={'arrowstyle':'->','color':'#367c8c'})
ax.text(.5,.34,'Bootstrap: RNA19 + matched RNA17/AUC17 | 5,000 donor resamples per scheme',ha='center',fontsize=11,fontweight='bold')
ax.text(.5,.16,'MacAM: 14 eligible donors (9 Ctrl / 5 CWP) | MonoMac: 15 (7 / 8) | MacWide: union',ha='center')
ax.text(.5,.01,'* CD4state: 1 Ctrl / 6 CWP at >=20 cells; disease models insufficient. Counts depend on the declared scope.',ha='center',fontsize=9)
save(fig,'R01_flow.png')
# Full immune-panel directions in fixed Cells20 Source model, no re-fitting.
d=df('stage1_models');d=d[(d.kind=='program')&(d.roster=='Cells20')&(d.model=='Source')]
sc=['Primary','Core','Subtype:CD4-CCR7','Subtype:CD4-GZMA','Subtype:CD4-IL7R','Subtype:Treg','Subtype:CD4state'];fs=['score_TAct','score_TReg','score_Inhib','score_Naive'];a=d.pivot(index='scope',columns='metric',values='effect').reindex(index=sc,columns=fs);ps=d.pivot(index='scope',columns='metric',values='P').reindex(index=sc,columns=fs)
fig,ax=plt.subplots(figsize=(9,4.3));im=ax.imshow(a,cmap='RdBu_r',vmin=-.16,vmax=.16,aspect='auto');ax.set_xticks(range(4),['TAct','TReg','Inhib','Naive']);ax.set_yticks(range(7),[s.replace('Subtype:CD4-','').replace('Subtype:','') for s in sc])
for i in range(7):
 for j in range(4):
  v=a.iloc[i,j];pv=ps.iloc[i,j];ax.text(j,i,'NA' if not np.isfinite(v) else f'{v:+.3f}'+ (' *' if pv<.05 else ''),ha='center',va='center',fontsize=10,color='black')
ax.set_title('Human fixed immune programs | Source model | Cells20',loc='left',fontweight='bold');fig.colorbar(im,ax=ax,label='CWP - Ctrl adjusted score difference');fig.tight_layout();save(fig,'R02_CD4_programs.png')
# Bootstrap: original effect and exact saved percentile intervals, independent three association panels.
d=df('bootstrap_effect_summary');cols={'RNA19':'#667d91','RNA17':'#b17a4a','AUC17':'#136f7d','AUC29':'#568f98','AUC43':'#8cb6bc'}
for axis in ['IL18','TNF']:
 adj=['B'] if axis=='IL18' else ['B','B_Treg'];rows=[(m,b) for b in adj for m in cols];fig,axs=plt.subplots(1,3,figsize=(11.5,4.0 if axis=='IL18' else 6.3),sharey=True)
 for ax,edge,title in zip(axs,['ligand_TF','TF_phenotype','ligand_phenotype'],['Ligand vs TF','TF vs Effector score','Ligand vs Effector score']):
  for i,(m,b) in enumerate(rows):
   r=d[(d.scheme=='stratified')&(d.axis==axis)&(d.modality==m)&(d.adjustment==b)&(d.edge==edge)].iloc[0];ax.hlines(i,r.CI95_low,r.CI95_high,color=cols[m],lw=3);ax.plot(r.original_effect,i,'o',color=cols[m],ms=5);ax.plot(r.bootstrap_median,i,'|',color='black',ms=9)
  ax.axvline(0,color='#aaa',ls='--');ax.set_xlim(-1.02,1.02);ax.set_xlabel('Rank correlation');ax.set_title(title,fontsize=11);ax.grid(axis='x',alpha=.15)
 axs[0].set_yticks(range(len(rows)),[m+' / '+b for m,b in rows]);axs[0].invert_yaxis();fig.suptitle(axis+' | Stratified donor bootstrap | 5,000 draws',fontweight='bold');fig.tight_layout();save(fig,'R03_bootstrap_'+axis+'.png')
# Actual source-adjusted NicheNet corrected AUPR, all fixed candidates, negative scores retained.
d=df('human_NN_main');q=d[(d.model=='Source')&(d.test=='Nominal05')&(d.direction=='Up')];lig=['IL15','TNF','IL18','ICAM2','VEGFB'];sc=['Primary','Core','CD4-CCR7','CD4-GZMA','CD4-IL7R','Treg'];a=q.pivot(index='test_ligand',columns='receiver',values='aupr_corrected').reindex(index=lig,columns=sc);rk=q.pivot(index='test_ligand',columns='receiver',values='rank').reindex(index=lig,columns=sc)
fig,ax=plt.subplots(figsize=(10,4));im=ax.imshow(a,cmap='RdBu_r',vmin=-.016,vmax=.016,aspect='auto');ax.set_xticks(range(6),[s.replace('CD4-','') for s in sc]);ax.set_yticks(range(5),lig)
for i in range(5):
 for j in range(6):
  v=a.iloc[i,j];rank=rk.iloc[i,j];ax.text(j,i,'NA' if pd.isna(v) else f'{v:+.4f}\nrank {int(rank)}',ha='center',va='center',fontsize=8.5)
ax.set_title('Human NicheNet | All / Source / nominal Up targets',loc='left',fontweight='bold');fig.colorbar(im,ax=ax,label='Corrected AUPR');fig.tight_layout();save(fig,'R04_human_NicheNet.png')
# External transfer retains both concordant and discordant times.
d=df('S4_direction_replication');d=d[(d['filter']=='All')&d.scope.isin(['Core','ActCD4'])];rows=[(sc,tf) for sc in ['Core','ActCD4'] for tf in ['IKZF2','STAT1']];columns=[(day,seed) for day in [7,56] for seed in [17,29,43]];a=np.empty((4,6))
for i,(sc,tf) in enumerate(rows):
 for j,(day,seed) in enumerate(columns):a[i,j]=d[(d.scope==sc)&(d.feature==tf)&(d.day==day)&(d.seed==seed)].iloc[0].delta
fig,ax=plt.subplots(figsize=(10,3.4));im=ax.imshow(a,cmap='RdBu_r',vmin=-.009,vmax=.009,aspect='auto');ax.set_yticks(range(4),[sc+' / '+tf for sc,tf in rows]);ax.set_xticks(range(6),[f'D{day} / s{seed}' for day,seed in columns]);
for i in range(4):
 for j in range(6):ax.text(j,i,f'{a[i,j]:+.4f}',ha='center',va='center',fontsize=9)
ax.set_title('SCENIC4 transfer | GSE183682 | 1 library per condition/time',loc='left',fontweight='bold');fig.colorbar(im,ax=ax,label='Silica - NS mean AUC');fig.tight_layout();save(fig,'R05_SCENIC4.png')
print('FIGURES_COMPLETE')
