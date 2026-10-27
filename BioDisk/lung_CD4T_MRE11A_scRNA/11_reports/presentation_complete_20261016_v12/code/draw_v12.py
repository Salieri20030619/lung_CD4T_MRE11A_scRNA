from pathlib import Path
import json,hashlib,shutil,re
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, FormatStrFormatter
from matplotlib.colors import Normalize
from scipy.stats import spearmanr
R=Path('C:/codex_test');O=R/'output/silicosis_v12';T=R/'tmp/silicosis_v12';V4=R/'tmp/silicosis_v9/v4';I=T/'inputs_R2';OLD=R/'output/silicosis_v11'
for folder in ['assets','figures','plot_data','source_tables','supplement','provenance','code']:(O/folder).mkdir(parents=True,exist_ok=True)
for folder in ['assets','figures','plot_data','source_tables','supplement']:shutil.copytree(OLD/folder,O/folder,dirs_exist_ok=True)
shutil.copytree(I,O/'source_tables/presentation_inputs',dirs_exist_ok=True)
def rd(f):return pd.read_csv(V4/f,sep='\t')
def dump(d,n):d.to_csv(O/'supplement'/n,sep='\t',index=False,encoding='utf-8')
plt.rcParams.update({'font.family':'Microsoft YaHei','font.size':11,'axes.unicode_minus':False,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','savefig.facecolor':'white'})
BLUE='#2677a3';RED='#c55b43';TEAL='#167e78';INK='#17384b';GRAY='#6c7b86'
def save(fig,n):
 for ext in ['png','svg','pdf']:fig.savefig(O/'figures'/f'{n}.{ext}',dpi=220,bbox_inches='tight')
 plt.close(fig)
score=rd('01_phenotype/sample_program_scores.tsv');stats=rd('01_phenotype/01_CD4_death_group_tests.tsv');a=score[score['filter']=='All'];primary=stats[stats['filter']=='All']
defs=pd.read_csv(I/'module_definition.tsv',sep='\t');coverage=pd.read_csv(I/'module_coverage.tsv',sep='\t');controls=pd.read_csv(I/'score_controls.tsv.gz',sep='\t')
dump(score,'S18_sample_program_scores.tsv');dump(stats,'S18_death_program_statistics.tsv');dump(defs[defs.module.isin(['PyroExec','ApoExec','Inflam'])],'S18_death_module_genes.tsv');dump(coverage[(coverage.dataset=='GSE314943')&coverage.module.isin(['PyroExec','ApoExec','Inflam'])],'S18_death_module_coverage.tsv');dump(controls[(controls.dataset=='GSE314943')&controls.module.isin(['PyroExec','ApoExec','Inflam'])],'S18_death_matched_controls.tsv')
for module,num,cn,name in [('PyroExec',1,'焦亡相关转录评分','V12_01_pyroptosis_scores'),('ApoExec',2,'凋亡执行相关转录评分','V12_02_apoptosis_scores'),('Inflam',None,'炎症小体相关转录评分','S18_inflammasome_scores')]:
 fig,axes=plt.subplots(1,4,figsize=(15,4.8),sharey=True)
 fig.suptitle((f'图{num}  ' if num else '补充图  ')+cn,x=.055,ha='left',y=.98,fontsize=20,fontweight='bold',color=INK)
 fig.text(.055,.875,f'{module}  |  每点为一只小鼠；对照3只、硅尘3只；横线为组均值',fontsize=11,color=GRAY)
 for ax,scope in zip(axes,['Core','Naive','ActCD4','Treg']):
  s=a[a.scope==scope];r=primary[(primary.scope==scope)&(primary.module==module)].iloc[0]
  for x,group,color,prefix in [(0,'Saline',BLUE,'C'),(1,'Silica',RED,'S')]:
   b=s[s.group==group].sort_values('sample');xs=x+np.array([-.12,0,.12]);ax.scatter(xs,b[module],color=color,s=62,zorder=3,edgecolor='white',linewidth=.6)
   for k,(xx,(_,row)) in enumerate(zip(xs,b.iterrows())):ax.annotate(prefix+row['sample'].split('_')[1],(xx,row[module]),xytext=[(-14,9),(0,-15),(9,7)][k],textcoords='offset points',fontsize=9,color=color)
   ax.hlines(b[module].mean(),x-.23,x+.23,color=color,lw=2)
  ax.set_title(f'{scope}\nΔ={r.difference:+.5f}\nP={r.P:.2g}；BH={r.BH:.2g}',fontsize=12,pad=12);ax.set_xticks([0,1],['对照','硅尘']);ax.set_xlim(-.4,1.45);ax.axhline(0,color='#d6dfe5',lw=.8,zorder=0)
 axes[0].set_ylabel('样本平均评分（相对匹配对照基因）',fontsize=11)
 fig.subplots_adjust(left=.07,right=.99,top=.64,bottom=.22,wspace=.18)
 genes='、'.join(defs.loc[defs.module==module,'gene']);fig.text(.055,.085,'基因集：'+genes,fontsize=11);fig.text(.055,.025,'Δ为硅尘减对照；负评分表示低于匹配对照基因的平均表达，不代表负的细胞死亡率。',fontsize=10,color=GRAY)
 save(fig,name)
 # Independently recover the displayed effect from frozen sample values.
 for scope,s in a.groupby('scope'):
  r=primary[(primary.scope==scope)&(primary.module==module)].iloc[0]
  assert np.isclose(s.loc[s.group=='Silica',module].mean()-s.loc[s.group=='Saline',module].mean(),r.difference,atol=1e-13)

# Replace the abstract coefficient-joining line with actual sample coordinates.
act=rd('03_receiver/TF_all_regulon_sample_activity.tsv.gz');act=act[(act.scope=='Core')&(act.seed==17)]
observed=act.pivot(index='sample',columns='regulon',values='mean_auc').join(a[a.scope=='Core'].set_index('sample')[['ApoExec','MitoROS']])
fixed=rd('06_support/Core_fixed_composition_TF_inputs.tsv').set_index('sample').join(rd('06_support/Core_fixed_composition_program_inputs.tsv').set_index('sample'))
pooled=rd('03_receiver/all_TF_program_associations.tsv.gz');pooled=pooled[(pooled.scope=='Core')&(pooled.seed==17)&(pooled.model=='Pooled')]
fstats=rd('06_support/Core_fixed_composition_focal.tsv');fstats=fstats[fstats.model=='Pooled'];weights=rd('06_support/Core_fixed_composition_weights.tsv')
fig,axes=plt.subplots(2,2,figsize=(16,6.7));fig.suptitle('图9  亚群比例统一后，同一样本的关联如何变化',x=.06,ha='left',fontsize=20,fontweight='bold',color=INK,y=.99)
fig.text(.06,.93,'两列均使用同一批6只小鼠。左：各样本实际亚群比例；右：统一Naive、ActCD4、Treg比例。',fontsize=11,color=GRAY)
coords=[];statrows=[]
for i,(tf,module,yname) in enumerate([('Ikzf2(+)','ApoExec','凋亡执行评分'),('Stat1(+)','MitoROS','线粒体抗氧化应答评分')]):
 xall=np.r_[observed[tf],fixed[tf]];yall=np.r_[observed[module],fixed[module]];dx=np.ptp(xall)*.14;dy=np.ptp(yall)*.18
 for j,(data,lab) in enumerate([(observed,'实际亚群比例'),(fixed,'统一亚群比例')]):
  ax=axes[i,j];row=pooled[(pooled.x==tf)&(pooled.y==module)].iloc[0] if j==0 else fstats[(fstats.regulon==tf)&(fstats.module==module)].iloc[0]
  assert np.isclose(spearmanr(data[tf],data[module]).statistic,row.rho,atol=1e-12)
  for sample,z in data.iterrows():
   color=BLUE if sample.startswith('Saline') else RED;tag=('C' if sample.startswith('Saline') else 'S')+sample.split('_')[1]
   offset=(-18,-16) if tag=='C3' else (5,-16) if tag=='S2' else (5,7)
   ax.scatter(z[tf],z[module],s=68,color=color,edgecolor='white',zorder=3);ax.annotate(tag,(z[tf],z[module]),xytext=offset,textcoords='offset points',fontsize=10,color=color)
   coords.append(dict(sample=sample,group='Saline' if sample.startswith('Saline') else 'Silica',composition=lab,tf=tf,module=module,tf_auc=z[tf],program_score=z[module]))
  ax.set_xlim(xall.min()-dx,xall.max()+dx);ax.set_ylim(yall.min()-dy,yall.max()+dy);ax.set_xlabel(tf.replace('(+)','').upper()+' 调控子活性（AUC）',fontsize=12);ax.set_ylabel(yname,fontsize=12)
  ax.xaxis.set_major_locator(MaxNLocator(5));ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'));ax.yaxis.set_major_locator(MaxNLocator(4))
  ax.set_title(f'{lab}  |  ρ={row.rho:+.3f}，P={row.P:.4g}，BH={row.BH:.4g}',fontsize=12,pad=10)
  statrows.append(dict(tf=tf,module=module,composition=lab,rho=row.rho,P=row.P,BH=row.BH))
fig.subplots_adjust(left=.08,right=.98,top=.83,bottom=.2,hspace=.64,wspace=.3)
fig.text(.06,.09,'C1-C3：对照；S1-S3：硅尘。统一权重：Naive 46.66%，ActCD4 31.00%，Treg 22.34%。',fontsize=11)
fig.text(.06,.035,'IKZF2-凋亡关联由正转负；STAT1-MitoROS仍为正。两列均为样本间Spearman关联，不是组别调整或因果效应。',fontsize=10,color=GRAY)
save(fig,'V12_09_composition_scatter');dump(pd.DataFrame(coords),'S18_composition_scatter_coordinates.tsv');dump(pd.DataFrame(statrows),'S18_composition_association_statistics.tsv');dump(weights,'S18_fixed_composition_weights.tsv')

# Explicit NicheNet activities with target-set sizes, beyond ordinal ranks.
rank=rd('02_sender/source_specific_activity_ranks.tsv');activity=rd('02_sender/ligand_activity_all_states.tsv');elig=rd('02_sender/ligand_activity_evaluability.tsv');members=rd('02_sender/receiver_target_background_members.tsv')
dump(rank,'S19_NicheNet_all_source_rankings.tsv');dump(activity,'S19_NicheNet_all_activities.tsv');dump(elig,'S19_NicheNet_evaluability.tsv');dump(members,'S19_NicheNet_target_background_members.tsv')
fig,axes=plt.subplots(1,3,figsize=(15,5.6));fig.suptitle('图14  NicheNet预测：哪些配体更能解释接收端上调基因',x=.055,ha='left',fontsize=19,fontweight='bold',color=INK,y=.99)
fig.text(.055,.91,'AMlike来源表达支持的配体；各亚群按corrected AUPR选前8名。颜色仅区分当前重点配体。',fontsize=11,color=GRAY)
plotr=[]
for ax,scope in zip(axes,['Naive','ActCD4','Treg']):
 s=rank[(rank.sender=='AMlike')&(rank.scope==scope)&(rank.direction=='Up')].sort_values('source_specific_rank').head(8).copy();plotr.append(s);s=s.iloc[::-1];r=elig[(elig.scope==scope)&(elig.direction=='Up')].iloc[0]
 ax.barh(s.ligand,s.aupr_corrected,color=[RED if x in ['CXCL16','TNF','IL15','IL18','C3'] else BLUE for x in s.ligand],height=.65)
 for y,(_,row) in enumerate(s.iterrows()):ax.text(row.aupr_corrected+max(s.aupr_corrected)*.025,y,f'{row.aupr_corrected:.4f}',va='center',fontsize=10)
 ax.set_xlim(0,max(s.aupr_corrected)*1.38);ax.set_title(f'{scope}\n上调目标 {r.targets}；背景 {r.background}',fontsize=12,pad=15);ax.set_xlabel('corrected AUPR（预测指标）');ax.tick_params(axis='y',length=0)
fig.subplots_adjust(left=.08,right=.99,top=.72,bottom=.21,wspace=.5)
fig.text(.055,.09,'目标集：同亚群差异表达原始P<0.05且log2FC>0；背景：受检、可映射至先验且有表达输入的基因交集。',fontsize=10)
fig.text(.055,.035,'corrected AUPR不是P值；各面板目标集不同，不把跨面板数值差异解释为细胞间信号强弱。完整上/下调结果见S19。',fontsize=10,color=GRAY)
save(fig,'V12_14_nichenet_activity');dump(pd.concat(plotr),'S19_NicheNet_activity_figure_values.tsv')

# Overlay unchanged prior weights and measured ActCD4 DE, without refitting death-target NicheNet.
prior=pd.read_csv(I/'nichenet_death_candidate_prior_weights.tsv.gz',sep='\t');assert len(prior)==33462
focal=['CXCL16','TNF','IL15','IL18','C3'];tops=rank[(rank.sender=='AMlike')&(rank.scope=='ActCD4')&(rank.direction=='Up')].sort_values('source_specific_rank').head(5).ligand.tolist();ligands=list(dict.fromkeys(tops+focal))
pyro=['PYCARD','CASP1','GSDMD','IL1B','IL18'];apo=['BAX','BAK1','APAF1','CASP9','CASP3','CASP7','CYCS'];genes=pyro+apo
matrix=prior.pivot(index='ligand',columns='target',values='weight').reindex(index=ligands,columns=genes)
de=rd('03_receiver/ActCD4_DE_all_genes.tsv.gz').set_index('canonical');bg=members[(members.scope=='ActCD4')&(members.direction=='Up')].set_index('gene')
status=[2 if g not in bg.index else 1 if bool(bg.loc[g,'in_target']) else 0 for g in genes]
assert [g for g,s in zip(genes,status) if s==1]==['GSDMD']
genevalues=[]
for g,s in zip(genes,status):
 genevalues.append(dict(gene=g,module='PyroExec' if g in pyro else 'ApoExec',target_membership=['Background_only','Up_target','Outside_DE_background'][s],logFC=float(de.loc[g,'logFC']) if g in de.index else np.nan,P=float(de.loc[g,'PValue']) if g in de.index else np.nan,BH=float(de.loc[g,'FDR']) if g in de.index else np.nan))
gv=pd.DataFrame(genevalues);sub=prior[prior.ligand.isin(ligands)&prior.target.isin(genes)].copy();dump(sub,'S19_NicheNet_death_weight_figure_values.tsv');dump(gv,'S19_ActCD4_death_genes_DE_and_target_membership.tsv');shutil.copy2(I/'nichenet_death_candidate_prior_weights.tsv.gz',O/'supplement/S19_NicheNet_all_death_candidate_prior_weights.tsv.gz')
fig=plt.figure(figsize=(16,6.5));fig.suptitle('图15  NicheNet配体-死亡程序基因先验与本队列证据',x=.055,ha='left',fontsize=19,fontweight='bold',color=INK,y=.99)
fig.text(.055,.925,'配体范围：ActCD4 / AMlike预测前5名，加CXCL16、TNF、IL15、IL18、C3；不按死亡基因权重重新挑选配体。',fontsize=10,color=GRAY)
gs=fig.add_gridspec(3,1,height_ratios=[.38,5,1],left=.12,right=.86,top=.81,bottom=.25,hspace=.18)
ax=fig.add_subplot(gs[0]);ax.set_xlim(-.5,11.5);ax.set_ylim(-.5,.5);ax.axis('off')
for j,(g,s) in enumerate(zip(genes,status)):ax.add_patch(plt.Rectangle((j-.48,-.4),.96,.8,color=['#c8d8e5',RED,'#e4e6e8'][s]));ax.text(j,0,['背景','上调目标','未入背景'][s],ha='center',va='center',fontsize=8)
ax.text(2,1.0,'焦亡相关模块',ha='center',color=RED,fontweight='bold');ax.text(8.1,1.0,'凋亡执行模块',ha='center',color=BLUE,fontweight='bold')
ax=fig.add_subplot(gs[1]);im=ax.imshow(matrix,aspect='auto',cmap='YlGnBu',vmin=0,vmax=float(matrix.max().max()));ax.set_yticks(range(len(ligands)),ligands);ax.set_xticks([]);ax.tick_params(length=0);ax.axvline(4.5,color='white',lw=3)
for i in range(len(ligands)):
 for j in range(len(genes)):ax.text(j,i,f'{matrix.iloc[i,j]:.3f}',ha='center',va='center',fontsize=10,color='white' if matrix.iloc[i,j]>matrix.max().max()*.6 else 'black')
cax=fig.add_axes([.89,.39,.013,.27]);fig.colorbar(im,cax=cax,label='原始先验权重（无符号）')
ax=fig.add_subplot(gs[2]);vals=gv.logFC.to_numpy()[None];lim=float(np.nanmax(abs(vals)));im=ax.imshow(np.ma.masked_invalid(vals),aspect='auto',cmap='RdBu_r',vmin=-lim,vmax=lim);ax.set_facecolor('#e4e6e8');ax.set_yticks([0],['ActCD4\n疾病log2FC']);ax.set_xticks(range(12),genes,rotation=45,ha='right');ax.tick_params(length=0);ax.axvline(4.5,color='white',lw=3)
for j,g in enumerate(genes):
 if np.isfinite(gv.iloc[j].logFC):ax.text(j,0,f'{gv.iloc[j].logFC:+.2f}'+('*' if gv.iloc[j].P<.05 else ''),ha='center',va='center',fontsize=9,color='white' if abs(gv.iloc[j].logFC)>lim*.7 else 'black')
cax=fig.add_axes([.89,.25,.013,.09]);fig.colorbar(im,cax=cax,label='log2FC')
fig.text(.055,.09,'上方条带：ActCD4上调目标集仅纳入GSDMD，其余为背景基因或未进入差异背景。下方*表示原始P<0.05。',fontsize=10)
fig.text(.055,.035,'先验权重不等于本队列效应或显著性，也无激活/抑制方向。本图是既有预测的模块映射，未运行死亡程序定向预测。',fontsize=10,color=GRAY)
save(fig,'V12_15_nichenet_death_targets')

# Re-number the six previously added data figures using their original builder.
code=(R/'tmp/silicosis_v9/build_v9.py').read_text(encoding='utf-8').split('# Preserve original evidence')[0]
code=code.replace("OUT=ROOT/'output/silicosis_v9'","OUT=ROOT/'output/silicosis_v12'")
for old,new,txt in [('01','12','Macrophage-to-CD4'),('02','13','Predicted ligand'),('03','07','CD4 regulon'),('04','08','Regulons connect'),('05','10','STAT1 prior'),('06','16','The ligand-to-STAT1')]:code=code.replace(old+'  '+txt,new+'  '+txt)
exec(compile(code,'renumber_frozen_figures','exec'),{})
(O/'provenance/presentation_update_summary.json').write_text(json.dumps({'new_inferential_tests':0,'death_module_sizes':{'PyroExec':5,'ApoExec':7,'Inflam':5},'matched_control_sizes':{'PyroExec':100,'ApoExec':140,'Inflam':100},'ActCD4_death_genes_in_Up_target':['GSDMD'],'prior_export_rows':len(prior),'scatter_coordinate_rows':len(coords),'nichenet_ligands_shown':ligands},indent=2,ensure_ascii=False),encoding='utf-8')
shutil.copy2(__file__,O/'code/draw_v12.py');shutil.copy2(T/'export_inputs.R',O/'code/export_inputs.R')
print('PLOTS_COMPLETE',len(coords),'scatter rows;',len(prior),'prior rows')
