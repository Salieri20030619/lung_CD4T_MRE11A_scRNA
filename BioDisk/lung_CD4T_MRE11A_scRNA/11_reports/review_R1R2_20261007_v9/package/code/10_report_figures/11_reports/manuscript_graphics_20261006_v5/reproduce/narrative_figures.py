from pathlib import Path
import json,pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import fontManager
B=Path(__file__).parent;D=B/'narrative_delivery';O=D/'figures'
E=json.load(open(B/'revision_delivery/report_evidence_aggregate.json',encoding='utf8'));F=json.load(open(D/'added_evidence_aggregate.json',encoding='utf8'))
fontManager.addfont('C:/Windows/Fonts/msyh.ttc');plt.rcParams.update({'font.family':'Microsoft YaHei','font.size':9,'axes.titlesize':10,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':200,'axes.unicode_minus':False})
blue='#287B9B';red='#B85350';grey='#808C91'
def panel(ax,letter,title):ax.set_title(title,loc='left',pad=15);ax.text(-.13,1.06,letter,transform=ax.transAxes,fontweight='bold',fontsize=13)
fig,axs=plt.subplots(1,2,figsize=(8.4,3.6),gridspec_kw={'wspace':.55})
ax=axs[0];d=pd.DataFrame(E['human_TF_disease']).query("readout=='RNA' and scope=='CD4-GZMA' and tf=='IKZF2'")
for j,mo in enumerate(['Base','Source','SourceDep']):
 r=d[d.model==mo].iloc[0];ax.errorbar(j,r.effect,yerr=[[r.effect-r.CI_low],[r.CI_high-r.effect]],fmt='o',color=blue,capsize=3);ax.text(j,r.CI_high+.002,'P='+format(r.P,'.3g'),ha='center',fontsize=8)
ax.axhline(0,color=grey,lw=.7);ax.set_xticks(range(3),['基础','来源校正','来源与深度校正']);ax.set_ylabel('RNA group effect (95% CI)');panel(ax,'A','人类 GZMA 亚群 IKZF2 RNA')
ax=axs[1];d=pd.DataFrame(E['old_TF_group']).query("species=='Mouse' and filter=='All' and variant=='DropOverlap' and model=='Group'")
for tf,col in [('IKZF2',blue),('STAT1',red)]:
 for scope,sty in [('Core','-'),('ActCD4','--')]:
  rr=d[(d.tf==tf)&(d.scope==scope)].sort_values('seed');ax.plot(range(3),rr.coef,sty+'o',color=col,ms=4,label=tf+' / '+scope)
ax.axhline(0,color=grey,lw=.7);ax.set_xticks(range(3),['17','29','43']);ax.set_xlabel('Random seed');ax.set_ylabel('Δ regulon AUC');ax.legend(fontsize=7,frameon=False);panel(ax,'B','小鼠调控子疾病差异')
fig.savefig(O/'candidate_expression.png',bbox_inches='tight',facecolor='white');plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(8.4,4.1),gridspec_kw={'wspace':.55});d=pd.DataFrame(F['all_TF_program_associations']).query("species=='Human' and filter=='All' and variant=='DropOverlap' and scope=='Core' and model=='SourceDep'")
mods=['OXPHOS','MitoISR','MitoROS','Mitophagy','Inflam','PyroExec','ApoExec']
for ax,tf,letter in zip(axs,['IKZF2','STAT1'],['A','B']):
 arr=[];marks=[]
 for mod in mods:
  rr=[d[(d.tf==tf)&(d.seed==seed)&(d.feature=='score_'+mod)].iloc[0] for seed in [17,29,43]];arr.append([r.rho for r in rr]);marks.append(['**' if r.FDR<.05 else '*' if r.p<.05 else '' for r in rr])
 a=np.asarray(arr);im=ax.imshow(a,cmap='RdBu_r',vmin=-1,vmax=1,aspect='auto');ax.set_xticks(range(3),['17','29','43']);ax.set_yticks(range(7),mods);ax.set_xlabel('Random seed')
 for i in range(7):
  for j in range(3):ax.text(j,i,f'{a[i,j]:.2f}'+marks[i][j],ha='center',va='center',fontsize=8,color='white' if abs(a[i,j])>.72 else '#222222')
 plt.colorbar(im,ax=ax,shrink=.8,pad=.025,label='Partial Spearman ρ');panel(ax,letter,'核心 CD4：'+tf+'–功能程序')
fig.savefig(O/'candidate_function.png',bbox_inches='tight',facecolor='white');plt.close(fig)
print('Recomposed two figures from unchanged source values')
