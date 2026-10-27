from pathlib import Path
import json, pandas as pd, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties, fontManager
B=Path(__file__).parent;D=B/'manuscript_delivery';O=D/'figures';O.mkdir(exist_ok=True)
E=json.loads((B/'revision_delivery/report_evidence_aggregate.json').read_text(encoding='utf8'))
F=json.loads((D/'added_evidence_aggregate.json').read_text(encoding='utf8'))
de=lambda k:pd.DataFrame(E[k]);df=lambda k:pd.DataFrame(F[k])
fontManager.addfont('C:/Windows/Fonts/msyh.ttc')
plt.rcParams.update({'font.family':'Microsoft YaHei','font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':200,'axes.unicode_minus':False})
blue='#287B9B';red='#B85350';grey='#808C91'
labels=[]
def panel(ax,letter,title):
 ax.set_title(title,loc='left',pad=15);ax.text(-.13,1.06,letter,transform=ax.transAxes,fontweight='bold',fontsize=13);labels.extend([title,letter])
def save(fig,name):
 fig.savefig(O/(name+'.png'),bbox_inches='tight',facecolor='white');plt.close(fig)
def heat(ax,arr,x,y,vmax=1,mark=None,fmt='.2f',label=''):
 a=np.asarray(arr,dtype=float);im=ax.imshow(a,cmap='RdBu_r',vmin=-vmax,vmax=vmax,aspect='auto');ax.set_xticks(range(len(x)),x);ax.set_yticks(range(len(y)),y)
 for i in range(a.shape[0]):
  for j in range(a.shape[1]):
   if np.isfinite(a[i,j]):ax.text(j,i,format(a[i,j],fmt)+(mark[i][j] if mark is not None else ''),ha='center',va='center',fontsize=7,color='white' if abs(a[i,j])>.72*vmax else '#222222')
   else:ax.text(j,i,'NA',ha='center',va='center',fontsize=7,color='#666666')
 plt.colorbar(im,ax=ax,shrink=.7,pad=.025,label=label);labels.extend(list(x)+list(y)+[label]);return im
def sig(p,q):return '**' if pd.notna(q) and q<.05 else '*' if pd.notna(p) and p<.05 else ''

# Cohort and subtype summaries; no patient-level records.
fig,axs=plt.subplots(2,2,figsize=(8.0,6.6),gridspec_kw={'hspace':.62,'wspace':.65})
ax=axs[0,0];names=['Epi','Endo','Myeloid','T','Stromal','NK','B','Mast'];vals=[85473,40770,29503,14853,8626,7954,4854,542]
ax.barh(names[::-1],vals[::-1],color=blue);ax.set_xlabel('Cell count');panel(ax,'A','肺细胞图谱（20名供者；描述性）')
ax=axs[0,1];x=np.arange(3);ax.bar(x,[2,8,0],color=blue,label='Ctrl');ax.bar(x,[0,2,7],bottom=[2,8,0],color=red,label='CWP');ax.set_xticks(x,['G135','G227','G288']);ax.set_ylabel('Donors');ax.legend(frameon=False);panel(ax,'B','疾病与来源分布（19名供者）')
ax=axs[1,0];names=['CCR7','GZMA','IL7R','Treg','CD4state'];v=[927,1893,670,272,22];w=[245,661,709,328,796];ax.barh(names,v,color=blue,label='Ctrl');ax.barh(names,w,left=v,color=red,label='CWP');ax.invert_yaxis();ax.set_xlabel('Cell count');panel(ax,'C','CD4 亚群细胞覆盖')
ax=axs[1,1];heat(ax,[[6,3],[7,5],[7,9],[5,5],[1,6]],['Ctrl','CWP'],names,10,fmt='.0f',label='Donors ≥20 cells');panel(ax,'D','亚群供者覆盖（每端 ≥20 细胞）');save(fig,'fig1_population')

fig,axs=plt.subplots(2,2,figsize=(8.4,7.2),gridspec_kw={'hspace':.63,'wspace':.75})
terms=['KEGG_OXIDATIVE_PHOSPHORYLATION','GOBP_OXIDATIVE_PHOSPHORYLATION','GOBP_REGULATION_OF_EXECUTION_PHASE_OF_APOPTOSIS','GOBP_EXECUTION_PHASE_OF_APOPTOSIS','KEGG_APOPTOSIS','KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY']
short=['OXPHOS (KEGG)','OXPHOS (GO)','ApoExec regulation','ApoExec','Apoptosis (KEGG)','TCR (KEGG)']
g=df('human_GSEA_focused');s=df('human_GSVA_focused').rename(columns={'pathway_id':'ID'})
for ax,dd,val,pcol,qcol,letter,title,lim in [(axs[0,0],g,'NES','pvalue','p.adjust','A','人类通路排序富集',3),(axs[0,1],s,'logFC','P.Value','adj.P.Val','B','人类供者通路评分差异',.7)]:
 arr=[];marks=[]
 for t in terms:
  a=[];m=[]
  for scope in ['Primary','Core']:
   r=dd[(dd['ID']==t)&(dd.analysis==scope)].iloc[0];a.append(r[val]);m.append(sig(r[pcol],r[qcol]))
  arr.append(a);marks.append(m)
 heat(ax,arr,['全体CD4','核心CD4'],short,lim,marks,label='NES' if letter=='A' else 'Δ GSVA');panel(ax,letter,title)
ax=axs[1,0];d=df('human_mitochondrial_pathway_stability').query("analysis=='Core'");cols=['NES_Full','NES_NoCytoRP','NES_NoRP+MT'];xx=np.arange(3)
for _,r in d.iterrows():ax.plot(xx,[r[c] for c in cols],marker='o',ms=3,label=r.pathway_short)
ax.set_xticks(xx,['完整排序','去胞质核糖体','再去线粒体\n核糖体及编码基因']);ax.set_ylabel('NES');ax.legend(fontsize=6,ncol=2,frameon=False,loc='lower left');ax.set_ylim(1.8,3.0);panel(ax,'C','核心 CD4 呼吸链富集的敏感性')
ax=axs[1,1];m=df('mouse_GSEA_focused');arr=[];marks=[]
for t in ['GO:0006119','GO:0097194','GO:0042110']:
 a=[];ss=[]
 for scope in ['Core','Naive','ActCD4','Treg']:
  r=m[(m.pathway==t)&(m.specification==scope+'_Full_Standard_GO_BP')].iloc[0];a.append(r.NES);ss.append(sig(r.pval,r.padj))
 arr.append(a);marks.append(ss)
heat(ax,arr,['Core','Naive','ActCD4','Treg'],['OXPHOS','ApoExec','T-cell activation'],2,marks,label='NES');panel(ax,'D','小鼠通路富集（3+3文库）');save(fig,'fig2_pathways')

fig,axs=plt.subplots(2,2,figsize=(8.4,7.3),gridspec_kw={'hspace':.55,'wspace':.58})
mods=['OXPHOS','MitoISR','MitoROS','Mitophagy','Inflam','PyroExec','ApoExec']
for ax,model,letter,title in [(axs[0,0],'Source','A','表达匹配评分：来源校正'),(axs[0,1],'SourceDep','B','表达匹配评分：来源与深度校正')]:
 d=df('human_all_12_programs');arr=[];marks=[]
 for mod in mods:
  a=[];ms=[]
  for sc in ['Primary','Core','Subtype:CD4-CCR7','Subtype:CD4-GZMA','Subtype:CD4-IL7R','Subtype:Treg']:
   r=d[(d['filter']=='All')&(d.model==model)&(d.scope==sc)&(d.metric=='score_'+mod)].iloc[0];a.append(r.coefficient);ms.append(sig(r.p_nominal,r.FDR))
  arr.append(a);marks.append(ms)
 heat(ax,arr,['全体','核心','CCR7','GZMA','IL7R','Treg'],mods,.15,marks,label='Δ score');panel(ax,letter,title)
ax=axs[1,0];d=df('human_UCell_module_group_tests').query("analysis=='Core'");mods2=['MitoResp','MitoISR','MitoROS','Mitophagy','PyroCore','PyroSens','PyroExec','PyroNC','GSDMAlt'];arr=[];marks=[]
for mod in mods2:
 rr=[d[(d.module==mod)&(d.model==mm)].iloc[0] for mm in ['SourceAdj','QCAdj']];arr.append([r.estimate for r in rr]);marks.append([sig(r.p_value,r.fdr_bh) for r in rr])
heat(ax,arr,['来源校正','来源与质量校正'],mods2,.035,marks,fmt='.3f',label='Δ UCell');panel(ax,'C','核心 CD4 的 UCell 疾病差异')
ax=axs[1,1];d=df('mouse_UCell_score_group_tests');arr=[];marks=[]
for mod in ['OxPhos','MitoOrg','OxStress','Inflam','PyroResp','InflamRsp','DDR','PyroExec']:
 rr=[d[(d.module==mod)&(d['set']==sc)].iloc[0] for sc in ['Core','Naive','ActCD4','Treg']];arr.append([r.difference for r in rr]);marks.append([sig(r.PValue,r.FDR) for r in rr])
heat(ax,arr,['Core','Naive','ActCD4','Treg'],['OXPHOS','MitoOrg','OxStress','Inflam','PyroResp','InflamRsp','DDR','PyroExec'],.025,marks,fmt='.3f',label='Δ UCell');panel(ax,'D','小鼠 UCell 评分（精确置换）');save(fig,'fig3_programs')

fig,axs=plt.subplots(2,2,figsize=(8.2,7.0),gridspec_kw={'hspace':.38,'wspace':.35})
d=de('network_display_edges')
for ax,(sp,tf,letter) in zip(axs.flat,[('Human','IKZF2','A'),('Human','STAT1','B'),('Mouse','IKZF2','C'),('Mouse','STAT1','D')]):
 dd=d[(d.species==sp)&(d.tf==tf)].sort_values('weight',ascending=False).head(12);angles=np.linspace(0,2*np.pi,len(dd),endpoint=False);xy=np.c_[np.cos(angles),np.sin(angles)]
 for (x,y),(_,r) in zip(xy,dd.iterrows()):
  ax.plot([0,x],[0,y],color='#A1BAC4',lw=.7+2*r.weight/dd.weight.max(),zorder=1);ax.scatter(x,y,s=90,color='#DFEAF0',edgecolor=blue,lw=.5,zorder=2);ax.text(x*1.20,y*1.20,str(r.target),ha='center',va='center',fontsize=7)
 ax.scatter(0,0,s=1100,color=blue if tf=='IKZF2' else red,zorder=3);ax.text(0,0,tf,color='white',ha='center',va='center',fontsize=10,fontweight='bold');ax.set_xlim(-1.5,1.5);ax.set_ylim(-1.5,1.5);ax.axis('off');panel(ax,letter,('人类' if sp=='Human' else '小鼠')+' '+tf+' 推断调控网络')
save(fig,'fig4_networks')

fig,axs=plt.subplots(2,2,figsize=(8.4,7.4),gridspec_kw={'hspace':.7,'wspace':.65})
ax=axs[0,0];d=de('human_TF_disease').query("readout=='RNA' and scope=='CD4-GZMA' and tf=='IKZF2'");xx=np.arange(3)
for j,mo in enumerate(['Base','Source','SourceDep']):
 r=d[d.model==mo].iloc[0];ax.errorbar(j,r.effect,yerr=[[r.effect-r.CI_low],[r.CI_high-r.effect]],fmt='o',color=blue,capsize=3);ax.text(j,r.CI_high+.002,'P='+format(r.P,'.3g'),ha='center',fontsize=7)
ax.axhline(0,color=grey,lw=.7);ax.set_xticks(xx,['基础','来源校正','来源与深度校正']);ax.set_ylabel('RNA group effect (95% CI)');panel(ax,'A','人类 GZMA 亚群 IKZF2 RNA')
ax=axs[0,1];d=de('old_TF_group').query("species=='Mouse' and filter=='All' and variant=='DropOverlap' and model=='Group'");xx=np.arange(3)
for tf,col in [('IKZF2',blue),('STAT1',red)]:
 for scope,sty in [('Core','-'),('ActCD4','--')]:
  rr=d[(d.tf==tf)&(d.scope==scope)].sort_values('seed');ax.plot(xx,rr.coef,sty+'o',color=col,ms=4,label=tf+' / '+scope)
ax.axhline(0,color=grey,lw=.7);ax.set_xticks(xx,['17','29','43']);ax.set_xlabel('Random seed');ax.set_ylabel('Δ regulon AUC');ax.legend(fontsize=6,frameon=False);panel(ax,'B','小鼠调控子疾病差异')
d=df('all_TF_program_associations').query("species=='Human' and filter=='All' and variant=='DropOverlap' and scope=='Core' and model=='SourceDep'")
for ax,tf,letter in [(axs[1,0],'IKZF2','C'),(axs[1,1],'STAT1','D')]:
 arr=[];marks=[]
 for mod in mods:
  rr=[d[(d.tf==tf)&(d.seed==seed)&(d.feature=='score_'+mod)].iloc[0] for seed in [17,29,43]];arr.append([r.rho for r in rr]);marks.append([sig(r.p,r.FDR) for r in rr])
 heat(ax,arr,['17','29','43'],mods,1,marks,label='Partial Spearman ρ');panel(ax,letter,'核心 CD4：'+tf+'–程序关联')
save(fig,'fig5_regulons')

fig,axs=plt.subplots(1,3,figsize=(8.6,4.9),gridspec_kw={'wspace':.6});d=de('bootstrap_effect_summary');modality=['RNA19','RNA17','AUC17','AUC29','AUC43']
# Exact modality and adjustment identifiers are taken from archived summaries.
print('BOOT IDs',d.modality.unique(),d.adjustment.unique())
for ax,axis,adjust,letter,title in [(axs[0],'IL18','B','A','IL18–IKZF2'),(axs[1],'TNF','B','B','TNF–STAT1'),(axs[2],'TNF','B_Treg','C','TNF–STAT1；加入Treg比例')]:
 dd=d[(d.scheme=='stratified')&(d.axis==axis)&(d.adjustment==adjust)&(d.edge=='ligand_TF')];
 for j,(_,r) in enumerate(dd.iterrows()):ax.errorbar(r.original_effect,j,xerr=[[r.original_effect-r.CI95_low],[r.CI95_high-r.original_effect]],fmt='o',color=blue if axis=='IL18' else red,capsize=3)
 ax.set_yticks(range(len(dd)),[str(v).replace('AUC_seed','AUC ').replace('AUCseed','AUC ') for v in dd.modality]);ax.invert_yaxis();ax.axvline(0,color=grey,lw=.8);ax.set_xlim(-1,1);ax.set_xlabel('Partial Spearman ρ (95% CI)');panel(ax,letter,title)
save(fig,'fig6_bootstrap')

fig,axs=plt.subplots(2,2,figsize=(8.4,7.2),gridspec_kw={'hspace':.63,'wspace':.63})
ax=axs[0,0];x=np.arange(6);ax.bar(x-.18,[7,3,7,10,0,4],width=.36,label='FDR targets',color=blue);ax.bar(x+.18,[178,126,157,195,88,188],width=.36,label='Nominal targets',color='#ADCBD8');ax.set_xticks(x,['全体','核心','CCR7','GZMA','IL7R','Treg']);ax.set_ylabel('Upregulated target genes');ax.legend(frameon=False,fontsize=7);panel(ax,'A','人类来源校正后的疾病靶集')
ax=axs[0,1];h=de('human_NN_main');hh=h[(h['filter']=='All')&(h.model=='Source')&(h.test=='Nominal')&(h.direction=='Up')];
if hh.empty:hh=h[(h['filter']=='All')&(h.model=='Source')&(h.test=='Nominal05')&(h.direction=='Up')]
for lig,col in [('IL15',blue),('TNF',red),('IL18',grey)]:
 dd=hh[hh.test_ligand==lig].drop_duplicates('receiver').set_index('receiver');ax.plot(range(4),[dd.loc[k,'aupr_corrected'] if k in dd.index else np.nan for k in ['Core','CD4-GZMA','CD4-IL7R','Treg']],marker='o',ms=4,label=lig,color=col)
ax.axhline(0,color='#222222',lw=.7);ax.set_xticks(range(4),['Core','GZMA','IL7R','Treg']);ax.set_ylabel('Corrected AUPR');ax.legend(frameon=False);panel(ax,'B','人类配体活性：名义上调靶集')
ax=axs[1,0];m=de('mouse_NN');dd=m[(m.receiver=='Core')&(m['filter']=='All')&(m.test=='FDR05')&(m.direction=='Up')&(m.omit=='Full')];dd=dd[dd.test_ligand.isin(['Cd86','Cd274','Il15','Tnf','Il18'])].drop_duplicates('test_ligand').sort_values('rank');ax.barh(dd.test_ligand,dd.aupr_corrected,color=[blue if a=='Il15' else grey for a in dd.test_ligand]);ax.invert_yaxis();ax.set_xlabel('Corrected AUPR');panel(ax,'C','小鼠核心 CD4：FDR上调靶集')
ax=axs[1,1];x=np.arange(3);ax.bar(x-.18,[9,7,10],.36,color=blue,label='Ctrl');ax.bar(x+.18,[5,8,8],.36,color=red,label='CWP');ax.set_xticks(x,['MacAM','MonoMac','MacWide']);ax.set_ylabel('Donors ≥20 sender cells');ax.legend(frameon=False);panel(ax,'D','人类巨噬发送端覆盖');save(fig,'fig7_communication')

fig,axs=plt.subplots(1,2,figsize=(8.0,3.9),gridspec_kw={'wspace':.6});d=de('S4_direction_replication').query("filter=='All' and scope=='ActCD4'")
for ax,tf,letter in [(axs[0],'IKZF2','A'),(axs[1],'STAT1','B')]:
 for seed,col in zip([17,29,43],[blue,red,grey]):
  rr=d[(d.feature==tf)&(d.seed==seed)].sort_values('day');ax.plot(rr.day,rr.delta,marker='o',color=col,label=str(seed))
 ax.axhline(0,color='#222222',lw=.7);ax.set_xticks([7,56]);ax.set_xlabel('Day');ax.set_ylabel('Δ transferred AUC');ax.legend(title='Seed',frameon=False);panel(ax,letter,'独立时间点的 '+tf+' 调控子评分')
save(fig,'figS1_transfer');(D/'figure_label_audit.json').write_text(json.dumps(labels,ensure_ascii=False),encoding='utf8')
print('Figures saved',len(list(O.glob('*.png'))))
