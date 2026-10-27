from pathlib import Path
import json, hashlib, shutil, html, zipfile
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import FancyArrowPatch, Rectangle
from scipy.cluster.hierarchy import linkage, leaves_list
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import Paragraph
from reportlab.lib.styles import ParagraphStyle
import fitz

ROOT=Path('C:/codex_test'); SRC=ROOT/'tmp/silicosis_v9/v4'; OLD=ROOT/'tmp/silicosis_v9/v8'; OUT=ROOT/'output/silicosis_v9'
for d in ['figures','plot_data','provenance','code','baseline','supplement']: (OUT/d).mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white','svg.fonttype':'none'})
BLUE='#28678a'; RED='#bd4d50'; TEAL='#178b82'; GRAY='#7a8692'; SCOPES=['Naive','ActCD4','Treg']; MODS=['PyroExec','Inflam','ApoExec','OXPHOS','MitoISR','MitoROS','Mitophagy']
used={}; audit={}; pages=[]
def rd(name):
 p=SRC/name; used[name]=hashlib.sha256(p.read_bytes()).hexdigest(); return pd.read_csv(p,sep='\t')
def dump(df,name): df.to_csv(OUT/'plot_data'/f'{name}.tsv',sep='\t',index=False,encoding='utf-8')
def save(fig,name,title,text,method):
 fig.savefig(OUT/'figures'/f'{name}.png',dpi=190,bbox_inches='tight'); fig.savefig(OUT/'figures'/f'{name}.svg',bbox_inches='tight'); fig.savefig(OUT/'figures'/f'{name}.pdf',bbox_inches='tight'); plt.close(fig)
 pages.append(dict(name=name,title=title,text=text,method=method))
def title(fig,head,sub):
 fig.suptitle(head,x=.04,ha='left',fontsize=19,fontweight='bold',color='#163347',y=.98); fig.text(.04,.925,sub,fontsize=10,color=GRAY)
def heat(ax,arr,rows,cols,vmin,vmax,label,cmap='RdBu_r'):
 im=ax.imshow(np.ma.masked_invalid(np.asarray(arr,dtype=float)),aspect='auto',cmap=cmap,vmin=vmin,vmax=vmax)
 ax.set_yticks(range(len(rows)),rows); ax.set_xticks(range(len(cols)),cols,rotation=55,ha='right'); ax.tick_params(length=0)
 ax.set_facecolor('#e5e7eb'); return im
def stars(ax,df,rows,cols,rcol,ccol,p='P',q='BH'):
 lookup={(r[rcol],r[ccol]):r for _,r in df.iterrows()}
 for i,r in enumerate(rows):
  for j,c in enumerate(cols):
   z=lookup.get((r,c))
   if z is not None and z[p]<.05: ax.text(j,i,'*' if z[q]>=.05 else '**',ha='center',va='center',fontsize=9,color='black')

lr=rd('02_sender/all_prior_LR_expression_coverage.tsv.gz')
rank=rd('02_sender/source_specific_activity_ranks.tsv')
act=rd('03_receiver/TF_all_regulon_sample_activity.tsv.gz')
effects=rd('03_receiver/all_TF_activity_group_effects.tsv')
assoc=rd('03_receiver/all_TF_program_associations.tsv.gz')
prior=rd('03_receiver/focal_TF_prior_targets_and_modules.tsv')
de={s:rd(f'03_receiver/{s}_DE_all_genes.tsv.gz') for s in ['Core']+SCOPES}
sde={s:rd(f'02_sender/{s}_DE_all_genes.tsv.gz') for s in ['AMlike','AllMac']}
chain=rd('04_chains/all_sender_receiver_associations.tsv.gz')
complexes=rd('03_receiver/receptor_complex_codetection.tsv')

# 1. Focused multilayer communication network, separate source definitions.
pairs={'CXCL16':['CXCR6'],'TNF':['TNFRSF1A','TNFRSF1B'],'IL15':['IL2RB','IL2RG'],'IL18':['IL18R1','IL18RAP'],'C3':['C3AR1']}
ligands=list(pairs); receptors=sum(pairs.values(),[])
edges=lr[(lr.ligand.isin(ligands))&lr.scope.isin(SCOPES)].copy()
edges=edges[edges.apply(lambda r:r.receptor in pairs[r.ligand],axis=1)]
dump(edges,'F01_LR_coverage'); fig,axes=plt.subplots(1,2,figsize=(15,8)); title(fig,'01  Macrophage-to-CD4 communication architecture','Focused report candidates | separate AMlike and AllMac views | mouse silica libraries n=3')
for ax,sender in zip(axes,['AMlike','AllMac']):
 ax.set_xlim(-.25,3.4); ax.set_ylim(-.65,8.7); ax.axis('off'); ax.set_title(sender+(' (nested within AllMac)' if sender=='AMlike' else ' (inclusive macrophage pool)'),fontsize=12,pad=8)
 ly=dict(zip(ligands,np.linspace(7.6,.5,5))); ry=dict(zip(receptors,np.linspace(7.8,.2,len(receptors)))); sy=dict(zip(SCOPES,[6.7,4,1.3]))
 for x,lab in [(0,'LIGAND'),(1.5,'RECEPTOR'),(3,'CD4 STATE')]: ax.text(x,8.4,lab,ha='center',fontsize=10,color=GRAY)
 ed=edges[edges.sender==sender]
 for lig,rs in pairs.items():
  for rec in rs: ax.plot([0,1.5],[ly[lig],ry[rec]],color='#9aa9b3',lw=1,ls='--',zorder=0)
 for rec in receptors:
  for scope in SCOPES:
   sub=ed[(ed.receptor==rec)&(ed.scope==scope)]
   if sub.empty: continue
   z=sub.iloc[0]
   if z.receiver_case_detected>0:
    ax.plot([1.5,3],[ry[rec],sy[scope]],color=TEAL,alpha=.45,lw=.5+2.7*z.receiver_case_mean_pct/100,zorder=0)
 for lig,y in ly.items():
  zz=sde[sender].set_index('canonical'); fc=float(zz.loc[lig,'logFC']) if lig in zz.index else np.nan
  ax.scatter(0,y,s=540,c=[fc],cmap='RdBu_r',vmin=-2,vmax=2,edgecolor='white',zorder=3); ax.text(-.17,y,lig,ha='right',va='center',fontweight='bold',fontsize=10)
 for rec,y in ry.items(): ax.scatter(1.5,y,s=85,color='#e8eff4',edgecolor=BLUE,zorder=3); ax.text(1.5,y+.18,rec,ha='center',fontsize=8.5)
 for scope,y in sy.items(): ax.scatter(3,y,s=700,color='#e0f1ed',edgecolor=TEAL,zorder=3); ax.text(3,y,scope,ha='center',va='center',fontsize=9,fontweight='bold')
fig.text(.04,.045,'Ligand color: source log2FC [-2 blue, +2 red]. Dashed: prior pair. Teal width: receptor mean detection %.\nReceptor branches denote expression coverage, not communication probability; IL15/IL18 complexes require joint detection.',fontsize=10)
fig.subplots_adjust(top=.83,bottom=.15,wspace=.28)
save(fig,'F01_communication_network','图谱1｜巨噬细胞与CD4亚群的通讯结构','CXCL16、TNF、IL15、IL18和C3的配体受体候选被放回发送端与接收亚群结构中。AMlike和AllMac分开呈现，避免把包含关系误作独立细胞来源。','范围沿用v8重点候选；虚线为先验配体受体关系，绿色线宽为病例接收端平均检出比例，不代表通讯强度。单亚基检出不等于完整受体复合物。原始覆盖及样本数见F01数据表。')

# 2. Prediction landscape plus source effects and receptor-complex coverage.
up=rank[(rank.direction=='Up')&rank.scope.isin(SCOPES)]
top=up[up.source_specific_rank<=5].ligand.unique().tolist(); selected=sorted(set(top+ligands))
dump(up[up.ligand.isin(selected)],'F02_ligand_prediction')
fig=plt.figure(figsize=(15,9)); title(fig,'02  Predicted ligand activity and receptor availability','Union of top-5 Up-target ligands per source/state + five focal ligands; ranks are predictions, not P values')
gs=fig.add_gridspec(1,3,left=.09,right=.96,top=.84,bottom=.25,width_ratios=[1,1,1.6],wspace=.5)
for k,sender in enumerate(['AMlike','AllMac']):
 ax=fig.add_subplot(gs[k]); a=up[up.sender==sender]; mat=a.pivot(index='ligand',columns='scope',values='source_specific_rank').reindex(index=selected,columns=SCOPES)
 arr=-np.log10(mat.to_numpy(dtype=float)); im=heat(ax,arr,selected,SCOPES,-2.5,0,'rank','viridis')
 for i in range(len(selected)):
  for j in range(3):
   v=mat.iloc[i,j]
   if np.isfinite(v): ax.text(j,i,str(int(v)),ha='center',va='center',fontsize=8,color='white' if v>12 else 'black')
 ax.set_title(sender+' | activity rank'); ax.tick_params(axis='y',labelsize=8)
ax=fig.add_subplot(gs[2]); comp=complexes[(complexes.group=='Silica')&complexes.scope.isin(SCOPES)].copy(); comp['row']=comp.ligand+' ['+comp.rule.str.upper()+']\n'+comp.components
cm=comp.groupby(['row','scope']).joint_pct.mean().unstack().reindex(columns=SCOPES); dump(comp,'F02_receptor_joint_detection'); im=heat(ax,cm,cm.index,SCOPES,0,100,'joint %','YlGnBu'); ax.set_title('Receptor availability (ALL / ANY)'); ax.tick_params(axis='y',labelsize=8)
for i in range(len(cm)):
 for j in range(3):
  if pd.notna(cm.iloc[i,j]): ax.text(j,i,f'{cm.iloc[i,j]:.1f}',ha='center',va='center',fontsize=9,color='white' if cm.iloc[i,j]>55 else 'black')
fig.colorbar(im,ax=ax,shrink=.5,label='Mean receptor-rule detection (%)')
fig.text(.05,.07,'Rank panels: bright = better prediction rank; gray = no evaluable ranking in the supplied table.\nRight: average across silica libraries with available state; separate source/receiver coverage is retained in plot data.\nALL: all listed components detected; ANY: either TNF receptor detected. Transcript coverage does not prove signaling.',fontsize=10)
save(fig,'F02_ligand_landscape','图谱2｜配体预测排名与受体完整性共同约束通讯解释','左侧展示各来源、各CD4状态预测排名前5的并集及既有重点配体。右侧ALL为同细胞全部组分检出，TNF的ANY为任一受体检出。两类证据分开编码，比较预测排序与接收端可用性。','排名基于既有NicheNet结果、接收端名义上调基因集合及来源表达资格。排名不是显著性检验，不同亚群的目标集不同。右侧为病例样本等权平均规则检出率；TNFRSF1A与TNFRSF1B为替代受体，不要求共检出。灰色为不可评估。')

# 3. Sample-level regulon atlas: choose variance, without significance selection.
a=act[(act.seed==17)&act.scope.isin(SCOPES)].copy(); a['col']=a.scope+' | '+a['sample']; wide=a.pivot(index='regulon',columns='col',values='mean_auc')
sel=wide.var(axis=1).nlargest(22).index.tolist(); sel=list(dict.fromkeys(sel+['Stat1(+)','Ikzf2(+)'])); sel=[x for x in sel if x in wide.index]
samples=['Saline_1','Saline_2','Saline_3','Silica_1','Silica_2','Silica_3']; actual=sorted(a['sample'].unique()); samples=[x for x in samples if x in actual]+[x for x in actual if x not in samples]
cols=[s+' | '+x for s in SCOPES for x in samples]; w=wide.reindex(index=sel,columns=cols); assert w.notna().all().all(),w.columns.tolist()
z=w.sub(w.mean(axis=1),axis=0).div(w.std(axis=1),axis=0).fillna(0); order=leaves_list(linkage(z,method='average')); z=z.iloc[order]; sel=z.index.tolist()
dump(a[a.regulon.isin(sel)],'F03_sample_AUC'); dump(effects[(effects.seed==17)&effects.regulon.isin(sel)],'F03_disease_effects')
fig=plt.figure(figsize=(15,9)); title(fig,'03  CD4 regulon activity across biological samples','22 highest-variance regulons + STAT1 and IKZF2 | seed 17, overlap-removed AUC | same six mice across states')
gs=fig.add_gridspec(1,2,left=.1,right=.94,top=.83,bottom=.2,width_ratios=[4,1],wspace=.35); ax=fig.add_subplot(gs[0]); im=heat(ax,z,sel,[x.split(' | ')[1].replace('Saline_','C').replace('Silica_','S') for x in cols],-2,2,'z')
for j,scope in enumerate(SCOPES):
 ax.text(j*6+2.5,-1.2,scope,ha='center',fontweight='bold'); ax.axvline(j*6-.5,color='white',lw=3)
 for k in range(6): ax.add_patch(Rectangle((j*6+k-.5,-.9),1,.3,color=BLUE if k<3 else RED,clip_on=False))
fig.colorbar(im,ax=ax,shrink=.5,label='Row z-score (AUC)')
ax=fig.add_subplot(gs[1]); ef=effects[(effects.seed==17)&effects.scope.isin(SCOPES)]; m=ef.pivot(index='regulon',columns='scope',values='delta_mean_auc').reindex(index=sel,columns=SCOPES); limit=max(abs(m.min().min()),abs(m.max().max())); im=heat(ax,m,[],SCOPES,-limit,limit,'delta'); stars(ax,ef,sel,SCOPES,'regulon','scope'); ax.set_title('Disease - control'); fig.colorbar(im,ax=ax,shrink=.5,label='Delta mean AUC')
fig.text(.05,.055,'C1-C3: saline; S1-S3: silica. Row scaling reveals patterns, not comparable absolute TF activity.\n* P<0.05; ** original BH<0.05. Selection uses variance, not disease P; biological n=3+3, not 18.\nExact 3-vs-3 group tests have coarse resolution; an unstarred activity pattern is descriptive.',fontsize=10)
fig.subplots_adjust(bottom=.18)
save(fig,'F03_TF_activity_atlas','图谱3｜CD4调控子活性呈现亚群与样本层次','调控子热图以独立小鼠为列，展示Naive、ActCD4、Treg的整体活性模式；右侧保留同一调控子的疾病效应。STAT1和IKZF2与其他高变调控子一起呈现，避免只展示两个候选造成缺乏背景。','按全部18个样本-亚群组合的AUC方差选前22个，再加入STAT1和IKZF2，去重后聚类；同一6只小鼠在各亚群重复出现，不构成18个独立样本。颜色为行标准化AUC；效应和P/BH使用原表。')

# 4. TF-program matrix, retain all values for displayed TFs.
b=assoc[(assoc.seed==17)&(assoc.model=='Pooled')&assoc.scope.isin(SCOPES)].copy(); counts=b.assign(hit=b.P<.05).groupby('x').hit.sum().sort_values(ascending=False,kind='stable'); selected=counts.head(18).index.tolist(); selected=list(dict.fromkeys(selected+['Stat1(+)','Ikzf2(+)'])); b=b[b.x.isin(selected)]; dump(b,'F04_TF_program_associations')
fig,axes=plt.subplots(1,3,figsize=(16,9)); title(fig,'04  Regulons connect to distinct transcriptional programs','18 TFs with most nominal links across three states + STAT1 and IKZF2 | full original BH families retained')
for ax,scope in zip(axes,SCOPES):
 ss=b[b.scope==scope]; m=ss.pivot(index='x',columns='y',values='rho').reindex(index=selected,columns=MODS); im=heat(ax,m,selected,MODS,-1,1,'rho'); stars(ax,ss,selected,MODS,'x','y'); ax.set_title(scope); ax.tick_params(labelsize=9)
fig.subplots_adjust(left=.09,right=.91,top=.84,bottom=.23,wspace=.46); cax=fig.add_axes([.94,.35,.012,.35]); fig.colorbar(im,cax=cax,label='Spearman rho')
fig.text(.04,.055,'* nominal P<0.05; ** BH<0.05. Unmarked cells remain visible. Selection is exploratory, not external validation.\nPyroExec: death-related transcript program; Inflam: inflammasome; ApoExec: apoptosis; MitoISR: mitochondrial stress.\nScores do not measure cell death or metabolic flux. Six biological samples; pooled relationships may reflect disease grouping.',fontsize=10)
save(fig,'F04_TF_function_map','图谱4｜TF与死亡、线粒体和自噬程序的关系具有状态特异性','将TF与7类功能程序按亚群并列展示。图中保留入选TF的全部关系，包括未达名义显著的单元格，便于比较同一TF在不同亚群和程序间的方向变化。','探索性展示按3个亚群内P<0.05关联数排序的前18个TF，并加入STAT1和IKZF2。星号沿用原始P与完整家族BH；未对展示子集重新校正。此排序用于展示候选背景，不用于确立新的实验优先级。')

# 5. Prior target modules; absence is explicit, not fabricated.
pt=prior[prior.module.notna()].copy(); dump(pt,'F05_prior_module_edges'); targets=sorted(pt.target.unique()); modules=sorted(pt.module.unique())
fig=plt.figure(figsize=(15,8)); title(fig,'05  STAT1 prior targets overlap multiple functional modules','Prior-derived target-module map | expression overlay from mouse pseudobulk | not a reconstructed causal pathway')
gs=fig.add_gridspec(1,2,left=.05,right=.94,top=.84,bottom=.19,width_ratios=[2,1],wspace=.38); ax=fig.add_subplot(gs[0]); ax.set_xlim(-.2,4.5); ax.set_ylim(-1,8); ax.axis('off'); ty=dict(zip(targets,np.linspace(6.7,.6,len(targets)))); my=dict(zip(modules,np.linspace(6.7,.6,len(modules))))
ax.text(0,7.5,'TF',color=GRAY); ax.text(1.7,7.5,'PRIOR TARGET',color=GRAY); ax.text(3.4,7.5,'MODULE',color=GRAY)
for gene,y in ty.items():
 ax.plot([0,1.8],[4,y],ls='--',color=BLUE,alpha=.65); ax.scatter(1.8,y,s=330,color='#e5eef4',edgecolor=BLUE,zorder=3); ax.text(1.8,y,gene,ha='center',va='center',fontsize=9)
for _,r in pt.iterrows(): ax.plot([1.8,3.7],[ty[r.target],my[r.module]],color=TEAL,alpha=.65,lw=1.8)
ax.scatter(0,4,s=1300,color=BLUE,zorder=3); ax.text(0,4,'STAT1',ha='center',va='center',color='white',fontweight='bold')
for mod,y in my.items(): ax.text(3.7,y,mod,ha='center',va='center',bbox=dict(boxstyle='round,pad=.4',fc='#e2f2ed',ec=TEAL),fontsize=10)
ax.text(0,-.45,'IKZF2: no target-module overlap in this supplied prior table.\nGSDMD is not added as a direct STAT1 target.',fontsize=10,color=GRAY)
ax=fig.add_subplot(gs[1]); vals=pd.DataFrame(index=targets,columns=['Core']+SCOPES,dtype=float); records=[]
for scope,df in de.items():
 for gene in targets:
  rr=df[df.canonical==gene]
  if len(rr):
   r=rr.iloc[0]; vals.loc[gene,scope]=r.logFC; records.append(dict(gene=gene,scope=scope,logFC=r.logFC,P=r.PValue,BH=r.FDR))
dd=pd.DataFrame(records); dump(dd,'F05_target_DE'); im=heat(ax,vals,targets,vals.columns,-2,2,'logFC'); stars(ax,dd,targets,vals.columns,'gene','scope'); ax.set_title('Target RNA disease effect'); fig.colorbar(im,ax=ax,shrink=.5,label='log2FC (clipped +/-2)')
fig.text(.04,.045,'Dashed blue: prior TF-target relation. Solid teal: gene-set membership, not signaling.\n* P<0.05; ** original BH<0.05; gray: not tested. Prior confidence and mode of regulation are preserved in plot data.',fontsize=10)
save(fig,'F05_target_module_network','图谱5｜STAT1先验靶基因覆盖多个功能模块','现有先验记录中，STAT1与AIM2、BAX、EIF2AK1、FOXP3、HMOX1、IL1B和IL2RA等模块成员存在连接。右侧显示对应基因在不同CD4范围的疾病效应，区分先验联系与本队列表达支持。','虚线为先验TF-靶基因边，实线为基因集归属，均不代表本项目验证的因果关系。现有表中没有IKZF2与所列模块的重合记录，不据此判定IKZF2无调控作用；GSDMD未被补画为STAT1直接靶基因。')

# 6. Source-to-TF bridge with adjusted results next to pooled results.
c=chain[(chain.seed==17)&chain.scope.isin(SCOPES)&chain.x.isin(ligands)&chain.y.eq('AUC::Stat1(+)')].copy(); dump(c,'F06_sender_STAT1_bridge')
fig,axes=plt.subplots(1,2,figsize=(15,7)); title(fig,'06  The ligand-to-STAT1 bridge depends on source and model','Same-sample association | pooled vs disease-adjusted permutation model | displayed regardless of significance')
for ax,sender in zip(axes,['AMlike','AllMac']):
 ss=c[c.sender==sender].copy(); ss['col']=ss.scope+'\n'+ss.model.map({'Pooled':'pooled','ConditionAdjusted':'adjusted'}); cols=[s+'\n'+m for s in SCOPES for m in ['pooled','adjusted']]; mat=ss.pivot(index='x',columns='col',values='rho').reindex(index=ligands,columns=cols); im=heat(ax,mat,ligands,cols,-1,1,'rho'); ax.set_title(sender); stars(ax,ss,ligands,cols,'x','col')
 for i,lig in enumerate(ligands):
  for j,col in enumerate(cols):
   rr=ss[(ss.x==lig)&(ss.col==col)]
   if len(rr): ax.text(j,i+.27,f'P={rr.iloc[0].P:.3g}',ha='center',va='center',fontsize=7,color='white' if abs(rr.iloc[0].rho)>.7 else 'black')
 for v in [1.5,3.5]: ax.axvline(v,color='white',lw=3)
 ax.tick_params(axis='x',labelrotation=0,labelsize=9)
fig.subplots_adjust(left=.09,right=.89,top=.81,bottom=.23,wspace=.3); cax=fig.add_axes([.93,.3,.013,.4]); fig.colorbar(im,cax=cax,label='Association coefficient')
fig.text(.05,.06,'* P<0.05; ** BH<0.05; full original association families. Gray: unavailable. n=6 matched mice.\nAdjusted model uses within-condition permutations (36); pooled uses 720 permutations.\nDifferences between AMlike and AllMac prohibit treating a pooled source signal as an AMlike-specific mechanism.',fontsize=10)
save(fig,'F06_bridge_evidence','图谱6｜配体与STAT1的联系依赖来源定义和组别调整','把同一样本的配体表达与接收端STAT1活性联系展开为完整矩阵，合并分析与疾病组别调整并排显示。AMlike和AllMac中的差异直接可见，不再由单个相关系数代表整条通讯链。','每格沿用既有模型的系数、P和BH，不新增推断检验。组别调整模型使用36种组内置换，合并模型使用720种置换。该图评价样本间联系，不证明受体激活或配体直接驱动STAT1。')

# Preserve original evidence and assemble new atlas plus complete v8 report.
for p in (OLD/'supplement').iterdir(): shutil.copy2(p,OUT/'supplement'/p.name)
oldpdf=next(OLD.glob('*.pdf')); shutil.copy2(oldpdf,OUT/'baseline'/'report_v8.pdf')
shutil.copytree(OLD/'assets',OUT/'baseline'/'assets',dirs_exist_ok=True)
shutil.copy2(next(OLD.glob('*.html')),OUT/'baseline'/'report_v8.html')
pdfmetrics.registerFont(TTFont('CJK','C:/Windows/Fonts/ARIALUNI.ttf'))
style=ParagraphStyle('body',fontName='CJK',fontSize=11,leading=18,textColor='#243746')
small=ParagraphStyle('small',parent=style,fontSize=9,leading=14)
W,H=842,595
def para(cv,text,x,y,width,sty=style):
 p=Paragraph(html.escape(text),sty); _,h=p.wrap(width,800); p.drawOn(cv,x,y-h); return y-h
tmp=ROOT/'tmp/silicosis_v9/atlas.pdf'; cv=canvas.Canvas(str(tmp),pagesize=(W,H)); cv.setTitle('矽肺CD4研究：细胞通讯与TF图谱补充整合版 v9')
def header(n,label):
 cv.setFillColorRGB(.075,.18,.25); cv.rect(0,H-13,W,13,fill=1,stroke=0); cv.setFont('CJK',9); cv.drawString(35,20,'SILICOSIS / COMMUNICATION & REGULON ATLAS / v9'); cv.drawRightString(W-35,20,f'新增图谱 {n} / 7'); cv.setFont('CJK',19); cv.drawString(35,H-51,label)
header(1,'矽肺CD4研究｜通讯与TF图谱补充整合版')
y=para(cv,'本版用途：补足通讯结构、调控状态与功能模块，供后续研究走向讨论。基于现有冻结统计表重绘，未改变候选筛选阈值、样本名单或已有统计结果。',38,510,760)
y=para(cv,'新增六组图谱；后附v8完整报告。原版保留不覆盖，原统计图、患者方向差异、细胞身份与背景RNA问题均保留。',38,y-12,760)
for i,p in enumerate(pages,1): y=para(cv,f'{i}. '+p['title'].split('｜',1)[1],50,y-14,720)
y=para(cv,'阅读要点：配体受体预测、TF先验靶基因、功能程序相关和患者定位是不同证据。图中关联不升级为因果；排名不替代P值，名义显著不等同BH显著。',38,y-25,760)
para(cv,'数据基准：integrated_macrophage_CD4_20261015_v4；既有报告：client_report_20261015_v8。新增图件为英文标签，中文解读；所有绘图取值、来源SHA256及构建脚本随附。',38,y-16,760,small); cv.showPage()
for i,p in enumerate(pages,2):
 header(i,p['title']); img=fitz.Pixmap(str(OUT/'figures'/f"{p['name']}.png")); iw,ih=img.width,img.height; maxw,maxh=770,355; scale=min(maxw/iw,maxh/ih); dw,dh=iw*scale,ih*scale
 cv.drawImage(str(OUT/'figures'/f"{p['name']}.png"),(W-dw)/2,151+(maxh-dh)/2,width=dw,height=dh)
 y=para(cv,p['text'],38,139,766,small); para(cv,p['method'],38,y-7,766,small); cv.showPage()
cv.save()
doc=fitz.open(tmp); baseline=fitz.open(oldpdf); original_pages=len(baseline); doc.insert_pdf(baseline); doc.set_toc([[1,'新增：细胞通讯与TF图谱',1]]+[[2,p['title'],i+2] for i,p in enumerate(pages)]+[[1,'既有v8完整报告（原页保留）',8]])
pdf=OUT/'silicosis_communication_TF_v9.pdf'; doc.save(pdf); doc.close()
sections=['<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>矽肺通讯与TF图谱 v9</title><style>body{max-width:1300px;margin:40px auto;font:17px/1.8 Arial,"Microsoft YaHei",sans-serif;color:#163347;padding:20px}img{width:100%}section{margin:50px 0}a{color:#178b82}</style><h1>矽肺CD4研究｜通讯与TF图谱补充整合版 v9</h1><p>六组新增图谱，复用冻结统计结果，供后续研究走向讨论。</p><p><a href="baseline/report_v8.html">阅读既有v8完整报告</a> · <a href="silicosis_communication_TF_v9.pdf">下载整合PDF</a></p>']
for p in pages: sections.append(f'<section><h2>{p["title"]}</h2><img src="figures/{p["name"]}.png"><p>{p["text"]}</p><p>{p["method"]}</p></section>')
sections.append('</html>'); (OUT/'report_v9.html').write_text('\n'.join(sections),encoding='utf-8')
(OUT/'provenance'/'source_sha256.json').write_text(json.dumps(used,indent=2,ensure_ascii=False),encoding='utf-8')
(OUT/'provenance'/'figure_notes.json').write_text(json.dumps(pages,indent=2,ensure_ascii=False),encoding='utf-8')
shutil.copy2(__file__,OUT/'code'/'build_v9.py')
audit.update(figures=len(pages),new_pages=7,baseline_pages=original_pages,total_pages=original_pages+7,baseline_sha256=hashlib.sha256(oldpdf.read_bytes()).hexdigest(),TF_regulons=int(act[act.seed==17].regulon.nunique()),new_inferential_tests=0)
(OUT/'provenance'/'build_validation.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print(json.dumps(audit)); print('BUILD COMPLETE',pdf)
