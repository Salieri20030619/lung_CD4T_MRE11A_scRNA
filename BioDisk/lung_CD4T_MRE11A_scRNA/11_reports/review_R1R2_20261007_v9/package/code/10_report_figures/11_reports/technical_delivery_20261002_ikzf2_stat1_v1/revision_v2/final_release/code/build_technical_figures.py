from pathlib import Path
import json,shutil,hashlib
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
P=Path('technical_delivery'); O=P/'figures';O.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white','pdf.fonttype':42,'axes.titleweight':'bold'})
COL={'All':'#21618c','BG1':'#c56332','IKZF2':'#176b87','STAT1':'#ba4b43'}
caps=[]
def rd(n,folder='data'):return pd.read_csv(P/folder/n,sep='\t')
def pv(x):return 'NA' if pd.isna(x) else ('<0.0001' if x<.0001 else f'{x:.4f}')
def flag(p,q=np.nan):
 if pd.isna(p):return 'No P'
 return 'F' if pd.notna(q) and q<.05 else ('N' if p<.05 else '-')
def done(fig,name,title,caption,source):
 fig.savefig(O/(name+'.png'),dpi=220,bbox_inches='tight');fig.savefig(O/(name+'.pdf'),bbox_inches='tight')
 plt.close(fig);caps.append(dict(figure=name,title=title,caption_zh=caption,source=source))
# F04: original RNA count models, no re-fit and no cell-wise inference.
de=rd('existing_RNA_DE_focused.tsv')
comparisons=[('Mouse','Core','All'),('Mouse','ActCD4','All'),('Mouse','Core','BG1'),('Mouse','ActCD4','BG1'),('Mouse','Naive','All'),('Mouse','Naive','BG1'),('Human','Primary','All'),('Human','Core','All')]
fig,(ax,tabax)=plt.subplots(1,2,figsize=(13,5.8),gridspec_kw={'width_ratios':[1.1,1.4]})
table=[]
for i,(sp,scope,arm) in enumerate(comparisons):
 row=[f'{sp} / {scope} / {arm}']
 for j,tf in enumerate(['IKZF2','STAT1']):
  z=de[(de.species==sp)&(de.scope==scope)&(de['filter']==arm)&(de.tf==tf)]
  if len(z):
   r=z.iloc[0];y=i+(-.11 if j==0 else .11);formal=r.PValue<.05
   ax.scatter(r.logFC,y,s=65,marker='o' if j==0 else 's',facecolor=COL[tf] if formal else 'white',edgecolor=COL[tf],linewidth=1.5,zorder=3)
   row.append(f'{pv(r.PValue)} / {pv(r.FDR)} [{flag(r.PValue,r.FDR)}]')
  else:row.append('Not in tested-gene table')
 table.append(row)
ax.axvline(0,color='#999',ls='--',lw=1);ax.set_yticks(range(len(comparisons)),[f'{a} / {b} / {c}' for a,b,c in comparisons]);ax.invert_yaxis()
ax.set_xlabel('RNA log2 fold change: disease - control');ax.set_title('Original RNA count models')
ax.legend(handles=[Line2D([0],[0],marker='o',color='none',markerfacecolor=COL['IKZF2'],markeredgecolor=COL['IKZF2'],label='IKZF2'),Line2D([0],[0],marker='s',color='none',markerfacecolor=COL['STAT1'],markeredgecolor=COL['STAT1'],label='STAT1')],loc='lower right',frameon=False)
tabax.axis('off');t=tabax.table(cellText=table,colLabels=['Original model scope','IKZF2 P / FDR','STAT1 P / FDR'],loc='center',cellLoc='left',colWidths=[.36,.32,.32]);t.auto_set_font_size(False);t.set_fontsize(8.7);t.scale(1,2.1)
for (r,c),cell in t.get_celld().items():
 cell.set_edgecolor('#dddddd')
 if r==0:cell.set_facecolor('#eaf1f4');cell.set_text_props(weight='bold')
fig.suptitle('F04 | RNA abundance evidence: nominal P and FDR are separate',fontsize=15,y=1.01)
fig.text(.5,-.01,'Filled marker: nominal P < 0.05. Open: not nominal. N: nominal only; F: FDR < 0.05. Missing tested genes are not zeros.',ha='center',fontsize=9)
fig.tight_layout(w_pad=2)
done(fig,'F04_RNA_existing_DE','两候选RNA：原有计数模型',
'直接展示原edgeR结果，未重新拟合。Mouse：GSE314943，独立文库3 Saline/3 Silica，原全库offset/TMM，All/BG1及原Core/ActCD4/Naive范围分别保留。Human：原Primary/Core源校正全基因模型；其统计样本范围沿用原模型，不把下一图的SCENIC评分细胞子集当作重新拟合样本。N=名义P<0.05但FDR未通过；F=FDR<0.05；空心点=未达名义标准；缺失基因=未出现在该检验表，不能当零表达。RNA差异与regulon活性不是同一指标。','data/existing_RNA_DE_focused.tsv')
# F04b: purely descriptive detection and RNA levels in frozen SCENIC scoring cells.
rna=rd('RNA_sample_summary.tsv');fig,axes=plt.subplots(2,4,figsize=(14,7),sharex='col')
for col,(sp,scope,tf) in enumerate([(s,'ActCD4' if s=='Mouse' else 'Core',t) for s in ['Mouse','Human'] for t in ['IKZF2','STAT1']]):
 z=rna[(rna.species==sp)&(rna.scope==scope)&(rna.tf==tf)&rna.eligible20]
 groups=['Saline','Silica'] if sp=='Mouse' else ['Control','CWP']
 for row,measure in enumerate(['mean_CP10k','pct']):
  ax=axes[row,col]
  for arm,off in [('All',-.12),('BG1',.12)]:
   for j,g in enumerate(groups):
    v=z[(z['filter']==arm)&(z.group==g)].sort_values('sample')
    jitter=np.linspace(-.035,.035,len(v)) if len(v)>1 else np.zeros(len(v))
    ax.scatter(j+off+jitter,v[measure],s=27,color=COL[arm],alpha=.85,label=arm if j==0 else None)
  ax.set_xticks([0,1],groups);ax.grid(axis='y',color='#eee');ax.set_axisbelow(True)
  if measure=='pct':ax.set_ylim(-3,103);ax.set_ylabel('Detected cells (%)')
  else:ax.set_ylabel('Sample mean RNA (CP10k)');ax.set_title(f'{sp} {scope} | {tf}\n'+'; '.join(f'{g}: n={z[(z["filter"]=="All")&(z.group==g)]["sample"].nunique()}' for g in groups),fontsize=10)
  if col==0 and row==0:ax.legend(frameon=False,fontsize=9)
fig.suptitle('F04b | RNA levels and detection in frozen scoring cells (descriptive)',fontsize=15,y=1.02);fig.tight_layout()
done(fig,'F04b_RNA_level_detection','两候选RNA：样本水平表达与检出',
'每点为一个独立文库/受试者；样本/亚群均沿用至少20细胞要求。RNA值来自冻结SCENIC1 score.loom，按原定义expm1还原log1p(CP10k)，以原始全库总数为分母后取样本均值；检出为该样本该范围内表达>0的细胞百分比。仅作技术汇总，无新检验；两过滤臂为同一样本的敏感性视图，不是独立重复。鼠ActCD4和人Core不是等同亚群，不合并绝对值。','data/RNA_sample_summary.tsv; provenance/S1_methods.txt')
# F05: effect/LOO panels, baseline and stronger original adjustment shown together.
human=rd('S2_TF_group_effects.tsv');human=human[(human.species=='Human')&(human.scope=='Core')&(human.variant=='DropOverlap')]
mouse=rd('S3_mouse_TF_identity_adjustment.tsv','source_tables');mouse['tf']=mouse.tf.str.upper()
fig,axs=plt.subplots(2,2,figsize=(12,8))
pairs=[(arm,seed) for arm in ['All','BG1'] for seed in [17,29,43]]
for r,sp in enumerate(['Mouse','Human']):
 for c,tf in enumerate(['IKZF2','STAT1']):
  ax=axs[r,c];d=mouse if sp=='Mouse' else human
  models=['Group','DepthTreg'] if sp=='Mouse' else ['Source','SourceDep']
  for i,(arm,seed) in enumerate(pairs):
   for j,m in enumerate(models):
    z=d[(d.tf==tf)&(d['filter']==arm)&(d.seed==seed)&(d.model==m)]
    if z.empty:continue
    q=z.iloc[0];y=i+(-.12 if j==0 else .12)
    ax.plot([q.loo_min,q.loo_max],[y,y],color=COL[arm],alpha=.45,lw=1.7)
    ax.scatter(q.coef,y,s=40,marker='o' if j==0 else 's',facecolor=COL[arm] if q.p<.05 else 'white',edgecolor=COL[arm],zorder=3)
  ax.axvline(0,color='#888',ls='--',lw=1);ax.set_yticks(range(6),[f'{a} s{s}' for a,s in pairs]);ax.invert_yaxis()
  ax.set_title(f'{sp} '+('ActCD4' if sp=='Mouse' else 'Core')+f' | {tf}',fontsize=12)
  ax.set_xlabel('Disease coefficient (AUCell units)')
  ax.legend(handles=[Line2D([0],[0],marker='o',color='#777',label=models[0]),Line2D([0],[0],marker='s',color='#777',label=models[1])],frameon=False,fontsize=8)
fig.suptitle('F05 | Regulon activity: original effects and leave-one-out ranges',fontsize=15,y=1.01);fig.tight_layout()
done(fig,'F05_regulon_activity_effects','Regulon活性：效应与敏感性',
'DropOverlap既有模型，三seed分别展示。鼠ActCD4：Group与DepthTreg（深度+TregNoTF调整）；人Core：Source与SourceDep（来源+接收细胞深度调整）。线段为逐一移除样本后的系数范围，不是置信区间；实心仅表示该完整样本模型名义P<0.05，空心未达标。各P/FDR及df完整保留于原表。每个seed是算法敏感性，不是生物学重复；AUCell变化不得改称RNA上调。','source_tables/S3_mouse_TF_identity_adjustment.tsv; data/S2_TF_group_effects.tsv')
# F06: display-only selection from full original target table.
edges=rd('TF_network_all_edges.tsv');display=[]
fig,axs=plt.subplots(2,2,figsize=(11,8))
for row,sp in enumerate(['Mouse','Human']):
 for col,tf in enumerate(['IKZF2','STAT1']):
  ax=axs[row,col]
  e=edges[(edges.species==sp)&(edges['filter']=='All')&(edges.seed==17)&(edges.tf==tf)&edges.retained_DropOverlap].sort_values(['weight','target'],ascending=[False,True]).head(12)
  e=e.copy();e['display_reason']='All, seed17, top12 original edge weights after original DropOverlap; visualization only';display.append(e)
  ys=np.linspace(.05,.95,len(e));ax.scatter([.12],[.5],s=1600,color=COL[tf],zorder=5);ax.text(.12,.5,tf,ha='center',va='center',color='white',weight='bold')
  for y,(_,q) in zip(ys,e.iterrows()):
   ax.add_patch(FancyArrowPatch((.25,.5),(.91,y),arrowstyle='->',mutation_scale=8,color='#a8b4bc',linewidth=.9,connectionstyle='arc3,rad=0'))
   ax.scatter([.94],[y],s=25,color='#536f80');ax.text(.98,y,q.target,va='center',fontsize=8)
  ax.set(xlim=(-.1,1.65),ylim=(-.02,1.02));ax.axis('off');ax.set_title(f'{sp} | {tf} (All, seed17)',fontsize=12)
fig.suptitle('F06 | Motif-pruned regulon edges: display subset, not causal paths',fontsize=14,y=1.0);fig.tight_layout()
pd.concat(display).to_csv(P/'data/network_display_edges.tsv',sep='\t',index=False)
done(fig,'F06_regulon_network_display','两候选网络：完整边表的固定展示子集',
'源于既有SCENIC1 motif-pruned regulons。展示规则固定为All、seed17、按原边权重降序的12条边（同权按基因名），沿用原DropOverlap排除；只用于版面简化，不是新科学筛选阈值。完整All/BG1、seed17/29/43及排除标记保留于TF_network_all_edges.tsv。箭头表示预测TF-靶基因关系，不表示已证实结合、激活/抑制方向或巨噬细胞到TF的因果通路。','data/TF_network_all_edges.tsv; data/network_display_edges.tsv')
# F07: all frozen modules, every seed/filter, no P-based display selection.
a=rd('S2_TF_associations.tsv');features=['score_OXPHOS','score_MitoISR','score_MitoROS','score_Mitophagy','score_PyroExec','score_ApoExec','score_TAct','score_TReg','score_Inhib','score_Inflam','score_Naive','score_Cycle']
for sp,scope,model in [('Mouse','ActCD4','Group'),('Human','Core','SourceDep')]:
 fig,axs=plt.subplots(2,1,figsize=(13,7),sharex=True)
 for ax,tf in zip(axs,['IKZF2','STAT1']):
  z=a[(a.species==sp)&(a.scope==scope)&(a.tf==tf)&(a.variant=='DropOverlap')&(a.model==model)]
  vals=np.full((6,12),np.nan);txt={}
  for i,(arm,seed) in enumerate(pairs):
   for j,f in enumerate(features):
    q=z[(z['filter']==arm)&(z.seed==seed)&(z.feature==f)]
    assert len(q)<=1
    if len(q):
     q=q.iloc[0];vals[i,j]=q.rho
     sig=flag(q.p,q.FDR)
     txt[i,j]=f'{q.rho:.2f}'+('\n'+sig if sig in ['N','F'] else '')
  im=ax.imshow(np.ma.masked_invalid(vals),cmap='RdBu_r',vmin=-1,vmax=1,aspect='auto')
  ax.set_yticks(range(6),[f'{arm} s{seed}' for arm,seed in pairs]);ax.set_title(tf,loc='left',fontsize=12)
  for (i,j),t in txt.items():ax.text(j,i,t,ha='center',va='center',fontsize=8,color='white' if abs(vals[i,j])>.68 else '#202020')
  ax.set_xticks(range(12),[f.replace('score_','') for f in features],rotation=30,ha='right')
 fig.suptitle(f'F07 | {sp} {scope}: TF - transcript-module association ({model})',fontsize=14)
 fig.subplots_adjust(top=.89,bottom=.16,hspace=.28,left=.09,right=.86)
 cax=fig.add_axes([.89,.19,.014,.65]);fig.colorbar(im,cax=cax,label='Partial rank correlation')
 done(fig,'F07_'+sp.lower()+'_module_associations',sp+'：两候选与冻结模块的关联',
 f'{sp}，范围{scope}，原模型{model}，DropOverlap，按All/BG1和三seed完整展示12个冻结模块。色值为原偏秩相关rho；N=原名义P<0.05，F=原BH FDR<0.05，未重算BH或择优挑模型。'+('小鼠n=6原规则不提供相关P，因此全图仅描述，不作为正式名义显著关联。' if sp=='Mouse' else '人Core各比较沿用原来源/组别/接收深度调整及原完整检验家族。')+' MitoROS是抗氧化应答转录评分，PyroExec/ApoExec不是蛋白裂解或细胞死亡实验证据；模块之间可有原定义重叠。','data/S2_TF_associations.tsv; source_tables/M1_module_definition.tsv')
# F08: audit existing working labels, not a new annotation.
rev=rd('myeloid_working_annotation_audit.tsv').sort_values('cluster')
exp=rd('S3_human_cluster_marker_expression.tsv','source_tables')
markers=['LST1','TYROBP','C1QA','C1QB','FABP4','MARCO','PPARG','FCN1','S100A8','VCAN','SPP1','GPNMB','CD1C','FCER1A','CLEC9A','EPCAM','SFTPC','KRT19','CD3D','TRAC','COL1A1','DCN']
fig,(ax,bar)=plt.subplots(1,2,figsize=(15,8),gridspec_kw={'width_ratios':[6,1.4]},sharey=True)
maxcol=np.log1p(exp.mean_CP10k.max())
for i,c in enumerate(rev.cluster):
 for j,g in enumerate(markers):
  z=exp[(exp.r4==c)&(exp.gene==g)]
  if len(z):
   q=z.iloc[0];ax.scatter(j,i,s=max(q.pct,0)*1.25,c=[np.log1p(q.mean_CP10k)],cmap='viridis',vmin=0,vmax=maxcol,edgecolor='none')
ax.set_xticks(range(len(markers)),markers,rotation=65,ha='right');ax.set_yticks(range(len(rev)),[f'r4:{q.cluster} {q.label} | n={q.cells}' for q in rev.itertuples()]);ax.invert_yaxis();ax.set_xlim(-.7,len(markers)-.3)
ax.grid(axis='y',color='#ededed');ax.set_axisbelow(True);ax.set_title('Frozen marker expression | existing working annotation')
yy=np.arange(len(rev));bar.barh(yy-.16,rev.dominant_subject_fraction*100,height=.3,color='#d9904a',label='Top donor');bar.barh(yy+.16,rev.parent16_fraction*100,height=.3,color='#748da6',label='Parent16')
bar.set_xlim(0,103);bar.set_xlabel('Cells (%)');bar.tick_params(axis='y',labelleft=False);bar.set_title('Composition audit');bar.legend(frameon=False,fontsize=8,loc='lower right')
fig.suptitle('F08 | Myeloid working-label audit: original mapping retained',fontsize=15)
fig.text(.5,.01,'Dot size = detected-cell %. Color = log1p(cluster mean CP10k). All 15 clusters retained in audit; sender scopes follow the prior review.',ha='center',fontsize=9)
fig.subplots_adjust(top=.90,bottom=.18,left=.21,right=.92,wspace=.1)
cax=fig.add_axes([.95,.20,.012,.64]);fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0,maxcol),cmap='viridis'),cax=cax,label='log1p(mean CP10k)')
done(fig,'F08_myeloid_working_label_audit','人髓系工作注释：原映射审阅证据',
'完全沿用2026-10-01原reviewed_annotation.tsv；本轮不重新判群。点大小为已导出簇内检出百分比，颜色为log1p原簇均值CP10k；22基因为固定展示面板，完整marker及程序/来源证据均附表。右侧显示主导供者比例和原parent16比例。r4=2/11分别由CTL068约98.8%与CTL070约97.2%主导，不命名为疾病新亚型。MacAM=0/2/11/13，MonoMac=3/7，MacWide为两者并集；其他簇保留审阅但不进入原发送端模型。工作注释已用于既有模型，不等于用户逐标签最终确认。','data/myeloid_working_annotation_audit.tsv; source_tables/S3_human_cluster_marker_expression.tsv; prior_review_repaired/scenic3_review/results/reviewed_annotation.tsv')
# F09: preserve fixed main line and all stronger-adjustment outcomes.
main=rd('fixed_IL18_IKZF2_all_sensitivities.tsv')
fig,axs=plt.subplots(1,3,figsize=(13,5.3),sharex=True,sharey=True)
specs=[('All','Matched'),('All','NoDom'),('BG1','Matched')]
seedcols={17:'#286b8c',29:'#9a7655',43:'#b4544c'}
labels=[f'{m} s{s}' for m in ['Spearman','SourceDep','BothEpi'] for s in [17,29,43]]
for ax,(arm,ss) in zip(axs,specs):
 z=main[(main['filter']==arm)&(main.sample_set==ss)&(main.sender_filter=='All')]
 for i,(m,seed) in enumerate([(m,s) for m in ['Spearman','SourceDep','BothEpi'] for s in [17,29,43]]):
  q=z[(z.model==m)&(z.seed==seed)].iloc[0]
  if pd.notna(q.loo_min):ax.plot([q.loo_min,q.loo_max],[i,i],color=seedcols[seed],alpha=.4,lw=2)
  if pd.isna(q.p):ax.scatter(q.rho,i,marker='x',color=seedcols[seed],s=42,lw=1.7)
  else:ax.scatter(q.rho,i,marker='o',edgecolor=seedcols[seed],facecolor=seedcols[seed] if q.p<.05 else 'white',s=45,lw=1.4)
 ax.set_title(f'{arm} / {ss} | n={int(z.n.min())}',fontsize=12);ax.axvline(0,color='#888',ls='--');ax.set_xlim(-.05,1.04);ax.set_yticks(range(9),labels);ax.invert_yaxis();ax.set_xlabel('rho; line = leave-one-out range');ax.grid(axis='x',color='#eee')
fig.suptitle('F09 | Fixed patient line: IL18 (MacAM) - IKZF2 (CD4-IL7R)',fontsize=14,y=1.01)
fig.text(.5,-.015,'Filled: nominal P < 0.05. Open: not nominal. x: no valid P by original rule. No model here passes FDR < 0.05.',ha='center',fontsize=9)
fig.tight_layout()
done(fig,'F09_patient_fixed_axis_robustness','患者固定主线：IL18–IKZF2完整敏感性',
'仅整理原先已汇报的固定主线，不按最小P重新挑选。MacAM发送IL18，CD4-IL7R接收IKZF2 DropOverlap活性。完整保留Spearman、SourceDep、BothEpi及三seed；图示All/Matched、All/NoDom、BG1/Matched，完整NoParent16和其余组合见表。实心=名义P<0.05，空心=未达名义，叉号=原规则不给有效P。线为LOO相关范围，不是CI。All/Matched中SourceDep名义通过，BothEpi未通过且df由5降3；不能仅归因为功效，混杂敏感仍可能。BG1 n=9无相关P，BothEpi df=1且LOO仅少量可估计，近1的rho不得解读为稳健证据。所有FDR未通过。','data/fixed_IL18_IKZF2_all_sensitivities.tsv; prior_review_repaired/scenic3_review/results/reviewed_sender_TF_models.tsv')
(P/'figure_manifest.json').write_text(json.dumps(caps,ensure_ascii=False,indent=2),encoding='utf-8')
(P/'图注与缩写说明.txt').write_text('\n\n'.join(x['figure']+' | '+x['title']+'\n'+x['caption_zh']+'\n来源：'+x['source'] for x in caps)+'\n\n通用缩写：All=原主过滤臂；BG1=固定背景比例≤1%的敏感性；s17/s29/s43=随机种子，不是独立生物重复；CP10k=按原全库总数归一化到每万计数；AUCell=regulon活性评分；LOO=逐一留出受试者/文库；N=名义P<0.05，F=原BH FDR<0.05；NoDom=原定义移除配对集合中CTL068/CTL070；NoParent16=原父簇16排除敏感性；CWP=原人肺病例标签，不改写为新诊断。SourceDep=组别/来源/接收深度调整；BothEpi=另加发送深度和连续上皮负担。\n',encoding='utf-8')
print('FIGURES_CREATED='+str(len(caps)))
