from pathlib import Path
import json,shutil,numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch,Circle
O=Path(__file__).parent;H=O.parent;B=H/'lung_client_rebuild_20261017_v20';F=O/'figures';F.mkdir(exist_ok=True);T=O/'plot_data';T.mkdir(exist_ok=True);fp=font_manager.FontProperties(fname=r'C:\Windows\Fonts\ARIALUNI.TTF')
a=pd.read_csv(H/'lung_three_cohort_20261017_v16/data/tables/sample_inventory.tsv',sep='\t');core=a[(a.scope=='Core')&a.dataset.isin(['GSE314943','GSE183682'])].copy();assert len(core)==10 and core.groupby('group').size().to_dict()=={'Saline':5,'Silica':5};core.to_csv(T/'joint10_sample_roster.tsv',sep='\t',index=False)
def txt(ax,x,y,s,size=11,color='#193f54',ha='center'):return ax.text(x,y,s,fontproperties=fp,fontsize=size,color=color,ha=ha,va='center',linespacing=1.5)
def box(ax,x,y,w,h,s,fill='#e8f0f4',size=11):ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.03,rounding_size=.1',fc=fill,ec='#93aebc',lw=.8));txt(ax,x+w/2,y+h/2,s,size)
def arrow(ax,x,y,x2,y2):ax.annotate('',xy=(x2,y2),xytext=(x,y),arrowprops={'arrowstyle':'-|>','lw':1.3,'color':'#688797'})
def save(fig,id):
 for ext in ['png','pdf']:fig.savefig(F/(id+'.'+ext),dpi=300,bbox_inches='tight',pad_inches=.15)
 plt.close(fig)
fig,ax=plt.subplots(figsize=(7,8));fig.subplots_adjust(left=.02,right=.98,bottom=.03,top=.97);ax.set(xlim=(0,10),ylim=(0,12));ax.axis('off');txt(ax,5,11.65,'图谱范围与联合统计范围',15)
box(ax,.35,8.55,4.35,2.4,'主队列 GSE314943\nD140：3个对照 + 3个矽尘样本\n6个样本；Core CD4共2,311细胞\n用于全肺/CD4图谱与主队列分析',fill='#e3eff5',size=10.3)
box(ax,5.3,8.55,4.35,2.4,'补充队列 GSE183682\nD7：1个对照 + 1个矽尘样本\nD56：1个对照 + 1个矽尘样本\n4个样本；Core CD4共1,413细胞',fill='#f4ead7',size=10.3)
arrow(ax,2.5,8.5,4,6.8);arrow(ax,7.5,8.5,6,6.8)
box(ax,1.0,4.85,8.0,1.9,'联合疾病分析：6 + 4 = 10个样本\n5个对照 + 5个矽尘样本\n模型：队列-时间层 + 疾病组别',fill='#e1efdf',size=12)
txt(ax,5,4.25,'三个层：GSE314943-D140；GSE183682-D7；GSE183682-D56',9.6)
box(ax,.55,2.6,8.9,1.15,'图谱保留各队列的独立坐标\n联合模型按对应细胞状态汇总样本级表达',fill='#f0f2f5',size=10.5)
box(ax,.55,.7,8.9,1.3,'时间队列 GSE274372：D0/D7/D28/D56各3个样本\n共12个样本，用于时间背景分析，不计入上述10个样本\n联合模型包含主队列，不能作为其独立复现',fill='#f4f0e6',size=10.5)
save(fig,'S_A')
# Actual existing supplementary-cohort coordinates; validate states against sample inventory.
z=pd.read_csv(O/'tables/GSE183682_CD4_coordinates.tsv.gz',sep='\t');assert len(z)==1413;counts=z.groupby(['sample','state']).size().rename('cells').reset_index();expect=a[(a.dataset=='GSE183682')&a.scope.isin(['Naive','ActCD4','Treg'])][['sample','scope','cells']].rename(columns={'scope':'state'});pd.testing.assert_frame_equal(counts.sort_values(['sample','state']).reset_index(drop=True),expect.sort_values(['sample','state']).reset_index(drop=True))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42});fig,axs=plt.subplots(2,1,figsize=(7,8.2));fig.subplots_adjust(left=.14,right=.74,bottom=.11,top=.94,hspace=.43);colors={'Naive':'#4c9e9c','ActCD4':'#dc8451','Treg':'#8b74aa'}
for state,c in colors.items():
 zz=z[z.state==state];axs[0].scatter(zz.UMAP1,zz.UMAP2,s=6,c=c,alpha=.65,label=f'{state} ({len(zz)})',rasterized=True)
axs[0].set(title='GSE183682 | CD4 states (n=4 samples)',xlabel='UMAP 1',ylabel='UMAP 2');axs[0].legend(loc='upper left',bbox_to_anchor=(1.01,1),frameon=False,fontsize=9)
w=counts.pivot(index='sample',columns='state',values='cells').reindex(columns=list(colors));bottom=np.zeros(len(w));labs=[]
for state,c in colors.items():axs[1].bar(range(len(w)),w[state],bottom=bottom,color=c,label=state);bottom+=w[state].to_numpy()
for s in w.index:
 row=core[core['sample']==s].iloc[0];labs.append(s+'\nD'+str(int(row.day))+' / '+row.group)
axs[1].set_xticks(range(len(w)),labs,rotation=25,ha='right',fontsize=8);axs[1].set_ylabel('CD4 cells');axs[1].set_title('Four samples included in the joint model')
for i,n in enumerate(bottom):axs[1].text(i,n+9,str(int(n)),ha='center',fontsize=9)
axs[1].set_ylim(0,bottom.max()*1.18);save(fig,'S_B');z.to_csv(T/'S_B.tsv.gz',sep='\t',index=False);counts.to_csv(T/'S_B_sample_counts.tsv',sep='\t',index=False)
# Conceptual overlap, not area-proportional gene-set Venn counts.
fig,ax=plt.subplots(figsize=(7,8));fig.subplots_adjust(left=.02,right=.98,bottom=.03,top=.97);ax.set(xlim=(0,10),ylim=(0,12));ax.axis('off');txt(ax,5,11.6,'从差异基因出发的候选证据框架',15)
for xy,c in [((3.25,7.6),'#65a6ca'),((6.75,7.6),'#e7b45a'),((5,5.15),'#73b7a0')]:ax.add_patch(Circle(xy,2.6,facecolor=c,edgecolor=c,alpha=.32,lw=1.5))
txt(ax,2.35,8.35,'差异证据\n目标效应P < 0.05\n保留状态与对比',11)
txt(ax,7.65,8.35,'目标功能\n内源性凋亡/线粒体损伤\n或焦亡的注释与文献',10.4)
txt(ax,5,4.0,'巨噬来源入口\n同状态配体-受体及\nNicheNet靶先验',11)
txt(ax,5,6.9,'三项证据交集\n55个基因',12)
arrow(ax,5,2.5,5,2.05);box(ax,.6,.95,8.8,1.0,'结合补充候选形成90基因评价集合\n创新性、证据完整度和实验可行性用于后续分层',fill='#e9edf4',size=10.6)
txt(ax,5,.38,'圆面积不代表基因数量；详细计数与候选去向见后文筛选流程。',9)
save(fig,'S_C')
# Keep figure5's separately validated label rail; replace only its cohort title.
source=(H/'lung_screening_display_20261017_v24/plots.py').read_text(encoding='utf-8');part=source[:source.index("fp=font_manager.FontProperties")];part=part.replace("scope+' | joint disease model (n=10)'","scope+' | GSE314943 + GSE183682\\nStratified disease model (6 + 4 samples)'");ns={'__file__':str(O/'plots.py')};exec(part,ns)
# Update figure6's legend only, using its already separated text column.
source=(H/'lung_client_layout_20261017_v22/portrait_plots.py').read_text(encoding='utf-8');prefix=source[:source.index('# Atlas legend')];ns={'__file__':str(H/'lung_client_layout_20261017_v22/portrait_plots.py')};exec(prefix,ns);ns['O']=O;ns['F']=F;ns['T']=T
block=source[source.index("a=rd('F06')"):source.index('# GO and KEGG')];block=block.replace("'Original six'","'GSE314943 (n=6)'").replace("'Joint ten'","'GSE314943 + GSE183682 (n=10)'").replace("bbox_to_anchor=(.55,.99),ncol=2","bbox_to_anchor=(.55,.995),ncol=1");exec(block,ns)
(O/'cohort_QA.json').write_text(json.dumps({'primary_samples':6,'supplemental_samples':4,'joint_samples':10,'joint_group_counts':{'Saline':5,'Silica':5},'joint_strata':3,'primary_CD4_cells':2311,'supplemental_CD4_cells':1413,'supplemental_embedding':'existing X_umap; no joint embedding computed','supplemental_counts_match_original':True,'temporal_samples_excluded_from_joint':12,'venn_is_conceptual':True,'statistics_changed':False},ensure_ascii=False,indent=2),encoding='utf-8');print('PLOTS_READY')