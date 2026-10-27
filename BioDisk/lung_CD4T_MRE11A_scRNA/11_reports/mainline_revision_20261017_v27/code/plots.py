from pathlib import Path
import pandas as pd,numpy as np,json,shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,Patch,Circle,FancyBboxPatch
from matplotlib import font_manager
O=Path(__file__).parent;H=O.parent;F=O/'figures';T=O/'plot_data'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
a=pd.read_csv(H/'lung_chord_clarity_20261017_v23/plot_data/F20.tsv',sep='\t');a['interpretation']=np.where(a.ligand.str.startswith('H2-DM'),'antigen_presentation_annotation',a.route);a.to_csv(T/'F20.tsv',sep='\t',index=False)
rows=a[['ligand','receptor']].drop_duplicates();rows['dm']=rows.ligand.str.startswith('H2-DM');rows=rows.sort_values(['dm','ligand','receptor']).reset_index(drop=True)
fig,ax=plt.subplots(figsize=(7,8));fig.subplots_adjust(left=.36,right=.97,bottom=.25,top=.92)
for i,r in rows.iterrows():
 for j,g in enumerate(['IL1R2','TTC39C','SKIL']):
  z=a[(a.ligand==r.ligand)&(a.receptor==r.receptor)&(a.gene==g)];color='#eef0f2';lab='';hatch=None
  if g=='TTC39C':color='#fff9e6';lab='NA';hatch='//'
  elif len(z):
   rr=z.iloc[0];color='#ad94bd' if r.dm else '#d9994e' if rr.route=='direct_receptor' else '#397fa5';lab=('Prior '+str(int(rr.target_rank))) if pd.notna(rr.target_rank) else 'Decoy'
  ax.add_patch(Rectangle((j-.48,i-.46),.96,.92,facecolor=color,edgecolor='#c5cad0' if hatch else 'white',hatch=hatch,linewidth=1.5));ax.text(j,i,lab,ha='center',va='center',fontsize=10,color='white' if color=='#397fa5' else '#263238')
ax.axhline(5.5,color='#745386',lw=2);ax.set_xlim(-.5,2.5);ax.set_ylim(7.5,-.5);ax.set_xticks(range(3),['IL1R2\nTreg','TTC39C\nActCD4','SKIL\nTreg']);ax.set_yticks(range(8),[r.ligand+'\n('+r.receptor.replace(';',' + ')+')' for r in rows.itertuples()]);ax.tick_params(length=0,pad=8)
for sp in ax.spines.values():sp.set_visible(False)
ax.set_title('Candidate connections by evidence type',fontsize=12,pad=17)
fig.legend(handles=[Patch(color='#397fa5',label='NicheNet target prior (rank)'),Patch(color='#d9994e',label='Direct decoy binding'),Patch(color='#ad94bd',label='Antigen-presentation annotation only'),Patch(facecolor='#fff9e6',hatch='//',label='NA: candidate not included in this analysis'),Patch(color='#eef0f2',label='No supporting route in displayed set')],loc='lower center',bbox_to_anchor=(.51,.015),frameon=False,fontsize=9)
for ext in ['png','pdf']:fig.savefig(F/('F20.'+ext),dpi=300,bbox_inches='tight',pad_inches=.15)
plt.close(fig)
# Redraw screening figures using identical counts and revised labels only.
fp=font_manager.FontProperties(fname=r'C:\Windows\Fonts\ARIALUNI.TTF')
def text(ax,x,y,s,size=11):ax.text(x,y,s,fontproperties=fp,fontsize=size,color='#193e53',ha='center',va='center',linespacing=1.5)
def box(ax,x,y,w,h,s,color='#e6eff4',size=11):ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.03',fc=color,ec='#9fb7c4'));text(ax,x+w/2,y+h/2,s,size)
def arrow(ax,x,y,x2,y2):ax.annotate('',xy=(x2,y2),xytext=(x,y),arrowprops=dict(arrowstyle='-|>',color='#6a8798'))
def init():
 fig,ax=plt.subplots(figsize=(7,8));fig.subplots_adjust(left=.03,right=.97,bottom=.03,top=.97);ax.set(xlim=(0,10),ylim=(0,12));ax.axis('off');return fig,ax
def save(fig,n):
 for ext in ['png','pdf']:fig.savefig(F/(n+'.'+ext),dpi=300,bbox_inches='tight',pad_inches=.15)
 plt.close(fig)
fig,ax=init();text(ax,5,11.55,'候选基因筛选策略',16)
for xy,c in [((3.25,7.6),'#65a6ca'),((6.75,7.6),'#e7b45a'),((5,5.15),'#73b7a0')]:ax.add_patch(Circle(xy,2.6,fc=c,ec=c,alpha=.32))
text(ax,2.3,8.5,'差异证据\n目标效应P < 0.05\n保留状态与对比',10.5);text(ax,7.7,8.5,'目标功能\n内源性凋亡/线粒体损伤\n或焦亡注释与文献',10);text(ax,5,4.0,'巨噬来源信号依据\n同一接收状态的受体\n及配体-靶先验',10.5);text(ax,5,6.9,'三项交集\n55个基因',12)
arrow(ax,5,2.5,5,2.0);box(ax,.6,.8,8.8,1.15,'55个交集 + 42个补充 − 7个重叠 = 90个评价基因\n创新性、证据完整度与实验可行性用于分层',size=10.5);text(ax,5,.25,'圆面积不代表数量；完整阈值与去向见附录A。',9);save(fig,'S_C')
fig,ax=init();text(ax,5,11.55,'差异、功能与巨噬来源信号的交集',14)
items=[('完整受检集合\n9,091个基因；52,448条记录','#edf1f4'),('目标效应P < 0.05\n2,765个基因；6,126条记录','#e6eff4'),('差异 ∩ 目标功能\n111个基因；274条记录','#edf2de'),('差异 ∩ 目标功能 ∩ 巨噬来源信号\n55个基因；155条记录','#dbeaf2'),('与42个补充候选合并，重叠7个\n55 + 42 − 7 = 90个评价基因','#f4ead7')]
for i,(label,c) in enumerate(items):
 y=9.5-i*1.85;box(ax,.6,y,8.8,1.1,label,c)
 if i<4:arrow(ax,5,y-.06,5,y-.68)
text(ax,5,.6,'按基因-状态-对比匹配后去重；补充候选为并集扩展。',10);save(fig,'F21')
fig,ax=init();text(ax,5,11.55,'90个基因的证据分层与实验推荐',14)
labels=['名义差异：89个基因 / 344条记录\n1个未达到名义差异条件','目标功能：57个基因 / 172条记录\n32个功能证据待补','同状态信号：56个基因 / 168条记录\n1个信号证据待补（IL10）','细胞来源：54个基因 / 164条记录\n2个来源证据待补（GSDMD、S100A8）','作用角色：52个基因 / 160条记录\n2个作用定位待补（GZMB、TNFSF10）']
for i,label in enumerate(labels):
 y=9.4-i*1.55;box(ax,.6,y,8.8,1.15,label,size=10.5)
 if i<4:arrow(ax,5,y-.04,5,y-.36)
box(ax,.6,.55,8.8,1.6,'分层：52个计算支持 + 37个证据待补 + 1个未达条件\n9个重点基因来自全部90基因：5个计算支持、4个待补\n结合功能关联与文献，推荐IL1R2、TTC39C、SKIL',color='#e1efdf',size=10.5);save(fig,'F22')
print('FIGURES_DONE')