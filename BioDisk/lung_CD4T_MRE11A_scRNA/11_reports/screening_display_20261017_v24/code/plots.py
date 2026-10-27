from pathlib import Path
import pandas as pd,numpy as np,json,shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch,FancyBboxPatch
from matplotlib import font_manager
O=Path(__file__).parent;B=O.parent/'lung_client_rebuild_20261017_v20';F=O/'figures';F.mkdir(exist_ok=True);T=O/'plot_data';T.mkdir(exist_ok=True);plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
a=pd.read_csv(B/'plot_data/F05.tsv',sep='\t');fig,axs=plt.subplots(2,1,figsize=(7,8.6));fig.subplots_adjust(left=.13,right=.72,bottom=.07,top=.90,hspace=.43);focus={'Treg':['IL1R2','SKIL','FOXP1'],'ActCD4':['TTC39C','MDFIC','DDIT4','GSN']};texts=[];positions=[]
fig.legend(handles=[Patch(color='#bc4d42',label='FDR < 0.05'),Patch(color='#d6a342',label='P < 0.05 only'),Patch(color='#b9c1c7',label='Other')],loc='upper center',bbox_to_anchor=(.5,.995),ncol=3,frameon=False,fontsize=9)
for ax,scope in zip(axs,['Treg','ActCD4']):
 z=a[a.scope==scope];ax.scatter(z.logFC,z.neglogP,c=np.where(z.FDR<.05,'#bc4d42',np.where(z.PValue<.05,'#d6a342','#b9c1c7')),s=7,alpha=.6,rasterized=True);ax.axhline(-np.log10(.05),c='grey',ls='--',lw=.6);ax.axvline(0,c='grey',lw=.5)
 zz=z[z.gene.isin(focus[scope])].sort_values('neglogP',ascending=False)
 for y,r in zip(np.linspace(.87,.18,len(zz)),zz.itertuples()):
  ax.scatter([r.logFC],[r.neglogP],s=24,c='#173e58',zorder=4);ax.annotate('',xy=(r.logFC,r.neglogP),xycoords='data',xytext=(1.035,y),textcoords='axes fraction',arrowprops={'arrowstyle':'-','color':'#526d7e','lw':.7},annotation_clip=False);text=ax.text(1.055,y,r.gene,transform=ax.transAxes,va='center',ha='left',fontsize=10.5,color='#173e58',clip_on=False,bbox=dict(boxstyle='round,pad=.18',fc='white',ec='none'));texts.append(text);positions.append({'scope':scope,'gene':r.gene,'logFC':r.logFC,'neglogP':r.neglogP,'label_x_axes':1.055,'label_y_axes':y})
 ax.set(title=scope+' | joint disease model (n=10)',xlabel='log2 fold change (Silica / Saline)',ylabel='-log10 raw P')
fig.canvas.draw();renderer=fig.canvas.get_renderer();boxes=[t.get_window_extent(renderer).expanded(1.05,1.15) for t in texts]
assert all(not boxes[i].overlaps(boxes[j]) for i in range(len(boxes)) for j in range(i+1,len(boxes)))
assert all(fig.bbox.contains(b.x0,b.y0) and fig.bbox.contains(b.x1,b.y1) for b in boxes)
assert all(not b.overlaps(ax.get_window_extent(renderer)) for b in boxes for ax in axs)
for ext in ['png','pdf']:fig.savefig(F/('F05.'+ext),dpi=300,bbox_inches='tight',pad_inches=.15)
plt.close(fig);shutil.copy2(B/'plot_data/F05.tsv',T/'F05.tsv');pd.DataFrame(positions).to_csv(T/'F05_label_positions.tsv',sep='\t',index=False);(O/'F05_label_QA.json').write_text(json.dumps({'labels':len(texts),'pairwise_bbox_overlap':False,'outside_plot_axes':True,'within_figure':True,'data_unchanged':True},indent=2),encoding='utf-8')
fp=font_manager.FontProperties(fname=r'C:\Windows\Fonts\ARIALUNI.TTF')
def text(ax,x,y,s,size=11,co='#193e53',ha='center'):return ax.text(x,y,s,fontproperties=fp,fontsize=size,color=co,ha=ha,va='center',linespacing=1.45,bbox=dict(facecolor='white',edgecolor='none',pad=1.2) if (x==5 and y in [9.96,8.42,6.90,1.62,5.25]) else None)
def box(ax,x,y,w,h,s,fill='#eaf1f5',size=11):
 ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.03,rounding_size=.12',facecolor=fill,edgecolor='#91abb8',lw=.8));text(ax,x+w/2,y+h/2,s,size)
def arrow(ax,x,y,x2,y2,col='#6b8796',style='-'):ax.annotate('',xy=(x2,y2),xytext=(x,y),arrowprops=dict(arrowstyle='-|>',lw=1.2,color=col,linestyle=style))
fig,ax=plt.subplots(figsize=(7,8.6));fig.subplots_adjust(left=.02,right=.98,bottom=.03,top=.97);ax.set(xlim=(0,10),ylim=(0,12));ax.axis('off')
text(ax,5,11.6,'差异、功能与巨噬入口的交集构建',14)
box(ax,.7,10.3,8.6,.8,'完整主要受检集合\n9,091个基因；52,448条基因-状态-对比记录')
arrow(ax,5,10.25,5,9.7);text(ax,5,9.96,'目标效应原始P < 0.05',9.8)
box(ax,.7,8.8,8.6,.85,'差异证据集合\n2,765个基因；6,126条记录')
arrow(ax,5,8.75,5,8.13);text(ax,5,8.42,'与目标损伤/死亡功能集合相交',9.8)
box(ax,.7,7.25,8.6,.83,'差异 ∩ 功能\n111个基因；274条记录')
arrow(ax,5,7.20,5,6.58);text(ax,5,6.90,'再匹配同一接收状态的巨噬来源入口',9.8)
box(ax,.7,5.7,8.6,.83,'差异 ∩ 功能 ∩ 巨噬入口\n55个基因；155条记录',fill='#daeaf1')
text(ax,5,5.25,'按状态和对比匹配后，再对基因去重',9.5)
box(ax,6.2,3.8,3.1,.9,'补充候选集合\n42个基因',fill='#f5eddc',size=10.5)
arrow(ax,4.2,5.65,4.2,2.45);arrow(ax,7.75,3.75,5.85,2.45)
text(ax,3.2,3.9,'集合合并',11);text(ax,3.2,3.4,'两集合重叠7个',10.5)
box(ax,1.2,1.45,7.6,.95,'综合评价集合：90个基因\n55 + 42 - 7 = 90',fill='#dcefdc',size=12)
text(ax,5,.83,'补充集合按研究问题纳入，名单及每项依据见附表。',9.5)
text(ax,5,.39,'该合并为并集扩展；不产生新的联合P值或集合FDR。',9.5)
for ext in ['png','pdf']:fig.savefig(F/('F21.'+ext),dpi=300,bbox_inches='tight',pad_inches=.12)
plt.close(fig)
fig,ax=plt.subplots(figsize=(7,9.2));fig.subplots_adjust(left=.02,right=.98,bottom=.03,top=.97);ax.set(xlim=(0,10),ylim=(0,13));ax.axis('off');text(ax,5,12.65,'候选资格分层与重点实验方向选择',14)
box(ax,.4,11.45,5.7,.7,'90个基因进入综合评价',fill='#dcefdc')
steps=[(10.25,'名义差异：89个基因 / 344条记录','1个未达到当前\n名义差异条件','#f3e3e3'),(9.0,'目标功能：57个基因 / 172条记录','32个功能证据待补','#fff1d9'),(7.75,'同状态入口：56个基因 / 168条记录','1个入口证据待补\nIL10','#fff1d9'),(6.5,'细胞来源：54个基因 / 164条记录','2个来源证据待补\nGSDMD、S100A8','#fff1d9'),(5.25,'作用角色：52个基因 / 160条记录','2个作用定位待补\nGZMB、TNFSF10','#fff1d9')]
prev=11.45
for y,label,side,fill in steps:
 arrow(ax,3.25,prev-.04,3.25,y+.75);box(ax,.4,y,5.7,.7,label,size=10.5);arrow(ax,6.15,y+.35,6.65,y+.35);box(ax,6.75,y-.03,2.85,.76,side,fill=fill,size=9.5);prev=y
box(ax,.4,3.9,9.2,.75,'资格分层：52个计算支持 + 37个证据待补 + 1个未达名义差异',fill='#dcefdc',size=10.2)
text(ax,5,3.35,'37个待补 = 32 + 1 + 2 + 2；待补不等于没有研究价值',9.5)
box(ax,.4,1.92,9.2,.98,'重点分析从90个候选中选择\n5个计算支持候选 + 4个证据待补候选 = 9个基因 / 10个状态组合',fill='#e6edf8',size=10.5)
arrow(ax,5,1.88,5,1.35);text(ax,5,1.62,'结合程序关联、功能文献、创新空间及实验可行性',9.3)
box(ax,.6,.35,8.8,.9,'优先预实验：IL1R2；探索备选：TTC39C；状态备选：SKIL\n创新性用于分层，不作为统一二元淘汰门槛',fill='#f5eddc',size=10.5)
for ext in ['png','pdf']:fig.savefig(F/('F22.'+ext),dpi=300,bbox_inches='tight',pad_inches=.12)
plt.close(fig)
print('PLOTS_DONE')