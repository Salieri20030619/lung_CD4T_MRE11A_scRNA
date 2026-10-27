from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch
O=Path(__file__).parent;fp=font_manager.FontProperties(fname=r'C:\Windows\Fonts\ARIALUNI.TTF')
fig,ax=plt.subplots(figsize=(7.4,8.8));fig.subplots_adjust(left=.02,right=.98,bottom=.03,top=.97);ax.set(xlim=(0,10),ylim=(0,12));ax.axis('off')
def t(x,y,s,n=11):ax.text(x,y,s,fontproperties=fp,fontsize=n,ha='center',va='center',color='#193e53',linespacing=1.45)
def box(x,y,w,h,s,c='#e7eff4',n=11):ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.04',fc=c,ec='#9bb3c0'));t(x+w/2,y+h/2,s,n)
def ar(x,y,xx,yy):ax.annotate('',xy=(xx,yy),xytext=(x,y),arrowprops=dict(arrowstyle='-|>',color='#5e7b8c',lw=1.2))
t(5,11.6,'资格评价与定向提名分别记录',15);box(.7,10.1,8.6,.9,'55个交集 + 42个补充 − 7个重叠 = 90个评价基因',n=11)
ar(3,10.05,2.5,9.5);ar(7,10.05,7.5,9.5)
box(.4,6.7,4.2,2.7,'资格分层\n52个计算支持\n37个证据待补\n1个未达名义差异条件\n完整成员均保留',c='#e6efde',n=10.5)
box(5.4,6.7,4.2,2.7,'定向提名9个深入评价对象\n5个功能机制候选\nIL1R2 / DDIT4 / SKIL / FOXP1 / GSN\n4个补充探索对象\nAHNAK / TTC39C / MDFIC / ST6GALNAC3',c='#f5ead5',n=8.8)
t(5,6.2,'9个对象含5个计算支持与4个待补；提名与资格分层分别记录',9.8)
ar(7.5,6.65,7.5,5.6);box(.7,3.6,8.6,1.9,'在已选9个对象中比较后续结果\n疾病表达与样本定位 → 表型/损伤程序关联 → 调整与外部方向\n结合功能文献及目标细胞可测性，记录支持与主要缺口',n=10.5)
ar(5,3.55,5,2.9);box(.7,.9,8.6,1.9,'推荐三条预实验路线\nIL1R2：Treg与IL-1环境响应\nTTC39C：存活/凋亡与线粒体状态探索\nSKIL：Treg状态调控',c='#dceceb',n=11)
t(5,.35,'其他6个深入评价对象保留为备选或参照；推荐依据为证据互补性与可操作性。',9)
for ext in ['png','pdf']:fig.savefig(O/'figures'/('F22.'+ext),dpi=300,bbox_inches='tight',pad_inches=.15)
plt.close(fig)
p=O/'figure_manifest.json';m=json.loads(p.read_text(encoding='utf-8'))
for x in m:
 if x['id']=='F22':x['file']=str(O/'figures/F22.png');x['title']='90基因资格分层与9对象提名分别通向综合评价';x['caption']='90个基因的资格分层为52/37/1；定向深入评价由5个功能机制候选与4个补充探索对象组成，包含5个计算支持和4个证据待补。9个对象在追加关联分析前确定；随后综合结果和实验可行性推荐3条路线，其他6个保留。'
p.write_text(json.dumps(m,ensure_ascii=False,indent=2),encoding='utf-8');print('SELECTION_FLOW_DONE')