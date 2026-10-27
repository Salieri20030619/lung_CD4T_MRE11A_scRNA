from pathlib import Path
import json,shutil,datetime
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,Image,PageBreak
from reportlab.lib.pagesizes import A4
D=Path('stepwise_diagnostic');R=D/'results';F=D/'figures';F.mkdir(exist_ok=True);O=D/'report';O.mkdir(exist_ok=True)
def rd(name):
 d=pd.read_csv(R/name,sep='\t');
 if 'covariates' in d:d['covariates']=d.covariates.fillna('None')
 return d
b=rd('baseline_models.tsv');r=rd('receiver_stepwise_models.tsv');s=rd('sender_stepwise_models.tsv');trans=rd('filter_transition_subset_removals.tsv')
sf=s[(s.scope=='CD4-IL7R')&(s['filter']=='All')&(s.variant=='Matched')&(s.sample_policy=='fixed_complete_cases')&(s.analysis=='SENDER_TF_STEPS')]
palette=['#107D92','#CF7135','#6556AA'];plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
fig,axes=plt.subplots(2,2,figsize=(9.4,5),sharex=True)
for j,tf in enumerate(['IKZF2','STAT1']):
 for k,seed in enumerate([17,29,43]):
  z=sf[(sf.tf==tf)&(sf.seed==seed)&(sf.covariates=='None')].sort_values('threshold')
  axes[0,j].plot(range(4),z.rho,'o-',color=palette[k],label=f'seed {seed}');axes[1,j].plot(range(4),z.p,'o-',color=palette[k])
 axes[0,j].set_title(tf);axes[0,j].set_ylim(-1,1);axes[0,j].axhline(0,color='#cccccc',lw=.7);axes[1,j].set_yscale('log');axes[1,j].axhline(.05,color='#a33',ls='--',lw=1);axes[1,j].set_xticks(range(4),['>=1\nn=17','>=5\nn=15','>=10\nn=13','>=20\nn=11']);axes[1,j].set_xlabel('Minimum cells at both ends')
axes[0,0].set_ylabel('Spearman rho');axes[1,0].set_ylabel('Nominal P (log scale)');axes[0,0].legend(ncol=3,fontsize=8,loc='lower left');fig.tight_layout();fig.savefig(F/'filter_thresholds.png');plt.close(fig)
covs=['None','G','G+S','G+S+R','G+S+R+D','G+S+R+D+E']
fig,axes=plt.subplots(2,2,figsize=(9.4,5),sharex=True)
for j,thr in enumerate([1,20]):
 for k,seed in enumerate([17,29,43]):
  z=sf[(sf.tf=='IKZF2')&(sf.seed==seed)&(sf.threshold==thr)].set_index('covariates').loc[covs]
  axes[0,j].plot(range(6),z.rho,'o-',color=palette[k],label=f'seed {seed}');axes[1,j].plot(range(6),z.p,'o-',color=palette[k])
 axes[0,j].set_title('Fixed n='+str(17 if thr==1 else 11));axes[0,j].set_ylim(0,1);axes[1,j].set_yscale('log');axes[1,j].axhline(.05,color='#a33',ls='--',lw=1);axes[1,j].set_xticks(range(6),['0','G','GS','GSR','GSRD','GSRDE']);axes[1,j].set_xlabel('Added covariates')
axes[0,0].set_ylabel('Partial rank rho');axes[1,0].set_ylabel('Nominal P (log scale)');axes[0,0].legend(ncol=3,fontsize=8);fig.tight_layout();fig.savefig(F/'fixed_roster_adjustment.png');plt.close(fig)
mods=['OXPHOS','MitoISR','MitoROS','Mitophagy','PyroExec','ApoExec'];stages=['None','G','G+S','G+S+R']
fig,axes=plt.subplots(1,2,figsize=(9.4,5))
for ax,tf in zip(axes,['IKZF2','STAT1']):
 z=r[(r.species=='Human')&(r.scope=='CD4-IL7R')&(r['filter']=='All')&(r.tf==tf)&(r.evidence=='TF_MODULE')&(r.threshold==1)&(r.sample_policy=='fixed_complete_cases')]
 vals=np.array([[z[(z.feature==m)&(z.covariates==c)&(z.seed==seed)].iloc[0].estimate for c in stages for seed in [17,29,43]] for m in mods]);ps=np.array([[z[(z.feature==m)&(z.covariates==c)&(z.seed==seed)].iloc[0].p for c in stages for seed in [17,29,43]] for m in mods])
 im=ax.imshow(vals,vmin=-1,vmax=1,cmap='RdBu_r',aspect='auto');ax.set_yticks(range(6),mods);ax.set_xticks([1,4,7,10],['0','G','GS','GSR']);ax.set_title(tf+' / n=17');ax.set_xlabel('Within each block: seeds 17, 29, 43')
 for y in range(6):
  for x in range(12):
   if ps[y,x]<.05:ax.text(x,y,'*',ha='center',va='center',fontsize=12,color='white' if abs(vals[y,x])>.6 else 'black')
 for v in [2.5,5.5,8.5]:ax.axvline(v,color='white',lw=2)
fig.subplots_adjust(wspace=.3,left=.12,right=.89,bottom=.16,top=.9);cax=fig.add_axes([.92,.2,.02,.6]);fig.colorbar(im,cax=cax,label='Partial rank rho');fig.savefig(F/'same_scope_module_associations.png');plt.close(fig)
pdfmetrics.registerFont(TTFont('CN','C:/Windows/Fonts/msyh.ttc',subfontIndex=0))
W,H=A4;navy=colors.HexColor('#123D50');teal=colors.HexColor('#107D92');gray=colors.HexColor('#53616A')
styles={
 'title':ParagraphStyle('title',fontName='CN',fontSize=20,leading=28,textColor=navy,spaceAfter=14),
 'head':ParagraphStyle('head',fontName='CN',fontSize=14,leading=20,textColor=teal,spaceAfter=10),
 'body':ParagraphStyle('body',fontName='CN',fontSize=10,leading=16,spaceAfter=9),
 'small':ParagraphStyle('small',fontName='CN',fontSize=8,leading=12,textColor=gray,spaceAfter=7),
 'cell':ParagraphStyle('cell',fontName='CN',fontSize=8,leading=11),
 'th':ParagraphStyle('th',fontName='CN',fontSize=8,leading=11,textColor=colors.white)}
story=[]
def p(txt,style='body'):story.append(Paragraph(txt,styles[style]))
def title(n,t):p(f'{n:02d} / '+t,'title')
def table(rows,widths=None):
 data=[[Paragraph(str(x),styles['th' if i==0 else 'cell']) for x in row] for i,row in enumerate(rows)]
 t=Table(data,colWidths=widths,repeatRows=1,hAlign='LEFT');t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),navy),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.HexColor('#EFF5F7'),colors.white]),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),6),('RIGHTPADDING',(0,0),(-1,-1),6),('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6),('LINEBELOW',(0,0),(-1,0),.8,teal)]));story.append(t);story.append(Spacer(1,10))
def img(name,height):story.append(Image(str(F/name),width=499,height=height));story.append(Spacer(1,8))
def new():story.append(PageBreak())
def fmt(v):return 'NA' if pd.isna(v) else f'{v:.4f}' if abs(v)>=.001 else f'{v:.2e}'
def footer(c,doc):
 c.setStrokeColor(teal);c.line(48,42,W-48,42);c.setFont('CN',8);c.setFillColor(gray);c.drawString(48,28,'IKZF2 / STAT1 | 冻结输入的阶梯诊断 | 2026-10-02');c.drawRightString(W-48,28,str(doc.page))
title(1,'完整诊断结论')
p('现有数据支持部分调控状态关联，尚不足以建立“巨噬细胞 IL18 → 肺 CD4 T 调控 → 细胞死亡”的完整机制链。','head')
p('此次完成基础模型、预定细胞阈值、固定样本协变量阶梯、反向顺序、逐项加入/移除和样本影响审计。IKZF2、STAT1 分开报告；Human Core 与 CD4-IL7R、小鼠 Core 与 ActCD4 不拼接。MRE11A 未恢复为主检验对象，Step8 未变更。')
table([['证据层','最宽基础模型（All）','可支持的结论'],['Human CD4-IL7R RNA / AUC','RNA：IKZF2 P=0.6604；STAT1 P=0.3435。IKZF2 AUC 3/3 种子均不显著；STAT1 仅 seed17 P=0.0202。','无稳定一致的双 TF 疾病组间效应。'],['MacAM IL18 → IKZF2','n=17：三种子 P=0.0838 / 0.2278 / 0.2477；严格≥20细胞后 n=11：P=0.0186 / 0.0186 / 0.0426。','依赖样本构成的名义关联；不是最宽人群中的稳健阳性。'],['同一 IL7R 范围内 IKZF2 → 模块','OXPHOS 与 PyroExec 为负相关；ApoExec 不显著。来源+组别+深度后 OXPHOS / PyroExec 仍三种子名义显著。','支持相关状态，不能解释为促进焦亡/凋亡。'],['Mouse ActCD4','n=6（3+3）；IKZF2 AUC 基础 P 三种子均<0.05；STAT1 2/3。RNA STAT1 上升，IKZF2 不显著。','小样本探索；外部 SCENIC4 每时间点仅1+1文库，只作方向验证。']],[105,236,158])
p('“名义阳性”仅用于可计算的 P&lt;0.05，并与 FDR 分开。全部新模型及低自由度标记都保留；不通过删样本、换阈值或挑种子构造结论。全 sender 阶梯新 BH 最小 FDR=0.0683，未过0.05；历史模型的 FDR 未改写。','small')
p('原有免疫图谱与注释仍属描述性基础；CD4 组成和 NicheNet 候选证据沿用已交付技术报告，本诊断没有重跑或扩展配体筛选。','small')
new();title(2,'细胞阈值改变了谁进入比较')
img('filter_thresholds.png',265)
table([['门槛变化','保留人数 / 对照+病例','此次排除的样本及原因'],['≥1 → ≥5','17（8+9）→15（7+8）','CTL092：sender=4；CWP149：sender=3'],['≥5 → ≥10','15（7+8）→13（7+6）','CWP127：sender=6；CWP94：sender=6'],['≥10 → ≥20','13（7+6）→11（6+5）','CTL106：receiver=10；CWP89：sender=16']],[104,143,252])
p('图：Human All / CD4-IL7R / MacAM IL18，两个端点均达到对应细胞数。每点为独立受试者汇总后的相关，不把细胞或3个随机种子视为重复样本。虚线仅标记名义 P=0.05。')
p('阈值≥1、5、10、20 在诊断前预定；它们是测量充分性与纳入范围的敏感性检查，不是用于寻找最小 P 的候选阈值。所有样本ID、组别、来源、细胞数、检测率和模型自由度见结果及成员清单。','small')
new();title(3,'CTL106：有影响证据，但不是错误样本证据')
table([['从13人中移除','n','seed17 ρ / P','seed29 ρ / P','seed43 ρ / P'],['不移除','13','0.5000 / 0.0819','0.4011 / 0.1744','0.4011 / 0.1744'],['仅 CTL106','12','0.7273 / 0.0074','0.7203 / 0.0082','0.6573 / 0.0202'],['仅 CWP89','12','0.4755 / 0.1182','0.3566 / 0.2551','0.3636 / 0.2453'],['两者均移除','11','0.6909 / 0.0186','0.6909 / 0.0186','0.6182 / 0.0426']],[103,29,122,122,123])
p('反向加回提供同一证据：从严格11人中仅加回 CTL106，P 变为0.1182 / 0.2551 / 0.2453；仅加回 CWP89，P 为0.0074 / 0.0082 / 0.0202。因此，13→11这一局部转折主要受 CTL106 影响，CWP89 的移除并未独立增强信号。')
p('但从最宽17人中只移除 CTL106，三种子 P=0.0165 / 0.0441 / 0.0718，仅2/3过0.05。不能把整个17→11的变化都归因于一个样本。已执行全受试者逐一移除，以及每个相邻阈值所排除样本的所有子集检查，避免只挑 CTL106。')
table([['样本','来源/组别','receiver / sender细胞','IKZF2 RNA检测率','seed17 IKZF2 AUC'],['CTL106','G227 / Control','10 / 1790','0%','0.014930'],['CWP89','G288 / CWP','41 / 16','0%','0.019236']],[79,119,119,92,90])
p('低细胞数会增加样本均值、检测率与排序的不确定性；RNA 未检出并不等于 regulon 活性为零。当前没有证据认定 CTL106 标注错误、污染或应删除。仅靠剔除无法区分测量噪声与真实异质性。','head')
p('依据：filter_transition_subset_removals.tsv（264行、全部预定子集）和 sender_stepwise_models.tsv 中 SENDER_TF_ONE_SAMPLE_REMOVAL。所有影响分析仍属于探索性敏感性分析。','small')
new();title(4,'固定样本后：最后一步主要消耗自由度')
img('fixed_roster_adjustment.png',265)
table([['固定11人模型','df','seed17/29 ρ / P','seed43 ρ / P'],['无协变量','9','0.6909 / 0.0186','0.6182 / 0.0426'],['G + S','6','0.8845 / 0.0035','0.8678 / 0.0052'],['G + S + R','5','0.8900 / 0.0073','0.8662 / 0.0117'],['G + S + R + D','4','0.8755 / 0.0223','0.8393 / 0.0367'],['G + S + R + D + E','3','0.8742 / 0.0526','0.8382 / 0.0762']],[145,40,157,157])
p('G=疾病组别；S=来源分类（此处2个参数）；R=receiver文库深度；D=sender文库深度；E=sender上皮负荷。相关模型对连续变量取秩再残差化。固定11人的 E 步骤几乎不改变ρ，却将 df 从4减为3并使P越过0.05；不能据此说“上皮校正消除了生物信号”。')
p('这只是局部解释。最宽17人的三个种子并未一致支持关联：全模型 P=0.0618 / 0.1549 / 0.1810。全部固定完整病例与逐模型可用病例无丢失，故此处不是缺失协变量导致人数变化。','small')
new();title(5,'来源调整改变效应，但来源与组别不独立')
table([['来源','最宽17人：对照 / 病例','严格11人：对照 / 病例'],['G135','1 / 0','1 / 0'],['G227','7 / 2','5 / 2'],['G288','0 / 7','0 / 3']],[99,200,200])
p('只有 G227 同时含对照和病例。来源调整既吸收测量/队列差异，也可能吸收来源相关生物差异；疾病组别的跨来源可比性有限。不能把来源调整当作额外独立样本或独立重复验证。')
table([['固定11人 IKZF2','seed17/29 ρ / P','seed43 ρ / P','df'],['仅加 S','0.8380 / 0.0048','0.7908 / 0.0112','7'],['全模型去 S','0.6070 / 0.1483','0.5241 / 0.2273','5'],['全模型去 R','0.8717 / 0.0237','0.8417 / 0.0356','4'],['全模型去 D','0.8898 / 0.0175','0.8671 / 0.0253','4'],['全模型去 E','0.8755 / 0.0223','0.8393 / 0.0367','4']],[134,162,162,41])
p('来源调整比最后一步 E 更实质：ρ 从0.6909升至0.8380（seed17/29），而不是仅P变化。秩协方差分解显示，总交叉乘积76.00，来源内63.95，来源间12.05，来源间分量并非负值；来源调整减少 IL18 秩方差22.3%、AUC秩方差38.1%，分母下降相对更大，残差相关因此增强。不能将其简单讲成“纠正了负向混杂”。')
p('逐项去除 R、D 或 E 都可使P重新过0.05，表明低自由度和协变量组合共同作用，不能唯一归因于上皮负荷。反向顺序与单独加入结果完整保留。')
p('NoDom 去除优势样本后 n=10：未调整 IKZF2 P=0.0537 / 0.0537 / 0.1173，GSR后P=0.0060 / 0.0060 / 0.0101；全模型又回到0.0516 / 0.0516 / 0.0793。NoParent16 在此选定严格样本中与 Matched 相同。来源/样本敏感性仍然存在。','small')
new();title(6,'同一 CD4-IL7R 内的表型关联')
img('same_scope_module_associations.png',265)
p('蓝色为负相关，红色为正相关；* 表示该模型名义 P&lt;0.05，未表示 FDR 通过。图中所有分析均为 Human All / CD4-IL7R / n=17，同一细胞范围与同一受试者汇总。全12个冻结模块、BG1、Core和其他模型均保留在完整表，不仅限于图中6个关注模块。','small')
table([['IKZF2 模块','基础 P（17 / 29 / 43）','G+S+R 后 P','解释'],['OXPHOS','0.0018 / 0.00014 / 0.0022','0.0046 / 0.00047 / 0.0028','持续负相关；基础新FDR三种子均<0.05'],['PyroExec','0.0125 / 0.0108 / 0.0376','0.0061 / 0.0044 / 0.0136','负相关；基础新FDR不通过'],['ApoExec','0.2644 / 0.4450 / 0.4168','0.4699 / 0.6933 / 0.8301','不构成名义阳性']],[88,158,148,105])
p('PyroExec 的 seed43 在加来源时 P=0.0507，加入深度后为0.0136；校正影响不是单调“压没”信号。STAT1 的相关分散且依赖种子，未形成相同的焦亡/凋亡证据。模块分数是转录状态汇总，不是细胞死亡率或功能测定。')
p('IL18 受体链仍有限：原 Human All / IL7R ≥20细胞的16个受体样本中，IL18R1+IL18RAP 双检出中位数0%，范围0-2.13%。这是受体表的分母16，不是sender配对的11人。RNA稀疏不能证明受体蛋白缺失，也不能反过来证明有效通讯。','small')
new();title(7,'小鼠身份调整与小样本边界')
table([['Mouse All / ActCD4 AUC','基础P','加深度 R 后P','R+TregNoTF 后P'],['IKZF2 seed17','0.0384','0.0642','0.0033'],['IKZF2 seed29','0.0312','0.0445','0.0037'],['IKZF2 seed43','0.0071','0.0087','0.0150'],['STAT1 seed17','0.0614','0.0849','0.0216'],['STAT1 seed29','0.0465','0.0799','0.0058'],['STAT1 seed43','0.0157','0.0382','0.0133']],[176,107,107,109])
p('n=6（3对照+3硅肺文库）始终不变。T 为既有冻结 TregNoTF 身份分数，仅用于 ActCD4 TF组间模型；没有把T泛加到所有终点。R+T后残余df=2，调整可以增强信号，但估计对个别样本非常敏感，不可包装为强机制证据。')
p('小鼠 ActCD4 基础 RNA：IKZF2 P=0.9538；STAT1 P=0.0160且上升。加入 R+T 后分别P=0.8780和0.0196。六个关注的线粒体/焦亡/凋亡模块基础组间差异均不显著；不同TF、RNA与AUC层不能互相替代。')
p('小样本P的处理','head')
p('未调整 Spearman 在 n≤9 时使用全部排列的双侧精确 P；较大n用t近似。调整后的秩相关保留正自由度下的探索性t近似，并明确标记低n/低df。它并不是小样本精确条件检验。')
p('检测到11个唯一小鼠调整规格在 n=6、df=2或3时 |ρ|≈1；其t近似P可小至10^-23，受离散排序、残差几何及数值下限影响，不能按该数量级解读证据强度。全部列于 small_n_near_perfect_partial_rank_flags.tsv；原值保留，未用它们支撑核心结论。')
p('SCENIC4 的 GSE183682 第7/56天各只有1对照+1硅肺文库，只作冻结 regulon 方向核对，不给有效组间P、不将时间点或种子算作生物重复。','small')
new();title(8,'原 edgeR 与新 meanCP10k 检验并列')
table([['终点 / 方法','IKZF2：效应 / P / FDR','STAT1：效应 / P / FDR'],['Human Core 原 edgeR<br/>来源调整，log2FC','1.5140 / 0.0169 / 0.8388','-0.4410 / 0.4071 / 0.9977'],['Human Core 新基础OLS<br/>log1p平均CP10k系数','-0.0141 / 0.9469 / 0.9469','-0.0964 / 0.4461 / 0.9469'],['Mouse ActCD4 原 edgeR<br/>group，log2FC','-0.2068 / 0.7781 / 0.9703','1.0424 / 0.0037 / 0.4179'],['Mouse ActCD4 新基础OLS<br/>log1p平均CP10k系数','0.0229 / 0.9538 / 0.9538','0.4293 / 0.0160 / 0.0256']],[187,156,156])
p('以上新OLS数值应以完整表为准：FDR分别来自原全基因检验和新预定TF诊断家族，不能横向互换。Human原Primary也保留：IKZF2 log2FC=1.5623，P=0.0061，FDR=0.5083；不能用这个更宽范围替代IL7R的结果。','small')
p('为何两种检验会不同','head')
p('原 edgeR 在既有计数聚合、归一化及离散度框架下做准似然计数模型；保留其原始细胞/样本定义、模型与全基因FDR。新基础OLS对本次冻结评分细胞中的每个受试者原文库CP10k均值取log1p，再做等权受试者组间比较；人类基础模型未加来源，另有同范围来源/深度阶梯。')
p('两者的聚合方式、归一化、均值-方差假设、协变量、可用细胞范围和效应单位均可能不同。方向或P不一致不等于计算错误，也不能凭新均值OLS否定原edgeR。此次只并列审计，没有重新拟合edgeR或把其log2FC与OLS系数直接相减。')
p('若要解释差异，应在明确授权后固定完全相同的细胞与样本，审计原计数、offset/归一化、离散度及来源设计，再开展限定的计数模型复核；不以得到显著为目标。','small')
new();title(9,'复现、交付与可执行的下一步')
qa=json.loads((D/'provenance/full_reproduction_audit.json').read_text())
table([['复现项目','规格数','核验结果'],['基础模型 vs 阶梯起点','760','样本、效应、P、df全部一致'],['原 sender 调整 / 未调整','144 / 72','ρ最大差1.72×10^-15；原有效P最大差5.00×10^-16'],['原小鼠身份调整 / S2 AUC组间','48 / 96','全部效应、原有效P、n、df一致'],['原 S2 TF-模块 / M2模块组间','1152 / 96','全部效应、原有效P、n、df一致'],['阶梯全量','31744 receiver + 6720 sender','成员清单、阈值、设计秩、df、P范围和新BH核验PASS'],['冻结输入','22个源文件','服务器SHA256全部未变化']],[176,119,204])
p('全量包包含760条基础模型、38464条阶梯模型、264条相邻阈值子集记录、完整受试者成员表、22项冻结输入清单、逐项复现差值、方法、脚本及本报告。阶梯重复的同数据/设计在各自新FDR家族内只计一次（receiver 5936；sender含独立LOO家族共3456）。')
p('下一步选择（仅建议，未自动执行）','head')
p('1. 按现有五项交付目标形成客户版结论：展示免疫/CD4结构、候选通讯和分开的两TF状态证据，明确完整机制链尚未闭合。<br/>2. 如需弄清RNA检验差异，授权同细胞、同样本的原edgeR计数/设计核验，不换基因、不换数据集。<br/>3. 后续验证优先解决同一CD4范围内受体蛋白/功能响应与独立匹配样本，再讨论功能扰动和死亡读出；本次未开展新实验或机制扩展。')
p('服务器交付根目录：/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/stepwise_diagnostic_20261002_v1。本报告主要来源：results/baseline_models.tsv、sender_stepwise_models.tsv、receiver_stepwise_models.tsv、filter_transition_subset_removals.tsv、source_group_design.tsv；复现证据：provenance/full_reproduction_audit.json。原模型见包内reference_tables和逐项comparison表。','small')
pdf=O/'lung_IKZF2_STAT1_full_stepwise_diagnostic.pdf'
SimpleDocTemplate(str(pdf),pagesize=A4,leftMargin=48,rightMargin=48,topMargin=40,bottomMargin=56,title='IKZF2 STAT1 全量阶梯诊断',author='Research audit').build(story,onFirstPage=footer,onLaterPages=footer)
print(pdf)
