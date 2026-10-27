from pathlib import Path
import json,html,pandas as pd,numpy as np
from scipy import stats
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reportlab.platypus import SimpleDocTemplate,Paragraph,Table,TableStyle,Spacer,PageBreak,Image
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from pypdf import PdfReader,PdfWriter
O=Path(r'C:\Users\31937\lung_append_analysis_20261017_v18');T=O/'tables';V=O/'figures';V.mkdir(exist_ok=True);D=Path(r'C:\Users\31937\Downloads\线粒体焦亡_v17追加v18分析_20261017');D.mkdir(exist_ok=True)
r=pd.read_csv(T/'all_associations.tsv',sep='\t');score=pd.read_csv(T/'sample_gene_and_program_scores.tsv',sep='\t');plan=json.loads((O/'protocol_frozen.json').read_text());updates=pd.read_csv(O/'review_tables/相对v17的逐基因证据更新.tsv',sep='\t');source=pd.read_csv(T/'source_competition_summary.tsv',sep='\t');ef=pd.read_csv(T/'programme_group_effects.tsv',sep='\t')
def row(g,sc,y,co='Joint',v='Raw',model='JointStratum',kind=None,sender=None,x=None):
 z=r[r.cohort.eq(co)&r.variant.eq(v)&r.model.eq(model)&r.scope.eq(sc)&r.x.eq(x or 'G_'+g)&r.y.eq(y)]
 if kind:z=z[z.kind==kind]
 if sender:z=z[z.sender==sender]
 assert len(z)==1,(g,sc,y,co,v,model,len(z));return z.iloc[0]
def fmt(v):return 'NA' if pd.isna(v) else f'{v:.3g}'
def statline(z):return f'rho={z.rho:+.3f}；P={fmt(z.P)}；FDR={fmt(z.FDR)}'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
labels=[g+'/'+sc for g,sc in plan['focus_pairs']];programs=plan['phenotypes']+plan['damage']
for model,filename in [('JointStratum','primary_heatmap.png'),('JointStratumCondition','adjusted_heatmap.png')]:
 vals=np.array([[row(g,sc,'M_'+m,model=model).rho for m in programs] for g,sc in plan['focus_pairs']]);fig,ax=plt.subplots(figsize=(9,4.5));im=ax.imshow(vals,cmap='RdBu_r',vmin=-1,vmax=1,aspect='auto');ax.set_xticks(range(len(programs)),programs,rotation=35,ha='right');ax.set_yticks(range(len(labels)),labels)
 for i,(g,sc) in enumerate(plan['focus_pairs']):
  for j,m in enumerate(programs):
   z=row(g,sc,'M_'+m,model=model);mark='**' if z.FDR<.05 else '*' if z.P<.05 else '';ax.text(j,i,f'{z.rho:.2f}{mark}',ha='center',va='center',fontsize=7,color='white' if abs(z.rho)>.72 else 'black')
 fig.colorbar(im,ax=ax,label='Partial rank correlation');ax.set_title('Joint n=10 | '+model+' | * permutation P<0.05; ** FDR<0.05');fig.tight_layout();fig.savefig(V/filename);plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(8,3.8))
for ax,g in zip(axs,['TTC39C','MDFIC']):
 d=score[(score.variant=='Raw')&score.scope.eq('ActCD4')&score.dataset.isin(['GSE314943','GSE183682'])].sort_values(['dataset','sample']);lab=d.dataset.astype(str)+'_'+d.day.astype(str);C=np.column_stack([np.ones(len(d)),pd.get_dummies(lab,drop_first=True).to_numpy(float)]);M=np.eye(len(d))-C@np.linalg.pinv(C);xx=M@stats.rankdata(d['G_'+g]);yy=M@stats.rankdata(d.M_ApoExec)
 for group,color in zip(sorted(lab.unique()),['#277da8','#ef9c32','#27876d']):
  for cond,marker in [('Saline','o'),('Silica','^')]:
   take=(lab==group)&(d.group==cond);ax.scatter(xx[take],yy[take],c=color,marker=marker,s=50,label=group.replace('GSE','')+' '+cond)
 z=row(g,'ActCD4','M_ApoExec');ax.set_title(f'{g} / ActCD4\nrho={z.rho:.3f}, P={z.P:.4f}, FDR={z.FDR:.4f}');ax.set_xlabel('Candidate residual ranks');ax.set_ylabel('ApoExec residual ranks')
axs[1].legend(fontsize=6,loc='upper left',bbox_to_anchor=(1.0,1));fig.tight_layout();fig.savefig(V/'historical_scatter.png');plt.close(fig)
pdfmetrics.registerFont(TTFont('CN',r'C:\Windows\Fonts\ARIALUNI.TTF'));styles={k:ParagraphStyle(k,fontName='CN',fontSize=fs,leading=ld,spaceAfter=7,wordWrap='CJK',textColor=colors.HexColor('#185568') if k in ['h','h2'] else colors.black) for k,fs,ld in [('h',17,25),('h2',11.5,17),('b',9.8,15.5),('s',8,12),('c',8.1,11.7)]};story=[]
def p(t,k='b'):return Paragraph(html.escape(str(t)).replace('\n','<br/>'),styles[k])
def add(t,k='b'):story.append(p(t,k))
def page(t):story.append(PageBreak());add(t,'h')
def table(head,rows,widths):
 a=Table([[p(x,'c') for x in head]]+[[p(x,'c') for x in rr] for rr in rows],colWidths=widths,repeatRows=1);a.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#d7e9ec')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f2f6f7')]),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),5),('RIGHTPADDING',(0,0),(-1,-1),5),('TOPPADDING',(0,0),(-1,-1),4),('BOTTOMPADDING',(0,0),(-1,-1),4)]));story.append(a);story.append(Spacer(1,8))
def img(n,w=490):
 from PIL import Image as PI
 a=PI.open(V/n);story.append(Image(str(V/n),width=w,height=w*a.height/a.width));story.append(Spacer(1,8))
add('v18 追加分析：补齐连接与检验替代解释','h');add('接续原v17第1-9页 | 2026-10-17 | 供独立AI复审','s')
add('阅读说明：前9页保留原v17基准。以下记录新增分析及其对候选解释的影响，不能只读旧优先名单，也不能把本补充的阴性结果当作否定全部候选。新增结果未证明湿实验机制。')
add('本轮实际完成9基因、10个基因-状态组合、18条上游路径。新增2,544项关联（含主分析与敏感性）、221项程序组间检验、28项定向TF组间检验，以及500项共同基因覆盖敏感性。联合主分析的131项候选-程序/代谢检验中，30项原始P<0.05，1项FDR<0.05；不能把不同模型重复命中累加成独立验证。')
table(['对象','新增支持','需要调整的解释'],[['TTC39C、MDFIC','ActCD4内与ApoExec出现较强转录关联','值得提升为第二梯队损伤关联备选；上游和独立疾病方向仍不完整。'],['SKIL','与Treg活化程序的关联在组别调整后仍为名义显著','更支持活化/状态调节；p53及损伤连接未被补齐。'],['IL1R2','Treg ENO1 RNA关联通过联合模型相应家族FDR','不能由RNA正相关证明抑制ENO1酶活或抗焦亡功能。'],['FOXP1','多种T细胞状态程序及部分损伤程序有名义关联','Core-PyroExec关系在固定组成后明显减弱。'],['DDIT4、GSN','有部分表型/单段近端线索','本次没有补成所选巨噬上游到损伤的连续链。'],['AHNAK、ST6GALNAC3','分别有OXPHOS、MitoROS转录联系','OXPHOS和氧化应激响应不是损伤实测，继续保留待补。']],[92,190,208])
add('本补充不重新清空候选池，不用FDR、单次敏感性P或NicheNet排名作统一排除门槛。新增证据用于修正研究问题和验证优先级。','s')
page('补充1  方法先冻结，评分和调控活性去除重叠')
table(['环节','实际口径'],[['对象与范围','IL1R2/Treg；DDIT4/ActCD4；SKIL/Treg；FOXP1/Core和Treg；GSN/ActCD4；AHNAK/Naive；TTC39C、MDFIC、ST6GALNAC3/ActCD4。'],['表型/损伤程序','复用v14固定成员：TAct、TReg、Inhib、Naive、Cycle；MitoROS、MitoISR、ApoExec、PyroExec；OXPHOS、Inflam、Mitophagy为辅助。另固定GO:0006096鼠糖酵解成员。'],['新评分','原始计数按独立样本汇总；mean log1p(UMI/library×10000)。从全部评分中移除9候选；映射至少3基因且覆盖至少50%才可评分。此v18分数不替换v14匹配对照分数。'],['调控活性','原pySCENIC 0.12.1网络、seed17、原9,877基因排名输入；重算12个调控子AUCell（阈值0.05、无权重、无最大值归一化）。移除候选、全部评分成员和TF自身。'],['设计','疾病：原六小鼠及两队列分层联合n=10。第三队列TLike n=12单独时间分析，另排除8细胞样本作n=11敏感性。GSE183682每时间点1+1只给描述。'],['敏感性','疾病组别调整、接收端BG1、Core各状态1/3固定组成、留一方向、来源竞争；三种覆盖不同模块额外使用共同基因。']],[100,390])
add('TAct/TReg/Inhib等是转录程序，不等同功能性活化、抑制能力或耗竭实测。MitoROS是氧化应激响应，OXPHOS是相关转录表达；ApoExec/PyroExec不证明caspase/gasdermin已激活。糖酵解GO评分含调控成员，不是代谢通量。')
add('统计：对变量秩进行协变量残差化后求相关。联合主模型控制队列/时间层；组别调整模型另加疾病状态。置换原始秩后重新残差化，联合主模型层内2,880种，组别调整层内36种；原六样本720/36种。时间组内1,296种，n=11为432种；不分时间的n=12/11使用固定种子10,000次Monte Carlo。P按双侧绝对相关计算，近似t模型P附表单列。','s')
add('组间程序检验：六样本20种分配（最小双侧P=0.10），联合80种分层分配（最小0.025），时间总体用F模型。完整同类假设家族独立BH。置换依赖指定块内可交换性；细胞数和测量精度不等，小样本结果为探索性。联合组别调整时第二队列四个样本均为单例单元，不能贡献组内交换，也不算独立重复验证。','s')
page('补充2  主分析与组别调整呈现不同证据层次')
img('primary_heatmap.png');img('adjusted_heatmap.png')
add('图中每格为预设联合相关；*原始置换P<0.05，**该完整家族FDR<0.05。上图控制队列/时间，下图再控制疾病组别。后者用于评价共同疾病状态解释，不作为所有候选必过门槛。完整131项候选检验以及80项状态程序-损伤程序检验均保留。','s')
page('补充3  TTC39C与MDFIC获得新的项目内凋亡转录线索')
img('historical_scatter.png')
for g in ['TTC39C','MDFIC']:
 z=row(g,'ActCD4','M_ApoExec');zz=row(g,'ActCD4','M_ApoExec',model='JointStratumCondition');add(g+'：联合主模型 '+statline(z)+'；组别调整 '+statline(zz)+'。','b')
add('两者联合留一分析方向均保留；原六样本BG1中也保留正向关联。MDFIC组别调整后的效应仍大，但36种条件排列限制了P值分辨率；不应仅因P=0.0833把该方向删除。相关性是候选与凋亡转录程序的联系，不能声称基因引发凋亡。')
add('独立时间队列中，TTC39C与MitoISR在日调整后rho=0.733、P=0.0185、FDR=0.567；去掉8细胞样本后rho=0.856、P=0.0231。共同基因评分保持该方向。该队列为TLike，不能称为ActCD4复现；MDFIC-ApoExec时间日调整P=0.444，未获得同等支持。')
add('关键矛盾：MDFIC在GSE183682的D7/D56疾病方向均与原候选方向相反。TTC39C D7对聚合口径敏感：旧v14细胞聚合差=-0.0704；v16 TMM-logCPM差=+0.1398；本轮log1p-CP10k差=+0.0293。D56均为负。各口径单位不同，只比较方向；不能选择一种口径宣称完全重复。')
add('更新建议：将TTC39C、MDFIC列入损伤转录关联的第二梯队备选，提交独立评审；目前没有补齐TNF/TNFSF9到候选的上游链。旧v17 A2结论保留作历史基准，本页是新增证据更新。','s')
page('补充4  FOXP1更明确关联状态重塑，Core焦亡关系组成敏感')
rows=[]
for v,model,label in [('Raw','JointStratum','联合：控制层'),('Raw','JointStratumCondition','联合：再控制疾病'),('FixedComp','JointStratum','联合：固定组成'),('FixedComp','JointStratumCondition','联合：固定组成+疾病')]:
 z=row('FOXP1','Core','M_PyroExec',v=v,model=model);rows.append([label,f'{z.rho:+.3f}',fmt(z.P),fmt(z.FDR)])
table(['Core FOXP1-PyroExec模型','rho','P','FDR'],rows,[245,70,80,95])
add('原六样本的FOXP1-PyroExec rho=-1，P=0.00278；BG1后相同，但固定组成后rho=-0.657、P=0.175；再控制疾病rho=-0.164、P=0.833。背景子集保留并不解决组成解释。')
add('联合FOXP1与TAct、TReg、Inhib负相关，与Naive正相关；固定组成后TAct/Inhib/Naive仍保留方向，并在较小的预设敏感性家族内通过FDR。不能用该敏感性FDR重标原主模型。FOXP1-ApoExec联合P=0.0396，但控制疾病后rho约0、P=1。')
add('FOXP1 RNA与去重叠自身AUCell：原Core rho=0.829、P=0.0583；组别调整rho=0.250、P=0.722。只有BG1+固定组成的特定模型达到名义/FDR支持，不能概括为“表达下降稳定等于调控活性下降”。')
add('第三队列FOXP1-PyroExec pooled rho=-0.789、P=0.00400；日调整rho=-0.671、P=0.0571。至少10细胞敏感性日调整rho=-0.849、P=0.0208。保留此独立时间线索，但主分析和敏感性分开解释。')
add('PyroExec仍包含GSDMD及IL1B/IL18等低丰度信号，既往已提示背景RNA风险。本轮BG1子集不是环境RNA去污染或蛋白定位；不能由相关性确认纯CD4自身执行焦亡。')
add('更新定位：FOXP1仍适合研究T细胞状态和损伤易感性；不能把原Core相关性当作已经证明纯T细胞内源焦亡调控。SEMA3E/CD200的竞争来源也仍未解决。')
page('补充5  IL1R2补上代谢关联，SKIL更贴近Treg活化')
add('IL1R2 / Treg','h2');z=row('IL1R2','Treg','G_ENO1');zz=row('IL1R2','Treg','G_ENO1',model='JointStratumCondition');add('ENO1 RNA联合：'+statline(z)+'。再控制疾病：'+statline(zz)+'。原六样本及BG1的rho均为+1；条件排列P均为0.0556。')
add('IL1R2还与Inhib、MitoROS、ApoExec及糖酵解程序有联合名义关联，但疾病调整后多数减弱。这里观察的是RNA共同变化；已有文献中的IL1R2抑制ENO1酶活不能由正RNA相关确认或反驳。保护性反馈、活化伴随变化和共同疾病状态仍需区分。')
add('本轮IL1R2-PyroExec联合：'+statline(row('IL1R2','Treg','M_PyroExec'))+'。当前新增支持主要在代谢和状态联系，不能据此声称抗焦亡功能已被项目数据确认。')
add('IL1A/B-IL1R2按诱饵结合入口处理，未补造“受体激活NF-κB”模型。相应配体至候选的样本关联未形成名义支持的连续链。')
add('SKIL / Treg','h2');z=row('SKIL','Treg','M_TAct');zz=row('SKIL','Treg','M_TAct',model='JointStratumCondition');add('TAct联合：'+statline(z)+'。再控制疾病：'+statline(zz)+'。原六样本及BG1均保留强正方向。')
add('SKIL-ApoExec联合P=0.0747，疾病调整后rho=-0.094、P=0.917；SKIL-p53 AUCell原样本rho=0.086、P=0.919，调整后P=0.222。当前项目内证据更支持Treg活化联系，尚未补上p53/内源性凋亡段。')
add('CCL2-SKIL的名义关联来自AMlike（rho=0.943，P=0.0167，FDR=0.150）；AllMac对应P=0.136。AMlike组别调整后P=0.278。不能把全肺Mac的CCL2来源优势与AMlike的相关性拼成同一已证实链。')
page('补充6  DDIT4等候选的新增阴性与弱证据必须保留')
table(['对象','新增结果','解释'],[['DDIT4 / ActCD4','与Inhib、Naive、Cycle有名义关联；MitoISR调整rho=0.888、P=0.0556；PyroExec主rho=0.126、P=0.751','应激关联有方向线索；当前数据不支持其RNA与PyroExec评分关联。'],['TNF-DDIT4近端','AllMac TNF-DDIT4 rho=-0.257、P=0.658；RELB-DDIT4 rho=0.886、P=0.0333、FDR=0.267；配体-近端段未成连续名义链','单独RELB段不补成TNF→DDIT4→焦亡机制。'],['GSN / ActCD4','Inhib名义关联；ApoExec主rho=-0.459、P=0.264；PROS1-GSN仅AMlike pooled名义支持，调整P=0.278','保持备选；没有新增足够的损伤连接。'],['AHNAK / Naive','OXPHOS主rho=0.811、P=0.0417；调整P=0.0833；NFATC1/2关联未达名义支持','氧化磷酸化转录联系不等于线粒体受损；保留表型问题。'],['ST6GALNAC3 / ActCD4','MitoROS主rho=0.872、P=0.00972；调整rho=0.569、P=0.222','补上应激响应线索，但共同疾病状态解释仍在。']],[100,245,145])
add('本轮AUCell清理后，定向TF组间检验共28项；三对三20种分配的最小P=0.10，故不能以这些组间P未达0.05推断“没有调控响应”。原始效应、AUCell输入和完整家族随附。')
add('组间程序检验也有分辨率限制：联合80种分配中，部分T细胞表型和MitoROS/ApoExec程序P=0.025，但完整家族FDR未达0.05。可保留名义程序变化，不能写成已观察到死亡执行。')
page('补充7  上游逐段比较保持同一发送端与接收状态')
rows=[]
for gene,sc,lig in [('DDIT4','ActCD4','TNF'),('SKIL','Treg','CCL2'),('FOXP1','Treg','SEMA3E'),('GSN','ActCD4','PROS1')]:
 for send in ['AllMac','AMlike']:
  z=row(gene,sc,'G_'+gene,co='314',model='Pooled',kind='ligand_candidate',sender=send,x='L_'+lig);zz=row(gene,sc,'G_'+gene,co='314',model='Condition',kind='ligand_candidate',sender=send,x='L_'+lig);rows.append([lig+' → '+gene+' / '+sc,send,f'{z.rho:+.3f}',fmt(z.P),fmt(z.FDR),f'{zz.rho:+.3f} / {fmt(zz.P)}'])
table(['连接（非因果箭头）','发送端','rho','P','FDR','调整rho / P'],rows,[165,65,50,60,60,90])
add('表中是原六小鼠配对样本相关；AMlike是独立列明的敏感性来源范围，不与AllMac拼接。TNF-NF-κB和CCL2/CXCL10-AP1使用的是共享转录响应兼容性，不代表唯一受体机制。FOXP1和GSN路线未预设特异近端TF，缺失处保持NA。')
add('全部18条预设路线，包括历史候选的TNF/TNFSF9和AHNAK的H2-EB1入口，均有逐段表。H2-EB1/CD4仅为抗原呈递相关的先验入口，不提供抗原或TCR特异性，也不等同可溶性配体作用。')
add('当前未形成一条在相同来源和状态中有充分连续支持的“配体—近端TF—候选—损伤”链。这里报告具体缺段，不把所有候选整体判为无研究价值。')
page('补充8  受体检出增加可评估性，竞争来源仍影响命题')
table(['接收范围与受体','Silica三样本受体阳性/总细胞','阳性中T标记共检出'],[['ActCD4 TNFRSF1B','88/284；32/100；50/225','88；32；50'],['Treg CCR2','38/178；37/156；42/161','37；37；42'],['Treg CXCR3','25/178；13/156；23/161','25；13；23'],['Treg CD200R1','16/178；34/156；53/161','16；32；53'],['Treg PLXND1','18/178；23/156；13/161','17；23；13'],['ActCD4 AXL','16/284；4/100；13/225','16；4；13']],[145,205,140])
add('对216个路线-样本-过滤组合核查了全部必需受体亚基同细胞检出及T标记。RNA共检出有助确认可评估性，但不是蛋白定位或占据证明；掉落未检出也不作绝对否定。')
rows=[]
for lig in ['IL1A','IL1B','TNF','CCL2','CXCL10','SEMA3E','CD200','PROS1']:
 z=source[(source.gene==lig)&source.Broad.eq('Mac')].iloc[0];top=source[(source.gene==lig)&source.samples.eq(3)].sort_values('mean_CP10k',ascending=False).iloc[0];rows.append([lig,f'{z.mean_CP10k:.3g}',int(z.rank_three_samples),top.Broad,f'{top.mean_CP10k:.3g}'])
table(['配体','Mac均值CP10k','Mac排名','最高的三样本谱系','该谱系均值'],rows,[70,100,65,145,110])
add('来源排名限三个病例均有≥10细胞的全肺Broad标签；少样本标签另表保留。全肺Broad=Mac与通讯冻结AllMac范围不完全相同，来源表仅作描述，不能与不同范围P值拼接。IL1A、CCL2在此表靠前；IL1B、SEMA3E、CD200、PROS1均有明显竞争来源。','s')
page('补充9  跨队列支持与评分覆盖敏感性分别保存')
rows=[]
for g,m in [('FOXP1','PyroExec'),('TTC39C','MitoISR'),('MDFIC','ApoExec'),('SKIL','PyroExec')]:
 a=row(g,'TLike','M_'+m,co='274',model='TimePooled');b=row(g,'TLike','M_'+m,co='274',model='TimeDayAdjusted');c=row(g,'TLike','M_'+m,co='274_min10',model='TimeDayAdjusted');rows.append([g+' - '+m,f'{a.rho:+.3f}\n{fmt(a.P)} / {fmt(a.FDR)}',f'{b.rho:+.3f}\n{fmt(b.P)} / {fmt(b.FDR)}',f'{c.rho:+.3f}\n{fmt(c.P)} / {fmt(c.FDR)}'])
table(['TLike时间队列连接','全12样本rho\nP / FDR','日调整rho\nP / FDR','≥10细胞日调整rho\nP / FDR'],rows,[130,120,120,120])
add('第三队列没有后期同期盐水对照，TLike也不是纯ActCD4/Treg。上述结果是独立时间线索，不把全部时间相关都视为同一疾病机制复现。MDFIC和SKIL在时间队列与TReg/Inhib表型程序的部分pooled结果通过FDR，仍需区分时间/组成解释。')
add('覆盖检查发现MitoISR、Mitophagy、Glycolysis在队列间映射成员数不同。追加共同成员敏感性分别使用15、10、117基因，重算500个受影响对比。TTC39C-MitoISR日调整方向和P保留；该敏感性有自己的较小检验家族，不将其FDR替换主模型。')
add('第二队列两个时间点均为1+1，只提供分样本方向。归一化口径对照表同时保留v14细胞聚合、v16 TMM-logCPM和v18 log1p-CP10k；不同单位不能直接比较数值大小。TTC39C早期符号变化是归一化敏感性，MDFIC反向、其他候选弱变化等均保留。')
add('这些检查用于判断支持的适用范围；没有要求每个队列、每种敏感性都通过P<0.05，也没有因为某个分支更显著而替换预设主结果。')
page('补充10  独立评审入口与追加分析完成范围')
table(['复核内容','对应文件'],[['预设范围与统计方案','protocol_frozen.json；coverage_sensitivity_record.json'],['全部主分析及敏感性','tables/all_associations.tsv（2,544项）；programme_group_effects.tsv（221项）；focal_TF_disease_effects.tsv（28项）'],['共同成员复评','tables/common_gene_sensitivity.tsv（500项）及common_gene_member_audit.tsv'],['可复算输入','sample_gene_and_program_scores.tsv；sender_ligand_expression.tsv；clean_TF_sample_AUC.tsv；分样本表达、成员与覆盖表'],['逐段链条和版本更新','review_tables/上游连接逐段证据_不跨来源拼接.tsv；相对v17的逐基因证据更新.tsv；候选表型损伤_主分析与敏感性对照.tsv'],['核查与原版','tables/numerical_QA.tsv；QA.json；baseline_v17原报告/附件及哈希']],[115,375])
add('已完成176个新增/敏感性BH家族重算，4项关键置换P独立逐排列复算；候选与评分成员未进入新的TF签名；原v17报告哈希未变。完整代码、运行日志、输入清单和缺口均附上。原细胞矩阵及细胞AUCell对象保留服务器，不在下载目录重复复制。')
add('建议独立AI重点检查：①置换块和四个单例样本的解释；②去重叠及跨队列评分覆盖；③主分析/敏感性FDR是否混用；④FOXP1组成敏感性的结论；⑤IL1R2 RNA与酶活作用的区别；⑥AMlike与AllMac是否被混拼；⑦TTC39C归一化和MDFIC反向结果；⑧是否把转录程序写成真实死亡。')
add('交付前计算层面的新增证据已补齐到本轮预定范围。当前可以交付可检验、可否定的候选问题，但不能以这些相关性宣布B2机制成立。后续湿实验应区分候选作用与伴随标记，并检验实际来源、候选干预和目标终点；不把“尚未完成湿实验”当作筛选失败。','s')
appendix=O/'v18_appendix.pdf'
def foot(c,d):c.setFont('CN',8);c.setFillColor(colors.HexColor('#536b77'));c.drawString(48,26,'v18 追加分析 | 原v17第1-9页保持原文 | 内部复审 | 2026-10-17');c.drawRightString(547,26,str(d.page+9))
SimpleDocTemplate(str(appendix),pagesize=A4,leftMargin=48,rightMargin=48,topMargin=42,bottomMargin=42,title='v18追加分析').build(story,onFirstPage=foot,onLaterPages=foot)
old=Path(r'C:\Users\31937\Downloads\线粒体焦亡_双基准复评_v17_20261017\v17_双基准二次评估报告.pdf');writer=PdfWriter();writer.append(str(old),import_outline=False);writer.append(str(appendix),import_outline=False);writer.add_outline_item('v17基准原文（第1-9页）',0);writer.add_outline_item('v18追加分析（从第10页起）',9);writer.add_metadata({'/Title':'v17结果及v18追加分析 - 独立复审版','/Author':'项目分析'});dest=D/'v17结果_追加v18新分析.pdf'
with open(dest,'wb') as f:writer.write(f)
(O/'delivery.json').write_text(json.dumps({'pdf':str(dest),'appendix':str(appendix),'directory':str(D),'baseline_pdf':str(old),'baseline_pages':len(PdfReader(old).pages),'appendix_pages':len(PdfReader(appendix).pages)},ensure_ascii=False,indent=2),encoding='utf-8');print(dest,'APPENDIX',len(PdfReader(appendix).pages),'TOTAL',len(writer.pages))