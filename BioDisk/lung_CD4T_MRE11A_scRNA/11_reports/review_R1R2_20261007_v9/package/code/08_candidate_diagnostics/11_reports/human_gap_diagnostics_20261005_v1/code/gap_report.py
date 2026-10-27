from pathlib import Path
import os,json,hashlib,base64,html,zipfile,datetime
os.environ['OPENBLAS_NUM_THREADS']='1'
import pandas as pd,numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/human_gap_diagnostics_20261005_v1';R=O/'results';F=O/'figures';D=O/'deliverables';V=P/'11_reports/human_alignment_20261005_v1/results'
F.mkdir(exist_ok=True);D.mkdir(exist_ok=True)
def rd(n):return pd.read_csv(R/n,sep='\t',low_memory=False)
def wr(d,n):pd.DataFrame(d).to_csv(R/n,sep='\t',index=False,na_rep='NA')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for c in iter(lambda:f.read(8*1024*1024),b''):h.update(c)
 return h.hexdigest()
up=rd('all_candidate_TF_stepwise_models.tsv.gz');fin=up[up.step=='SourceBothDepth'];nom=up[up.P<.05]
assert len(up)==12600 and len(nom)==997 and (fin.P<.05).sum()==172 and (fin.FDR<.05).sum()==0
assert up.groupby(['receiver','tf','readout','variant','seed','sender','ligand']).n.nunique().max()==1
flow=rd('IL15_fixed_condition_attrition.tsv');meas=rd('IL15_measurements_with_descriptive_sampling_intervals.tsv.gz');con=rd('IKZF2_matched_reference_contrast_models.tsv');chains=rd('selected_nominal_chains_reference_contrast.tsv');influence=rd('excess_program_influence_summary.tsv');gz=rd('GZMA_TReg_independence_disease_models.tsv');comp=rd('GZMA_IKZF2_source_composition_summary.tsv');controls=rd('matched_control_specificity_models.tsv');ind=rd('IL7R_independence_state_diagnostics.tsv');state=rd('state_quartile_standardized_donor_models.tsv')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
colors={'Control':'#357d9e','CWP':'#c16b40'}
fig,axes=plt.subplots(1,2,figsize=(13,5.8),gridspec_kw={'width_ratios':[1.35,1]})
contexts=[('MacAM','IL15'),('MacAM','IL15RA'),('MonoMac','IL15'),('MonoMac','IL15RA'),('MacWide','IL15'),('MacWide','IL15RA'),('CD4-IL7R','IL2RB'),('CD4-IL7R','IL2RG')]
for i,(scope,gene) in enumerate(contexts):
 for j,g in enumerate(['Control','CWP']):
  q=meas[(meas.scope==scope)&(meas.gene==gene)&(meas.group==g)&(meas.cells>=20)].sort_values('pct');x=i+(j-.5)*.23+np.linspace(-.045,.045,len(q));axes[0].scatter(x,q.pct*100,s=22,color=colors[g],label=g if i==0 else None,alpha=.85)
axes[0].axhline(10,color='#999',ls='--',lw=.8);axes[0].axhline(5,color='#aaa',ls=':',lw=.8);axes[0].set_xticks(range(8),[s.replace('CD4-','')+'\n'+g for s,g in contexts],rotation=45,ha='right');axes[0].set_ylabel('Observed RNA-positive cells (%)');axes[0].set_title('Continuous donor measurements (>=20 cells)',loc='left',fontweight='bold');axes[0].legend(frameon=False)
stages=['PairedCellsPositive','BothEnds20','FourGenesObserved','FourMarginal5pct','FourMarginal10pct']
annotated=set()
for i,s in enumerate(['MacAM','MonoMac','MacWide']):
 q=flow[(flow.sender==s)&(flow.receiver=='CD4-IL7R')&(flow.group=='CWP')].set_index('stage').loc[stages];axes[1].plot(range(5),q.n,'o-',label=s,color=['#306c89','#b66537','#688259'][i])
 for x,y in enumerate(q.n):
  if (x,y) not in annotated:axes[1].annotate(str(y),(x,y),xytext=(4,7),textcoords='offset points',fontsize=8);annotated.add((x,y))
axes[1].set_xticks(range(5),['Paired\n>0 cells','Both ends\n>=20','4 genes\nobserved','4 marginal\n>=5%','4 marginal\n>=10%']);axes[1].set_ylim(-.5,10.5);axes[1].set_ylabel('CWP donors');axes[1].set_title('IL7R | Fixed-rule attrition',loc='left',fontweight='bold');axes[1].legend(frameon=False);fig.tight_layout();fig.savefig(F/'01_gap_IL15_measurement.png',dpi=170);plt.close(fig)
L='IL15 TNF CLSTN1 ADAM15 PLAU DSC2 CXCL2 ICAM2 VTN VEGFB GSTP1 RTN4'.split();summ=fin.groupby('ligand').agg(evaluated=('P','count'),nominal=('P',lambda x:sum(x<.05)),minP=('P','min'),minFDR=('FDR','min')).reindex(L).reset_index();wr(summ,'report_all_candidate_summary.tsv')
fig,axes=plt.subplots(1,2,figsize=(12,5.8),gridspec_kw={'width_ratios':[1.25,1]});axes[0].barh(L,summ.nominal,color='#357d9e');axes[0].invert_yaxis();axes[0].set_xlabel('Nominal P<0.05 specification rows');axes[0].set_title('All 12 candidates | Final model (0 BH<0.05)',loc='left',fontweight='bold')
for i,v in summ.iterrows():axes[0].text(v.nominal+.4,i,str(v.nominal),va='center')
steps=['Unadjusted','Group','Source','SourceSenderDepth','SourceBothDepth'];iu=up[(up.ligand=='IL15')&(up.readout=='AUC')];p=iu.groupby('step').agg(nominal=('P',lambda x:sum(x<.05)),computed=('P','count')).reindex(steps);axes[1].plot(range(5),p.nominal,'o-',color='#b66537');axes[1].set_xticks(range(5),['None','+Group','+Source','+Sender\ndepth','+Both\ndepths']);axes[1].set_ylabel('Nominal specification rows');axes[1].set_title('IL15 | Same-roster nested AUC models',loc='left',fontweight='bold');axes[1].set_ylim(-1,21)
for x,y in enumerate(p.nominal):axes[1].text(x,y+.7,str(y),ha='center')
fig.tight_layout();fig.savefig(F/'02_gap_candidate_nominal.png',dpi=170);plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(12,5.5));endpoints=['IL7R_nominal_original','IL7R_nominal_nonRiboMito','MatchedControl1','MatchedControl2','MatchedControl3'];labels=['Original','Non-RiboMito','Control 1','Control 2','Control 3']
for seed,col in zip([17,29,43],['#306c89','#b66537','#688259']):
 q=controls[(controls.variant=='DropOverlap')&(controls.seed==seed)].set_index('endpoint').reindex(endpoints);axes[0].plot(range(5),q.rho,'o-',label=str(seed),color=col)
 q=con[(con.variant=='DropOverlap')&(con.seed==seed)&(con.adjust=='SourceReceiverDepth')].set_index('endpoint').loc[['original','reference','excess']];axes[1].plot(range(3),q.rho,'o-',label=str(seed),color=col)
axes[0].set_xticks(range(5),labels,rotation=25,ha='right');axes[1].set_xticks(range(3),['Original','Mean controls','Original -\nmean controls'])
for ax in axes:ax.set_ylim(-1,0);ax.set_ylabel('Partial Spearman rho');ax.grid(axis='y',alpha=.2);ax.legend(title='Existing seed',frameon=False)
axes[0].set_title('IL7R | Disease program and matched backgrounds',loc='left',fontweight='bold');axes[1].set_title('IL7R | Fixed excess score remains nominal',loc='left',fontweight='bold');fig.tight_layout();fig.savefig(F/'03_gap_IL7R_specificity.png',dpi=170);plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(12,5.2));labels=[str(v.series)+'\n'+str(v.group) for _,v in comp.iterrows()];pct=comp.IKZF2_positive/comp.cells*100;axes[0].bar(labels,pct,color=[colors[x] for x in comp.group]);axes[0].set_ylabel('Pooled observed IKZF2-positive GZMA cells (%)');axes[0].set_title('GZMA | Source-concentrated rare detection',loc='left',fontweight='bold');axes[0].set_ylim(0,6)
for i,v in comp.iterrows():axes[0].text(i,pct.iloc[i]+.18,f'{v.IKZF2_positive}/{v.cells}\n{v.donors} donors',ha='center',fontsize=9)
q=gz[(gz.program=='TReg_minus_IKZF2')&gz.model.isin(['Base','Source','SourceDep'])].copy().reset_index(drop=True)
for i,v in q.iterrows():axes[1].errorbar(v.effect,i,xerr=[[v.effect-v.CI_low],[v.CI_high-v.effect]],fmt='o',capsize=3,color='#306c89' if v.roster=='PositiveCells' else '#b66537')
axes[1].set_yticks(range(len(q)),[f'{v.roster}/{v.model} (n={v.n})' for _,v in q.iterrows()]);axes[1].axvline(0,color='#999',lw=.8);axes[1].invert_yaxis();axes[1].set_xlabel('CWP - Control, score units (95% CI)');axes[1].set_title('TReg program excluding IKZF2 (7 genes)',loc='left',fontweight='bold');fig.tight_layout();fig.savefig(F/'04_gap_GZMA_composition.png',dpi=170);plt.close(fig)
gaps=[
 ['IL15测量','四基因在8/8个MacWide可配对病例均有RNA观察','固定10%覆盖不足；RNA不能证明蛋白复合物','原始映射、连续分布、固定门槛流失、深度分层、描述区间','真实捕获效率、蛋白转呈和接触需独立蛋白/空间证据'],
 ['候选→TF','最终172行名义关联；原BH全部未通过','规格行并非独立通路；IL15部分模型df不足','12候选12600模型、997名义行、CI/n/df、同名册阶梯','同来源独立供者/干预才能区分残余混杂和方向'],
 ['IL7R IKZF2→程序','去重叠/去CARD16/去核糖体及状态控制后仍负向','匹配背景也负向，疾病专属性受限','3组固定匹配、共同状态分层、同供者描述','共有背景之外仍有名义残余，需独立疾病程序验证'],
 ['IL7R固定残余','6版本P=0.00227–0.04014；96次删供者均不翻转','名义P受删供者影响；为选定后诊断','固定相减、背景协变量、16供者逐一剔除','无法用同一批供者证明疾病专属性或IKZF2因果抑制'],
 ['ICAM2/VEGFB同亚群链','原始程序与配体呈名义关联','固定残余关联P=0.234/0.292/0.299','相同供者原始/背景/残余直接比较','尚未连接到残余疾病程序；NicheNet对应活性亦不支持完整链'],
 ['GZMA IKZF2','G227病例16/369阳性；对照13/1891','G288病例0/292；来源仅病例且G227仅2病例','原始UMI、来源分解、少量标记共检出','缺G288对照及更多同来源病例，不能统计恢复该信息'],
 ['TReg模块','去IKZF2后≥20细胞名册Source仍名义支持','扩至所有阳性细胞名册Source不支持；去网络仅剩2基因','原始160背景精确重建；自包含/网络重叠分开','2基因读数及转录模块均不是Treg功能实验'],
 ['虚拟扰动','既有工具与counts可用','当前关联不能校准带方向的因果KO预测','已审查适用问题，未把未计算写成阴性','无本轮必需KO；后续若开展须独立固定功能终点和验证设计']]
ledger=pd.DataFrame(gaps,columns=['证据环节','已支持','断点','已完成补证','仍不可识别/下一步']);wr(ledger,'evidence_gap_ledger.tsv')
(O/'provenance/round_stop_assessment.json').write_text(json.dumps({'finite_round_completed':True,'followups_executed':['matched_background_specificity','common_state_standardization','fixed_reference_excess_contrast','fixed_excess_donor_influence'],'findings':'Broad shared covariation AND remaining nominal residual; selected upstream links do not connect residual.','remaining_not_identifiable':['capture versus biological absence','protein complex/trans-presentation/spatial contact','G288 case-control effect without controls','causal ligand-TF-program direction'],'experiments':'client decision'},ensure_ascii=False,indent=2))
def tab(d,cols=None):
 d=d[cols] if cols else d.drop(columns=[c for c in ['sample','samples','omitted','cell'] if c in d])
 return '<div class="table">'+d.to_html(index=False,border=0,escape=True,na_rep='未计算/NA',float_format=lambda x:f'{x:.5g}')+'</div>'
def image(n,c):return '<figure><img src="data:image/png;base64,'+base64.b64encode((F/n).read_bytes()).decode()+'"><figcaption>'+c+'</figcaption></figure>'
def section(h,b):return '<section><h2>'+h+'</h2>'+b+'</section>'
def fold(h,b):return '<details><summary>'+h+'</summary>'+b+'</details>'
parts=['<p class="eyebrow">HUMAN FIRST · FINITE GAP DIAGNOSTICS · 2026-10-05</p><h1>人类巨噬细胞→CD4证据缺口诊断</h1><p class="lead">IL7R的IKZF2–程序负相关包含广泛共同背景，也保留名义残余信号；目前仍未连接成巨噬配体→TF→疾病程序的完整因果链。</p><div class="note">新增有限诊断已经执行，不是待审批计划。全部正负及不可计算结果保存在osv01新目录；仅本报告与4张必要汇总图转移本地并交付Library。原始对象、单独供者/细胞表和复现ZIP保留服务器。旧报告不覆盖。</div>']
parts.append(section('1 · 本轮推进及边界','<p>连接osv01成功并先排查重复进程；复用原始RNA、既有供者聚合及3个seed的Original/DropOverlap AUC，未全量重聚类、未重跑SCENIC/Bootstrap或旧大网格。本轮在既有授权诊断基础上连续完成匹配背景、共同状态分层、固定残余以及逐供者影响检查。未按显著性改终点或搜索阈值，也未把虚拟KO作为必做任务。</p><p>RNA总体19供者：10 Control、9 CWP；G135=2 Control，G227=8 Control+2 CWP，G288=7 CWP。Primary包含CD4state，Core排除；既有AUC17名册不强加给RNA，IL7R本轮AUC可用16供者。CWP保持原疾病标签，不能自动等同其他尘肺分型。MacWide是MacAM/MonoMac union，不是独立重复；T/NK未替代客户的巨噬发送端。</p><p>前轮已完成全部受体/巨噬来源的疾病DE、通讯、NicheNet和鼠证据对照；本轮12候选为明确的缺口诊断名单，并不把12基因扩展成新发现全集。此前Core/IL7R五小面板30疾病比较无名义显著，以及CCR7 Inhib和GZMA TReg的名义提示，均保留原检验家族及范围，不能用配体相关替代疾病表型。</p>'))
parts.append(section('2 · IL15：RNA有观察，固定覆盖规则造成明显流失',image('01_gap_IL15_measurement.png','图1：左图每点为一个供者的观察阳性细胞比例，不显示供者标识；仅对应范围≥20细胞。Control=对照，CWP=病例；MacAM/MonoMac为复核巨噬来源，MacWide为union；IL7R为CD4-IL7R。右图为固定规则累计筛除：两端有细胞→两端≥20→四基因各有观察→四边际比例≥5%→≥10%；不是阈值优化。')+'<p>原始37643细胞核对、578项供者表达比较与既有聚合一致，最大比例误差&lt;5×10⁻¹⁶；IL15未出现基因计数或映射丢失。IL7R病例MacWide配对从9→8→8→2→0，Core从9→8→8→3→0。由此可见，10%不通过并不意味着四个基因完全未检出。四边际条件也不等于同细胞复合物。</p><p>深度四分位观察显示病例MacWide的IL15、IL15RA在高深度端检出比例中位数分别增加10.95、8.74个百分点；IL7R的IL2RB/IL2RG分别增加7.69、18.59个百分点。深度既包含捕获/复杂度，也可能关联生物状态，不能反推真实捕获效率。≥20细胞的相关测量中有9项零检出，其中3项的描述性Wilson上限仍≥10%；这些区间仅在细胞独立二项近似下描述采样，不是供者层面的生物推断或真实表达率置信区间。</p><p><a href="https://reactome.org/content/detail/R-HSA-8983432">Reactome IL15 signaling</a>及<a href="https://reactome.org/content/detail/R-HSA-449115">IL15:IL15RA与IL2RB/IL2RG结合事件</a>支持蛋白层面的转呈/顺式模式及β/γ链信号；该通路涉及JAK1/JAK3和STAT3/STAT5，不能自动替换为患者中已证明的STAT1。由机制与测量层级推断，固定10%单细胞RNA共检出并非证明蛋白复合物的必要或充分条件；现有RNA也无法判别患者究竟采用哪一模式。</p>'+fold('全部范围的固定规则流失',tab(flow))+fold('来源×疾病组连续测量汇总',tab(rd('IL15_source_measurement_summary.tsv')))))
parts.append(section('3 · 全部12候选→TF：保留名义结果与不可计算结果',image('02_gap_candidate_nominal.png','图2：左为最终Group+Source+双端深度模型的名义规格行数；每个配体有182项可计算RNA/AUC规格，不是独立生物重复。右为IL15 AUC固定供者名册阶梯，None/Group/Source/Sender depth/Both depths分别为未调整/疾病组/来源/发送深度/双端深度。')+'<p>共12600行，5步模型中997行P&lt;0.05；最终模型172行名义关联，0行通过其原BH家族。所有CI、n、病例/对照数、df与状态均保留。不能把“BH未通过”写成“没有任何关联”。RNA、3个网络seed和DropOverlap属于不同读数/规格，不是独立复制。最终模型与前轮相同名册的rho/P逐项复现。</p>'+tab(summ)+'<p>IL15在双端深度前有17条名义AUC规格，最终0条：其中6条CCR7–STAT1/MacAM在n=7时df从3降至2，属于不可计算，不能算阴性；其余11条保留供者，效应减弱且df下降。因此这里主要是条件关联和自由度变化，不是最后一步丢失了供者。完整五步模型固定同一名册，样本筛除与协变量校正严格分开。</p>'+fold('全部997条名义规格：方向、近似95%CI、P/BH、n/df',tab(nom,['receiver','tf','readout','variant','seed','sender','ligand','step','n','n_case','n_control','df','rho','CI_low','CI_high','P','FDR','status']))+fold('IL15全部五步规格（含未支持和不可计算）',tab(up[up.ligand=='IL15'],['receiver','tf','readout','variant','seed','sender','step','n','df','rho','CI_low','CI_high','P','FDR','status']))))
parts.append(section('4 · IL7R：稳定相关、共同背景和残余信号必须同时报告',image('03_gap_IL7R_specificity.png','图3：仅为清晰展示绘制DropOverlap的3个既有seed；完整表同时包含Original。左为原88基因程序、去7个核糖体基因后的81基因程序及3组各88基因固定匹配背景；右为原程序、3组背景均值、两者固定相减后的读数。纵轴为控制Group/Source/受体深度的偏Spearman相关，不是调控方向证据。')+'<p>去IKZF2网络重叠实际仅从88基因中去除CARD16；去RiboMito去除7个RPL/RPS基因；各自仍在3个DropOverlap seed保持名义负向。加入复杂度、线粒体比例、Naive或TAct单个协变量也未消除。共同Naive状态分箱等权后相关减弱但保留；同供者内部控制细胞深度/复杂度的描述性相关在15/16或16/16供者为负。这些不能把细胞数当独立样本量，亦不能排除相对表达/状态共变。</p><p>为检验专属性，用等供者权重参考的均值和检出比例、同RiboMito类别，确定3组互不重叠的88基因背景；264个匹配均完成。匹配距离中位数0.01789、最大1.376，匹配并不完美。匹配过程未使用AUC或疾病效应，但排除了既有疾病靶集，因此不是独立验证。18项背景–AUC规格全部呈名义负相关，说明稳定性并非疾病程序专属。</p><p><strong>固定扣除背景后，残余仍为负：</strong>6版本rho约−0.790至−0.598，P=0.00227–0.04014；加入背景为协变量也全部名义负向。96次逐供者剔除均未翻转，但各版本仅7–15/16次保持P&lt;0.05，最大P达0.291。因此当前证据支持“广泛共同背景＋名义残余”，不支持“全由背景解释”，也不足以认定残余为独立疾病机制。目标程序本就在同一批病例/对照中选出，选择依赖仍存在。</p>'+fold('全部匹配背景结果',tab(controls))+fold('全部固定残余与背景协变量结果',tab(con))+fold('固定残余逐供者影响汇总',tab(influence))+fold('全部去重叠及状态单变量诊断',tab(ind))+fold('共同状态分层标准化完整结果',tab(state))))
parts.append(section('5 · 同亚群链：上游尚未连接到残余疾病程序','<p>完整172条最终名义上游规格均在相同供者上补查TF–程序和配体–程序，共344行诊断。ICAM2/MonoMac、ICAM2/MacWide和VEGFB/MacWide在IL7R有原始程序名义关联，故执行预先记录的固定背景残余检查：三条原始程序P分别0.00778、0.02988、0.04734，残余程序P变为0.23370、0.29151、0.29890，且95%CI均跨0；不能称完整同亚群机制链。对应IL7R Source/Nominal-Up的NicheNet活动原本也不支持，不能隐藏与相关证据的分歧。</p>'+tab(chains)+'<p>TNF/MonoMac→IL7R IKZF2 RNA名义关联强（rho=−0.8886，P=0.000588，BH=0.0917），但同名册TF→程序P=0.0692、配体→程序P=0.246，仍有明确断点。预测配体活性、LR兼容性、真实疾病程序、同供者相关分别判断，不从其任一项推因果。</p>'+fold('所有名义上游的同名册程序检查',tab(rd('nominal_upstream_same_roster_program_checks.tsv')))))
parts.append(section('6 · GZMA和TReg模块：来源、少量细胞及自重叠',image('04_gap_GZMA_composition.png','图4：左为按来源/组汇总的GZMA观察IKZF2阳性比例，细胞分母只用于描述，不代替供者重复。右为固定原背景下去IKZF2的7基因TReg转录分数；PositiveCells=至少1细胞，Cells20=至少20细胞；SourceDep=来源＋深度；误差条为供者OLS 95%CI。')+'<p>GZMA共2554细胞仅29个IKZF2阳性：G227对照13/1891，G227病例16/369；G288六位具有GZMA细胞的病例共292细胞均未检出。G227病例仅2供者，来源校正的可分辨疾病对比主要依靠他们；G288缺对照的信息不能靠新增模型恢复。G227病例16个IKZF2阳性细胞中FOXP3/IL2RA/CTLA4/TIGIT分别共检出3/1/4/8，不能据此整体重标为Treg。</p><p>为排除模块自包含，固定原160背景基因重建TReg分数，去IKZF2后≥20细胞名册n=12的Source/SourceDep P=0.0210/0.0210；扩至17位有细胞供者则为0.207/0.171。WithinG227仍仅2病例（≥20时7对照）。去整个IKZF2网络后只剩IKZF4和LRRC32两基因，某些深度模型P较小也不能再称完整TReg功能或用其替换预定终点。</p><p>复现时发现原M2与NicheNet别名空间不同，5个固定背景符号在后者缺失/改名，初版重构误差0.005344导致断言停止。已返回原M2 canonical原始计数及完整160背景重建，逐细胞最大误差6.11×10⁻¹⁶、供者1.11×10⁻¹⁶；旧发布分数未改变。无效初版保留为INVALID，仅正确重构进入本报告。该修复不影响IL15原始测量核对。</p>'+tab(comp)+fold('全部TReg独立性模型及新诊断BH',tab(gz))))
parts.append(section('7 · 检验家族、有限规格与方法限制','<p>所有关联以供者为独立单位，部分秩相关控制固定协变量；样本数不因旧“每组n&lt;10”规则一律拒算。要求每组至少2、设计满秩且df≥3，不能计算者写明状态，FDR留NA。相关95%CI为Fisher近似，方差1/(df−1)，小样本下仅探索性；OLS分数效应亦为描述现有样本的条件估计。</p><ul><li>候选→TF：每readout/variant/seed/step跨全部12候选×3发送端×5受体×2TF作BH；仅可计算P计入，最终每版本m=312，所有家族标签随表保存。</li><li>TF→程序沿用每readout/variant/seed/step的完整预设范围。既有疾病DE/旧M1/M2和前轮BH不重定义、不因新小家族改变结论。</li><li>TReg去重叠诊断：每名册/模型3读数；IL7R去重叠/状态诊断：每variant/seed/adjust 4读数；匹配背景：每variant/seed 5终点；共同状态：每variant/seed/axis/mode 2读数。仅说明新有限诊断，不是对旧全程序家族的替代。</li><li>从名义上游选出的同名册链、固定背景残余和逐供者剔除均为选定后诊断，保留P/CI而不新增“发现BH”。seed及All/BG1敏感性不是独立验证。</li></ul><p>新规格按证据断点分轮冻结并留痕：GAP_SPEC、SPECIFICITY_FOLLOWUP_SPEC、CONTRAST_FOLLOWUP_SPEC、contrast_influence_spec。使用实际表达，不补造低检出值；5/10%为已有固定规则诊断，不搜索最显著阈值。</p>'))
parts.append(section('8 · 证据缺口与本轮停止边界',tab(ledger)+'<p>剩余问题是当前数据的可识别性，不是审批或算力阻断：G288没有对照；RNA没有蛋白转呈/接触观测；同一数据选出的程序和网络无法提供独立机制验证。继续在相同小样本叠加模型不能生成这些信息。本轮已完成有判别力的后续诊断后收束，不以未跑KO写成阴性。后续目的基因机制工作应先约定独立程序/功能读数及验证对象，实验由客户决定。</p><p>人鼠证据仍沿用前轮方向对照：GSE314943的3+3供者/文库条件、既有Core混合与ActCD4名义靶集分别保留；GSE183682的7/56天各组1文库仅方向。前轮human Source/All/FDR-Up同向共享靶基因计数为0，不能把跨物种未对齐的亚群或prior分数当独立重复。没有新增小鼠调控推断。</p>'))
parts.append(section('9 · 服务器交付、源码与复现','<p>目录：<code>'+str(O)+'</code>。结果含全12600模型、全部名义及未支持/不足状态、原始表达核对、靶集/背景完整成员、控制匹配、残余/影响诊断、证据缺口表、4图、中文报告、代码和SHA256清单。前轮human_first及human_alignment目录保留完整原始疾病/通讯/NicheNet/鼠对照资源版本与基因映射记录；本轮不改变prior。</p><p>复现包仅服务器：<code>deliverables/human_gap_diagnostics_20261005_reproducible_SERVER_ONLY.zip</code>。它收录本轮代码/表/日志/规格及必要前轮输入；原始大型计数和既有AUC对象保留原路径，以哈希清单定位。README明确原M2重构修复顺序。报告不含供者标识；完整供者/细胞明细保留服务器。</p>'))
css='body{font:16px/1.7 "Microsoft YaHei",Arial,sans-serif;color:#243442;background:#f3f5f6;margin:0}main{max-width:1180px;margin:36px auto;padding:48px;background:white}h1{font-size:36px;line-height:1.3;color:#183e54}h2{font-size:23px;color:#1c526c;margin-top:40px;border-bottom:1px solid #d5e2e6;padding-bottom:10px}.eyebrow{letter-spacing:2px;color:#548395;font-size:12px}.lead{font-size:21px}.note{padding:18px;background:#edf4f7;border-left:4px solid #347e9e}figure{margin:24px 0}img{width:100%;height:auto}figcaption{font-size:14px;color:#536974}table{border-collapse:collapse;font-size:12px;white-space:nowrap}th,td{padding:7px 10px;border-bottom:1px solid #dde6e9;text-align:left}th{background:#eaf1f4}.table{overflow:auto;max-height:650px;margin:15px 0}details{margin:20px 0;padding:12px;background:#f6f8f9}summary{cursor:pointer;font-weight:bold}code{overflow-wrap:anywhere;background:#eef2f4;font-size:13px}a{color:#126589}section{margin-bottom:30px}'
report=D/'human_gap_diagnostics_20261005_report_zh.html';report.write_text('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>人类通讯证据缺口诊断</title><style>'+css+'</style><main>'+''.join(parts)+'</main></html>',encoding='utf-8')
readme='''人类证据缺口有限诊断；2026-10-05；osv01
执行环境：Python batch_bench；R_4.5.2。仅按进程设置OPENBLAS_NUM_THREADS=1/OMP_NUM_THREADS=1。
原始目录变量P=/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA。
所有脚本采用固定P路径；复现须保持服务器原项目依赖。没有改变原始对象、全局配置或权限。
顺序：gap_export.R → gap_raw.py → gap_module_raw.py → gap_compute.py → gap_controls.R → gap_specificity.py → gap_contrast.py → gap_contrast_influence.py → gap_summary.py → gap_report.py。
关键：gap_export.R中的临时NicheNet空间TReg分数必须由gap_module_raw.py使用原M2 canonical计数和160背景覆盖后再做推断。首轮compute断言发现误差后停止；实际以gap_resume_independence.py恢复。修复后完整gap_compute.py可从头执行。INVALID_mapped_space_TReg_reconstruction为无效诊断，不是结果。
gap_report.py最终更新证据缺口表，覆盖gap_summary.py中的早期阶段措辞。
gap_compute.py读取前轮align_compute.py的已冻结OLS函数；gap_contrast等读取gap_compute.py的固定相关函数，未执行其主程序。
背景成员、靶集、模型全表、不可计算状态和日志在results/provenance；旧完整NicheNet结果和资源映射仍在human_first_20261004_v1，既有AUC在02_analysis。归档不复制大型原始对象。
仅HTML报告及4张汇总图允许转移本地/Library；本ZIP及单独供者/细胞表仅服务器。
'''
(D/'README_reproduction_zh.txt').write_text(readme,encoding='utf-8')
# Preserve small direct input tables and inherited source/resource manifests in the server-only archive.
dep=O/'dependencies';dep.mkdir(exist_ok=True)
for n in ['TF_readout_donors.tsv.gz','stage1_donor_design.tsv','alignment_gene_donors.tsv.gz','TF_excluded_program_donors.tsv.gz','fixed_TF_bridge_associations.tsv.gz','IL15_matched_donor_coverage.tsv.gz','existing_TF_edges_complete.tsv.gz']:
 p=V/n
 if p.exists():(dep/n).write_bytes(p.read_bytes())
sources=[]
for folder in [P/'02_analysis/mechanism_step2_20260927_154919_EsN4gk/input',P/'02_analysis/scenic3_20261001_123452_Kvp9h5/input']:
 for p in sorted(folder.glob('human*')):
  if p.is_file():sources.append(p)
for p in [P/'11_reports/human_first_20261004_v1/results/stage2_aggregates.rds',P/'02_analysis/mechanism_step2_20260927_154919_EsN4gk/tables/human_score_controls.tsv.gz']:
 if p.exists():sources.append(p)
for p in sorted(V.glob('*.tsv*')):sources.append(p)
manifest=[{'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sources]
(O/'provenance/source_sha256_manifest.json').write_text(json.dumps(manifest,indent=2))
localfiles=[p for folder in ['code','results','provenance','figures','dependencies'] for p in (O/folder).rglob('*') if p.is_file()]+[report,D/'README_reproduction_zh.txt']
audit=[{'path':str(p.relative_to(O)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(localfiles)]
(D/'SHA256_manifest.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2))
archive=D/'human_gap_diagnostics_20261005_reproducible_SERVER_ONLY.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in sorted(localfiles)+[D/'SHA256_manifest.json']:z.write(p,str(p.relative_to(O)))
with zipfile.ZipFile(archive) as z:assert z.testzip() is None
out={'status':'COMPLETE','report':str(report),'report_bytes':report.stat().st_size,'report_sha256':sha(report),'figures':[str(p) for p in sorted(F.glob('*.png'))],'archive':str(archive),'archive_bytes':archive.stat().st_size,'archive_sha256':sha(archive),'archived_files':len(localfiles)+1,'nominal_final':172,'BH_final':0,'upstream_rows':12600,'residual_influence_no_sign_flip':bool((influence.sign_flips==0).all()),'external_source_hashes':len(manifest)}
(D/'completion.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
