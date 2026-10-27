from pathlib import Path
import os,json,hashlib,base64,re,zipfile,html
os.environ['OPENBLAS_NUM_THREADS']='1'
import pandas as pd,numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/phenotype_candidate_discovery_20261005_v1';R=O/'results';F=O/'figures';D=O/'deliverables';H=P/'11_reports/human_first_20261004_v1/results'
rd=lambda n:pd.read_csv(R/n,sep='\t',low_memory=False)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
qa=json.loads((O/'provenance/final_audit_QA.json').read_text());assert qa['original_BH_reproduced'] and qa['program_reconstruction_verified']
ph=rd('all_12_program_phenotype_models.tsv');new=rd('complete_direct_phenotype_discovery_models.tsv.gz');sel=rd('shortlist_frozen.tsv');fm=rd('selected_fixed_roster_associations.tsv');inf=rd('selected_influence_summary.tsv');comp=rd('program_component_selected_associations.tsv');cd=rd('program_component_disease_models.tsv');pdisease=rd('selected_paired_roster_disease_models.tsv');cov=rd('selected_expression_coverage_summary.tsv');trace=rd('all_gene_endpoint_final_decisions.tsv.gz');old=rd('complete_231168_candidate_audit.tsv.gz');src=rd('selected_source_stratified_models.tsv');memb=rd('selected_complete_program_DE_members.tsv');weights=rd('selected_full_program_prior_weights.tsv');memb=memb.merge(weights,on=['ligand','gene'],how='left');memb.to_csv(R/'selected_complete_program_DE_prior_members.tsv',sep='\t',index=False,na_rep='NA')
leads=new.merge(sel,on=['ligand','receiver','metric','sender']);leads=leads[(leads.kind=='program')&leads.denominator.eq('within_scope')];leads.to_csv(R/'selected_discovery_models_original_FDR.tsv',sep='\t',index=False)
cards=[]
for rank,lig in enumerate(['HBEGF','ALCAM','HMGB1'],1):
 q=sel[sel.ligand==lig].iloc[0];a=leads[(leads.ligand==lig)&(leads.step=='Source')].iloc[0];d=ph[(ph.receiver==q.receiver)&(ph.metric==q.metric)&(ph['filter']=='All')&(ph.roster=='Cells20')&(ph.model=='Source')].iloc[0];raw=comp[(comp.ligand==lig)&(comp.component=='target_raw')&(comp.step=='Source')].iloc[0]
 cards.append(dict(review_priority=rank,gene=lig,role='MacAM sender/state candidate',receiver=q.receiver,phenotype=q.metric,phenotype_n=d.n,phenotype_CWP_effect=d.effect,phenotype_P=d.P,phenotype_FDR_expanded=d.FDR,paired_n=a.n,paired_controls=a.n_control,paired_cases=a.n_case,rho=a.rho,CI_low=a.CI_low,CI_high=a.CI_high,P=a.P,FDR_discovery=a.FDR,FDR_global=a.FDR_global,raw_target_rho=raw.rho,raw_target_P=raw.P,conclusion={'HBEGF':'conditional_source_adjusted_relative_program_hypothesis_receptor_unresolved','ALCAM':'contact_ligand_hypothesis_relative_score_background_sensitive','HMGB1':'lower_priority_relative_TReg_transcript_association_not_Treg_function'}[lig]))
cards=pd.DataFrame(cards);cards.to_csv(R/'three_independent_evidence_cards.tsv',sep='\t',index=False)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False});colors={'HBEGF':'#217d8e','ALCAM':'#b06c37','HMGB1':'#71599c'}
fig,ax=plt.subplots(2,2,figsize=(14.5,10.5));a=ax[0,0]
for i,(rec,metric) in enumerate([('CD4-IL7R','score_Mitophagy'),('CD4-GZMA','score_TReg')]):
 z=ph[(ph.receiver==rec)&(ph.metric==metric)&ph['filter'].eq('All')&ph.roster.eq('Cells20')&ph.model.eq('Source')].iloc[0];a.errorbar(z.effect,i,xerr=[[z.effect-z.CI_low],[z.CI_high-z.effect]],fmt='o',capsize=4,color='#217d8e');a.text(-.19,i+.22,f'n={z.n}; P={z.P:.4g}; expanded q={z.FDR:.3f}',fontsize=9)
a.set_yticks([0,1],['IL7R / Mito','GZMA / TRegProg']);a.invert_yaxis();a.axvline(0,color='#888',lw=.8);a.set_xlim(-.2,.13);a.set_ylim(1.55,-.5);a.set_title('A | Actual disease phenotype (Source)',loc='left',fontweight='bold');a.set_xlabel('CWP - Control effect (95% CI)')
a=ax[0,1]
for i,lig in enumerate(['HBEGF','ALCAM','HMGB1']):
 z=cards[cards.gene==lig].iloc[0];a.errorbar(z.rho,i,xerr=[[z.rho-z.CI_low],[z.CI_high-z.rho]],fmt='o',capsize=4,color=colors[lig]);a.text(-.98,i+.24,f'n={z.paired_n}; P={z.P:.4g}; discovery q={z.FDR_discovery:.3f}',fontsize=9)
a.set_yticks(range(3),['HBEGF / IL7R-Mito','ALCAM / IL7R-Mito','HMGB1 / GZMA-TRegProg']);a.set_ylim(2.65,-.6);a.axvline(0,color='#888',lw=.8);a.set_xlim(-1.03,1.03);a.set_title('B | MacAM - receiver association (Source)',loc='left',fontweight='bold');a.set_xlabel('Partial rank correlation (approx. 95% CI)')
a=ax[1,0];steps=['Unadjusted','Group','Source','SourceBothDepth']
for lig in ['HBEGF','ALCAM','HMGB1']:
 z=leads[leads.ligand==lig].set_index('step').loc[steps];a.plot(range(4),z.rho,'o-',label=lig,color=colors[lig])
a.axhline(0,color='#888',lw=.8);a.set_xticks(range(4),['Raw','Group','+Source','+BothDepth']);a.set_ylim(-1.05,1.05);a.set_ylabel('Partial rank correlation');a.set_title('C | Fixed donor roster; nested models',loc='left',fontweight='bold');a.legend(frameon=False,loc='lower left',fontsize=9)
a=ax[1,1];component=['original_score','target_raw','background_raw'];markers=['o','s','^']
for i,lig in enumerate(['HBEGF','ALCAM','HMGB1']):
 for j,c in enumerate(component):
  z=comp[(comp.ligand==lig)&comp.component.eq(c)&comp.step.eq('Source')].iloc[0];a.plot(z.rho,i+(j-1)*.18,markers[j],color=colors[lig],ms=8,label=['Original relative score','Raw target mean','Fixed background mean'][j] if i==0 else None)
a.axvline(0,color='#888',lw=.8);a.set_xlim(-1.03,1.03);a.set_yticks(range(3),['HBEGF','ALCAM','HMGB1']);a.set_ylim(2.7,-.6);a.set_title('D | Score decomposition (Source)',loc='left',fontweight='bold');a.set_xlabel('Partial rank correlation');a.legend(frameon=False,fontsize=8,loc='lower left',bbox_to_anchor=(0,-.27))
fig.suptitle('Phenotype-centered discovery | Independent hypotheses, same human dataset',x=.025,ha='left',fontsize=15,fontweight='bold');fig.tight_layout(rect=[0,.04,1,.96]);fig.savefig(F/'independent_phenotype_candidates_summary.png',dpi=145);plt.close(fig)
# Required raw/stratified plots with donor labels remain on server only.
donors=rd('selected_fixed_donor_inputs_SERVER_ONLY.tsv.gz');donors=donors[donors.followup_roster=='All_eligible'];sc={'G135':'#d9983c','G227':'#297e96','G288':'#9c5b90'}
fig,axes=plt.subplots(1,3,figsize=(16,5))
for a,lig in zip(axes,['HBEGF','ALCAM','HMGB1']):
 z=donors[donors.ligand==lig]
 for (series,group),v in z.groupby(['series','group']):
  a.scatter(np.log1p(v.x),v.y,c=sc[series],marker='o' if group=='Control' else '^',label=series+'/'+group)
  for rr in v.itertuples():a.annotate(rr.sample,(np.log1p(rr.x),rr.y),fontsize=6,xytext=(3,3),textcoords='offset points')
 a.set_title(lig+' / '+z.receiver.iloc[0]);a.set_xlabel('log1p(sender mean CP10k)');a.set_ylabel(z.metric.iloc[0]);a.legend(fontsize=6,frameon=False)
fig.tight_layout();fig.savefig(F/'selected_raw_donor_source_plots_SERVER_ONLY.png',dpi=160);plt.close(fig)
def tbl(d,cols=None):
 if cols:d=d[cols]
 d=d.drop(columns=['sample','samples','omitted','definition_source'],errors='ignore');return '<div class="tablewrap">'+d.to_html(index=False,na_rep='NA',border=0,float_format=lambda x:f'{x:.6g}')+'</div>'
def para(s):return '<p>'+s+'</p>'
def detail(title,d,cols=None):return '<details><summary>'+title+'</summary>'+tbl(d,cols)+'</details>'
img=base64.b64encode((F/'independent_phenotype_candidates_summary.png').read_bytes()).decode();count=json.loads((O/'provenance/compute_complete.json').read_text())
parts=['<!DOCTYPE html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>表型中心独立候选发现：完整审计与三方向证据卡</title><style>body{font-family:Segoe UI,Microsoft YaHei,sans-serif;max-width:1280px;margin:32px auto;padding:0 24px;color:#173548;line-height:1.75;background:#f7fafc}h1{font-size:30px}h2{margin-top:32px;border-bottom:2px solid #2b7c8c}h3{color:#216b7b}.callout{background:#e4f0f3;border-left:5px solid #2b7c8c;padding:18px}.limit{background:#fff1df;padding:14px}.tablewrap{overflow:auto;margin:14px 0}table{border-collapse:collapse;background:white;font-size:12px;width:100%}td,th{padding:7px;border-bottom:1px solid #d7e1e7;text-align:left;vertical-align:top}th{background:#e8f0f4}img{width:100%}figure{margin:20px 0;padding:10px;background:white}figcaption{font-size:13px;color:#506774}details{padding:12px;margin:14px 0;background:#eef4f7}summary{cursor:pointer;font-weight:600}code{overflow-wrap:anywhere;font-size:12px}a{color:#176b85}</style><h1>表型中心的独立候选发现</h1><p>2026-10-05 · 人类患者优先 · 全候选池审计 · 最多三个独立方向</p>',
'<div class="callout"><b>可交给客户的结论：</b>本轮找到三条新的、可被实验或独立患者数据检验的名义关联假说。HBEGF–MacAM/IL7R 相对 Mitophagy 程序的统计稳定性较好；ALCAM 的 CD6 接触配对解释较明确，但相对评分与背景的关系明显；HMGB1–GZMA/TReg 程序为较低优先级。三者各自独立，不要求连接 IKZF2、STAT1 或其他候选。它们尚不能写成已证实改变线粒体自噬、Treg 功能或导致疾病。</div>',
para('这不是将前15候选重新排序。完整原池包含688个先验配体、231168条候选规格和17596条同供者关联规格；新增直接表型审计覆盖全池，三方向均为此前未重点核验的基因。名义显著保留为探索关联，完整发现FDR同时披露。候选筛选及后续拆分均发生在同一数据集，不能作为独立重复验证。'),
'<figure><img alt="Four-panel independent phenotype candidate summary" src="data:image/png;base64,'+img+'"><figcaption>图1。A：原≥20细胞供者名单中的来源校正疾病效应及95%区间。B：具体MacAM发送端与接收亚群原定义相对程序的来源校正偏秩相关，区间为近似Fisher区间；q为完整发现家族FDR。C：固定同一供者名单的四层模型，Raw未调整，Group控制疾病组，Source另加数据来源，BothDepth另加发送及接收深度秩。D：原相对评分、未扣背景的目标均值与固定背景均值的分解，点不代表独立重复。Mito=原Mitophagy转录程序；TRegProg=原TReg转录程序，均非功能实验。IL7R/GZMA=CD4-IL7R/CD4-GZMA；MacAM=已复核巨噬群。所有标签为英文缩写，供者原始散点及来源分层图只保存在服务器。</figcaption></figure>',
'<h2>1. 样本、实际表型与本轮搜索规模</h2>',
para('RNA19为10 Control、9 CWP：G135为2对照；G227为8对照、2病例；G288为7病例、无对照。来源和疾病明显混杂；Source模型可计算不等于来源效应已被可靠分离。AUC17或前轮固定16不用于强行限制本轮RNA分析。全部五类CD4注释保留：CCR7 1172、GZMA 2554、IL7R 1379、Treg 600、CD4state 818个细胞。≥20细胞供者分别为9、12、16、10、7；CD4state仅1位对照，常规两组程序推断不足。Primary包含CD4state，Core排除它；比例分母一直单独标记。'),
para('原12个程序均保留：OXPHOS、MitoISR、MitoROS、Mitophagy、Inflam、PyroExec、ApoExec、TAct、TReg、Inhib、Naive、Cycle。新表型模型共1344项，其中1140项可计算、62项名义显著，旧504项P准确复现。12程序扩展与全局BH和原12程序小家族BH是不同检验家族，不能互换。原human_first只展开四功能程序，IL7R/Mitophagy是这次补齐的真实遗漏线索。'),
tbl(pd.DataFrame([{'层次':'原先验基因','规格数':688,'说明':'基因，不是独立验证数'},{'层次':'原候选表','规格数':231168,'说明':'含不支持、非正活性与数据不足'},{'层次':'原同供者关联','规格数':17596,'说明':'15520可计算；原BH已重现'},{'层次':'12程序疾病模型','规格数':1344,'说明':'1140有限P；含来源/深度/固定细胞名单'},{'层次':'固定表型端点','规格数':15,'说明':'5比例分母规格+10亚群程序'},{'层次':'新直接关联规格','规格数':82560,'说明':'688配体×15端点×2分离来源×4固定模型'},{'层次':'新直接关联有限P','规格数':44282,'说明':'未计算项保留原因；全局BH<0.05='+str(qa['new_global_BH05'])},{'层次':'基因—端点—来源审计','规格数':20640,'说明':'每行有保留/未深入理由'}])),
para('新直接关联的2个来源为MacAM与MonoMac。MonoMac是混合单核/巨噬状态，仅作来源上下文；MacWide是两者并集，不增加独立证据。原候选池三种发送范围仍完整保留。比例关联的接收细胞数门槛作用于明确的Primary/Core分母，不按某稀有亚群非零细胞数筛除零比例供者。'),
detail('15个冻结端点',rd('frozen_nominal_phenotype_endpoints.tsv')),
detail('全部亚群≥20细胞名单的名义疾病程序（包含All/BG1与不利结果）',ph[(ph.P<.05)&ph.roster.eq('Cells20')&~ph.receiver.isin(['Primary','Core'])],['receiver','filter','model','metric','n','n_case','n_control','effect','CI_low','CI_high','P','FDR','FDR_global']),
detail('全部组成疾病模型：分母不可混用',rd('all_composition_phenotypes.tsv'),['scope','denominator','model','n','effect','CI_low','CI_high','P','FDR','status']),
'<h2>2. 冻结策略与选择理由</h2>',
para('先于新候选关联计算保存DISCOVERY_FROZEN_SPEC.txt，依据所有已定义程序的任一名义疾病结果确定15端点。固定All/BG1、≥20细胞名单与四层相关模型；不追加检测阈值、seed或有利终点。可计算条件为至少2位每组、设计可识别；相关残差df≥3，疾病OLS df≥2。不执行旧每组n≥10禁算规则。所有名义阴性、未计算和BH不支持均保留。'),
para('三方向是完整多轴审计后的研究者判断，非预注册自动排名、非最小P冠军。ALCAM具有清晰接触配对和同向嵌套模型；HBEGF提供同一真实疾病程序的另一独立、来源条件下正向假说；HMGB1覆盖另一个有疾病信号的GZMA程序。NUCB2保留为候补状态相关线索，其ERAP1边不建立经典表面信号受体；NAMPT及其他新旧方向完整保留。稀疏受体不被判成生物不存在，也不单独作为剔除条件。三方向随后冻结，只做供者/来源/影响点与真实靶程序/背景两组核验，不因结果不理想换候选。'),
'<h2>3. 三个独立方向的证据卡</h2>',tbl(cards,['review_priority','gene','receiver','phenotype','phenotype_n','phenotype_CWP_effect','phenotype_P','phenotype_FDR_expanded','paired_n','paired_controls','paired_cases','rho','CI_low','CI_high','P','FDR_discovery','FDR_global','raw_target_rho','raw_target_P']),
'<h3>HBEGF：优先检验的来源条件下巨噬状态/发送假说</h3>',
para('角色：MacAM表达的候选发送因子或巨噬状态标记；接收端CD4-IL7R；终点为原11基因Mitophagy均值减固定背景。RNA16中CWP效应−0.094864、Source P=0.008151；同11人配对名单中疾病效应−0.083394、P=0.003893。发送端—程序Source ρ=0.847658、P=0.007860、完整发现FDR=0.406863；SourceBothDepth ρ=0.873744、P=0.022905、FDR=0.789354。未调整及仅疾病组模型P约0.20，说明关联依赖来源条件。'),
para('11次删一供者的Source相关均为正、P范围0.00180–0.0370；G227七人、两病例控制组别后ρ=0.915494、P=0.010410，仍是同一数据且样本小。BG1同九人Source ρ=0.992442、P约0.000085；这是筛选后的敏感性，不能用它的后选择小家族FDR替代完整发现FDR。经典EGFR/ERBB4在接收细胞检测稀疏；旧10%门槛的支持来自CD44，其原始证据类型实际为ppi_prediction。不能因此写成已建立HBEGF→CD44信号通路。'),
para('原程序未扣背景目标均值也随疾病下降，P=0.015429；HBEGF与该均值Source ρ=0.556625、P=0.151865，深度模型ρ=0.787897、P=0.062711，不能称为目标均值显著关联。可回答客户的问题是：已复核巨噬细胞的HBEGF状态是否标记IL7R细胞疾病相关的相对自噬转录状态，直接作用还是共同微环境变化仍待区分。HBEGF自身Source DE P=0.333、FDR约0.994，不作为自动否决，也不宣称病变由其差异表达驱动。'),
'<h3>ALCAM：受体解释较清楚，但背景依赖更明显</h3>',
para('角色：MacAM接触相关候选；配对CD6在IL7R病例9/9供者观察到RNA、6/9达到既有10%门槛。ALCAM–CD6原先验来自KEGG CAMs，并有直接细胞黏附原始研究背景。Source ρ=−0.795132、P=0.018329、完整发现FDR=0.597377；SourceBothDepth ρ=−0.889515、P=0.017636、FDR=0.727162。控制疾病组模型P=0.003336，未调整P=0.160，不能将两者合并写成无条件关联。'),
para('11次删一供者Source ρ仍为−0.858至−0.686，但P最高0.0889。All同九人Source P=0.0510；BG1同九人P=0.2306，说明筛选会减弱证据；BG1深度模型df不足不计算。未扣背景的11基因均值Source ρ=0.042583、P=0.920254；固定背景均值ρ=0.697372、P=0.054515。该结果支持相对程序关联，尚不能证明ALCAM特异地控制Mitophagy基因。'),
para('可检验假说：在MacAM–IL7R接触环境中，ALCAM/CD6标记的状态与IL7R相对Mitophagy转录读出相关。受体背景清晰使其具有验证价值，但不能据接触配对补上程序特异性。MacAM ALCAM在CWP方向为降低（Base log2FC−1.201、P=0.0147；Source−0.973、P=0.0510），而其供者内关联为负，这不能直接拼接成“下降ALCAM导致下降Mitophagy”的中介链。'),
'<h3>HMGB1：GZMA相对TReg转录程序的较低优先级线索</h3>',
para('角色：MacAM的RNA状态/候选胞外信号来源；接收端CD4-GZMA，终点为原TReg转录程序，不能将该亚群直接改名为功能性Treg。RNA12中Source疾病效应+0.038841、P=0.013119；配对9人中+0.038945、P=0.038197。Source ρ=0.792202、P=0.033684、完整发现FDR=0.765531；深度模型P=0.158074。'),
para('逐供者删除Source方向9/9一致，但P可到0.134；深度模型删除后仅1/9可计算，剩余为df不足。All/BG1相关秩相同不等于独立重复。未扣背景TReg目标均值的疾病Source P=0.465611，与HMGB1的Source P=0.881532。经典AGER/TLR4及LY96覆盖不足；支持旧门槛的ITGA4/FAS边不建立完整HMGB1信号。RNA不等于胞外释放、氧化还原形式或蛋白活性。该方向保留为相对转录程序的名义关联，低于前两方向，不能给出Treg抑制功能增强结论。'),
'<h2>4. 固定供者、来源与测量边界</h2>',tbl(rd('selected_paired_source_cell_summary.tsv')),
para('MacAM总9771细胞来自复核C1qAM、AM_Dom和AM_MARCO标签。病例687、对照9084个细胞高度不均衡；不将细胞作为独立样本。≥20细胞的MacAM供者为9对照/5病例。主关联只按两端合格供者独立赋权。AM_Dom含CTL068/CTL070主导簇；原簇身份及供者主导表原样保存。三个基因表达阳性细胞多有成熟巨噬标记，但单核标记可共同存在，这不是自动纯度证明。没有将排除的单核群替代MacAM。'),
detail('所有固定供者/筛选关联规格',fm,['ligand','roster','step','n','n_case','n_control','df','rho','CI_low','CI_high','P','FDR_followup_specs','status']),
detail('来源分层：没有对照的来源不伪造疾病推断',src,['ligand','receiver','source','step','n','n_case','n_control','rho','P','FDR','status']),
detail('全部删一供者影响摘要',inf),
detail('受体检测：零或低检出是测量边界',cov[cov.key.isin(['CD4-IL7R_All','CD4-GZMA_All'])&cov.gene.isin(['CD6','CD44','EGFR','ERBB4','AGER','TLR4','LY96','ITGA4','FAS'])]),
detail('实际原先验配对来源',rd('selected_LR_prior_original_provenance.tsv')),
detail('MacAM具体簇表达汇总，低细胞簇单列',rd('selected_cluster_expression_summary.tsv').query("scope=='MacAM'")),
detail('表达阳性细胞身份标记共现，不重新注释',rd('selected_positive_identity_summary.tsv')),
'<h2>5. 实际目标、背景和NicheNet支持到哪里</h2>',
para('原程序使用供者等权参考表达、25个表达分位箱、每个目标匹配20个固定控制基因并取并集，再计算目标平均log1pCP10k减控制均值。此次原始规范符号空间重建每细胞评分，最大误差<1e−10。完整目标与控制成员留服务器；不能以NicheNet别名映射后的背景替换原评分控制。'),tbl(rd('original_program_reconstruction.tsv')),
para('IL7R/Mitophagy的11原成员均保留，Source接收者DE可测试成员为SQSTM1、OPTN、MAP1LC3B：log2FC分别约−1.460、−1.311、−1.117，P=0.0470、0.0524、0.0904，基因FDR均>0.62；只有SQSTM1进入原名义下调靶集。GZMA/TReg的Source实际背景内成员更有限，不能因NicheNet对宽泛疾病集合有活性就声称预测到了完整TReg程序。'),
tbl(memb[(memb['filter']=='All')&memb.model.eq('Source')],['ligand','program','gene','logFC','PValue','FDR','in_original_background','Nominal05_Up','Nominal05_Down','weight']),
para('完整未删选配体—疾病靶成员、背景成员、所有原模型活性和相同供者关联均在服务器。NicheNet权重无正负干预含义；活性排行没有推断P/FDR。HBEGF/ALCAM在IL7R Source名义Down集合的活性为正，但对该宽泛集合的原同供者关联P约0.763/0.965；不能挪用相对Mitophagy程序关联替它补证据。'),
detail('评分拆分：全部疾病模型',cd,['receiver','program','filter','component','model','n','effect','CI_low','CI_high','P','FDR']),
detail('评分拆分：全部候选关联模型',comp,['ligand','receiver','component','step','n','rho','CI_low','CI_high','P','FDR']),
'<h2>6. 其他候选、旧排名及明确不支持</h2>',
para('TGFB1/SPP1只作审计示例，没有被强制入选：TGFB1对GZMA/TReg的Source直接关联P=0.877，对IL7R/Mitophagy P=0.816；SPP1分别为0.720、0.311。它们的较高NicheNet排名不能替代实际表型关联。SPP1在MacAM达到旧10%门槛的病例只有1位；此次仍计算可用表达关联并保留数据，不因旧LR门槛直接封闭候选。NUCB2的IL7R/Mitophagy相关为名义线索（Source P=0.00592），但ERAP1边与受体功能解释、独立表型验证均有缺口，本轮不添加为第四个深入方向。'),
detail('示例候选在所有已进入程序端点的结果',trace[trace.ligand.isin(['TGFB1','SPP1','NUCB2'])&trace.kind.eq('program')&trace.sender.eq('MacAM')],['receiver','metric','ligand','n_Source','rho_Source','P_Source','FDR_Source','rho_SourceBothDepth','P_SourceBothDepth','FDR_SourceBothDepth','observed_donors','median_sender_detection','shortlist_status']),
para('全体688基因摘要、20640个基因—端点—来源决定及231168条逐规格证据表均已保存。每个候选保留原正/非正活性、LR边界、实际目标成员、同供者P/FDR和本轮不深入理由；未深入不是证明其生物学不存在。319条新MacAM名义线索是依赖规格中的发现记录，不能称为319个可靠靶点。大量稀疏表达或零值并列的名义结果没有被用来选择“冠军”。'),
'<h3>全池复核补充：SERPINA1强混合来源信号与TXLNA精确诊断</h3>',
para('本轮全局复核发现SERPINA1与IL7R/Mitophagy在MonoMac中Source ρ=0.967027、P=2.904×10⁻⁷，原全局FDR=0.004287。不能仅凭混合标签把它当成单核细胞全部排除。为避免身份偷换，在结果出来前追加冻结审计说明，只检查原固定簇3/7，不重注释、不改变端点。簇3的配对8人Source ρ=0.855217、P=0.029926，但成熟巨噬和单核标记广泛共存；簇7仅6人、4病例2对照，Group P=0.223，Source df不足。MacAM对照Source ρ=−0.056625、P=0.894055。由此保留<strong>强混合髓系状态关联假说</strong>，明确其不是已确认的纯MacAM信号，也不是阴性；本轮三个MacAM详细证据卡不因此被偷偷替换。'),
tbl(rd('SERPINA1_source_identity_associations.tsv'),['sender','step','n','n_case','n_control','df','rho','CI_low','CI_high','P','FDR','status']),
detail('SERPINA1原完整发现家族结果',rd('SERPINA1_original_discovery_context.tsv'),['sender','step','n','rho','P','FDR','FDR_global']),
detail('SERPINA1固定来源细胞身份：标记共现不等于纯度',rd('SERPINA1_source_identity_summary.tsv')),
para('TXLNA–Treg/OXPHOS在6人中完全负秩相关，近似t公式产生10⁻²³至10⁻³¹的数值P，不能当作可靠的小样本精度。全部720种排列的未调整精确P=0.002778；48种组内排列的组别条件精确P=0.020833。二者仍为名义关联，但放回完整发现家族后不通过FDR；原近似P/FDR保存并标记，不隐藏。全局原3条数值BH通过记录经该诊断修正后剩1条，即SERPINA1混合来源；其全局校正FDR约0.01286。此排列为有限审计诊断，依赖固定供者可交换性，不是独立验证。'),
tbl(rd('TXLNA_exact_small_sample_diagnostic.tsv'),['ligand','receiver','metric','step','n','baseline_rho','permutations','as_or_more_extreme','exact_P']),
tbl(rd('all_global_hits_exact_diagnostic_comparison.tsv'),['ligand','receiver','sender','step','n','rho','P','P_exact_corrected','FDR_global','FDR_global_exact_corrected']),
'<h2>7. 既有小鼠证据：方向上下文，不冒充验证</h2>',
para('只复用既有GSE314943 3+3：ActCD4中Mitophagy原相对程序差异为−0.038711（All），与人IL7R方向相同；BG1为−0.046397。小鼠Treg则可能反向，保留全部结果。ActCD4与人IL7R没有建立一一亚群等价，因此同向只是上下文。人v1与鼠v2 NicheNet prior不同，不能比较名次绝对值；Core混合与ActCD4分开、FDR靶集与名义靶集分开。未找到可用的某候选鼠活性行不等于阴性。GSE183682每时间点每组1文库，仅方向、不生成P。'),
detail('全部既有鼠程序方向，包含不一致和低细胞数结果',rd('selected_mouse_program_directions.tsv')),
detail('三候选已有鼠NicheNet规格',rd('selected_mouse_existing_activities.tsv')),
'<h2>8. P、BH家族与发现偏差</h2>',
para('原17596关联按receiver/filter/DE-model/target-test/direction内所有有限配体与发送源计算BH，本轮准确重现。原疾病基因DE FDR按各接收者/filter/model全部测试基因计算。新1344程序模型BH按filter/roster/model跨所有scope×12程序，另报告全局1140有限P的BH。新82560直接关联按kind/receiver/metric/denominator/step内全部配体×2发送源的有限P计算BH，另报告全局44282有限P的BH；入选三个基因不重置发现家族。'),
para('后续固定模型的小家族FDR只描述后选择敏感性，不能取代发现FDR，也不能当作确认性检验。重复BG1别名不重复计入followup_specs家族；同一配对疾病端点由ALCAM/HBEGF共享时按不同端点去重。拆分疾病测试每filter/model为6个固定端点（2程序×3成分），拆分关联每模型9项（3方向×3成分）。这些较小FDR不升级全池发现结论。来源分层BH按固定模型内有限方向/来源计算；删一供者为描述性诊断，不选择最有利删除。'),
para('名义结果是可报告的关联，不因FDR未过线而改写为“没有发现”；同时，按疾病端点再筛候选的两阶段发现偏差不能由普通BH完全消除。当前属于探索，不是严格控制整个自适应流程错误率的确认研究。未计算写NA及原因，不能写阴性。'),
detail('本轮所有新检验家族、有限P数量及名义/FDR计数',rd('all_new_BH_families.tsv')),
'<h2>9. 已补缺口、剩余问题与客户可讨论方向</h2>',
para('已补齐此前四程序范围遗漏、全池而非top-ranking审计、候选到实际表型的直接关联、固定名单/来源/深度、具体巨噬身份、真实靶成员与评分背景拆分。没有用“因果链不完整”自动淘汰所有名义方向，也没有用重复同数据分析冒充独立证据。HBEGF与ALCAM可各自提出具体可检验假说；HMGB1的程序/受体限制更大，应后置。'),
para('最有价值的下一步科学问题不同：HBEGF需要区分巨噬状态标记与直接接收信号；ALCAM需要区分接触状态及背景相关性与Mitophagy目标程序本身；HMGB1需要验证胞外蛋白与真实接收/功能，而不只是RNA相对评分。独立来源平衡患者、蛋白定位/接触及功能干预可回答现数据无法识别的问题；具体实验由客户决定。此处不提供成功的因果结论，也不要求生信先完成所有实验才能把关联假说交给客户。'),
para('本轮不继续增加有利阈值、模型、种子或替换新基因。现有公开BAL患者资料可用于后续定位，但此前受者注释/病例对照和配对不足，不能直接称为本三方向独立复现；本轮没有为凑独立证据下载新大队列。新增大队列与功能机制整合应先说明科学问题与资源，原目的基因阶段仍单独讨论。'),
'<h2>10. 来源、复现与交付</h2>',
para('服务器：<code>'+str(O)+'</code>。完整表、独立供者明细、原始散点、代码、冻结规格、输入哈希、评分定义、R环境信息与复现ZIP均留服务器。此次只有本报告和图1转本地/Library。新结果不覆盖旧目录。复现入口：DISCOVERY_FROZEN_SPEC.txt → discover_compute.py → DISCOVERY_SHORTLIST_SPEC.txt → discover_followup.py / discover_raw_context.R → discover_decompose.py → discover_finalize_audit.py → discover_report.py。完整原矩阵仍按manifest路径复用，ZIP不重复打包大矩阵。'),
para('原始资料：human_first_20261004_v1（候选、DE、靶集、来源）、technical_delivery_20261002（原504模型及评分）、mechanism_step2（原程序与固定背景）、mechanism_step5及human_first stage3（既有鼠NicheNet与正交映射）。已安装human NicheNet v1 prior，原Zenodo 3260758；nichenetr 2.2.1.1/edgeR 4.8.2源结果复用。实际prior文件和所有读取源表有SHA256记录。本轮无新NicheNet阈值、全SCENIC、Bootstrap或重聚类。'),
para('受体解释仅作基础生物学核对：<a href="https://pubmed.ncbi.nlm.nih.gov/16352806/">Zimmerman等：CD6/ALCAM细胞黏附原始研究</a>；<a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC2199889/">Bennett等：特定CD44变体呈递肝素结合生长因子</a>（本数据没有证明CD44变体或完整信号）；<a href="https://pubmed.ncbi.nlm.nih.gov/25559892/">Yang等：MD-2参与二硫键HMGB1的TLR4信号</a>（本数据没有测量该蛋白状态）。文献不是本患者候选效果的独立验证。'),
detail('完整冻结发现规格',pd.DataFrame({'specification':(O/'provenance/DISCOVERY_FROZEN_SPEC.txt').read_text().splitlines()})),
detail('完整三方向后续规格',pd.DataFrame({'specification':(O/'provenance/DISCOVERY_SHORTLIST_SPEC.txt').read_text().splitlines()})),
detail('全池末轮统计/身份审计的显式范围修正',pd.DataFrame({'specification':(O/'provenance/DISCOVERY_AUDIT_ADDENDUM.txt').read_text().splitlines()})),
detail('全部1344预定义程序模型（供者身份不随报告输出）',ph,['receiver','filter','roster','model','metric','n','n_control','n_case','effect','CI_low','CI_high','P','FDR','FDR_global','status']),'</html>']
content='\n'.join(parts);content=re.sub(r'>\s+<','><',content);content=re.sub(r'</(?:td|th|tr)>','',content).replace('<tbody>','').replace('</tbody>','');report=D/'phenotype_independent_candidates_20261005_zh.html';report.write_text(content,encoding='utf-8')
assert not re.search(r'CTL\d|CWP\d',content.replace('CTL068/CTL070',''))
gap={'status':'COMPLETE_FINITE_DISCOVERY','directions':cards.to_dict('records'),'completed_gaps':['full_688_gene_pool','12_program_scope','fixed_roster_and_source','real_target_vs_fixed_background','reviewed_macrophage_identity','original_LR_provenance','existing_mouse_direction_context'],'remaining':['independent source-balanced patient replication','direct receptor signaling for HBEGF/HMGB1','protein/secretion/contact and functional endpoint validation','relative-score specificity for ALCAM/HMGB1'],'further_same_data_grid':'not_justified','experiments':'client_decision'};(O/'provenance/gap_closure_assessment.json').write_text(json.dumps(gap,indent=2,ensure_ascii=False))
manifest=pd.DataFrame([dict(path=str(p.relative_to(O)),bytes=p.stat().st_size,sha256=sha(p)) for folder in ['code','results','provenance','figures'] for p in sorted((O/folder).glob('*')) if p.is_file() and p.suffix!='.lock']);manifest.to_csv(O/'provenance/output_manifest.tsv',sep='\t',index=False)
archive=D/'phenotype_candidate_discovery_SERVER_ONLY.zip'
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as zz:
 for folder in ['code','results','provenance','figures']:
  for f in sorted((O/folder).glob('*')):
   if f.is_file() and f.suffix!='.lock':zz.write(f,str(f.relative_to(O)))
 zz.write(report,str(report.relative_to(O)))
with zipfile.ZipFile(archive) as zz:assert zz.testzip() is None;countzip=len(zz.namelist())
receipt={'status':'SERVER_COMPLETE_LIBRARY_PENDING','report':{'path':str(report),'bytes':report.stat().st_size,'sha256':sha(report)},'figure':{'path':str(F/'independent_phenotype_candidates_summary.png'),'bytes':(F/'independent_phenotype_candidates_summary.png').stat().st_size,'sha256':sha(F/'independent_phenotype_candidates_summary.png')},'server_only_zip':{'path':str(archive),'bytes':archive.stat().st_size,'sha256':sha(archive),'files':countzip,'testzip':'PASS'},'QA':qa};(D/'SERVER_DELIVERY_RECEIPT.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt,indent=2));print('REPORT_TABLES',content.count('<table'),'FIGURES',content.count('<img'))
