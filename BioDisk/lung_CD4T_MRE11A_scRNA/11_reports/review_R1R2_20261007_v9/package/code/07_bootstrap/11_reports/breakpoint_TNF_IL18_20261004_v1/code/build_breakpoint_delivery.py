from pathlib import Path
import json, hashlib, html, shutil, zipfile, datetime, sys
import pandas as pd
import numpy as np

O=Path(sys.argv[1]) if len(sys.argv)>1 else Path('breakpoint_delivery')
R=O/'results'; P=O/'provenance'
def rd(n): return pd.read_csv(R/n,sep='\t',low_memory=False)
m=rd('aligned_chain_models.tsv'); matrix=rd('chain_breakpoint_matrix.tsv')
roster=rd('donor_roster_reasons.tsv'); cov=rd('chain_detection_summary.tsv')
inf=rd('influence_summary.tsv'); sign=rd('same_donor_chain_sign_audit.tsv')
verify=json.loads((P/'verification.json').read_text()); loo=json.loads((P/'missing_LOO_complete.json').read_text())
before=pd.read_csv(P/'aligned_models_before_LOO.tsv',sep='\t')
assert len(m)==1368 and len(matrix)==12 and verify['status']=='PASS'
assert loo['historical_models_with_LOO_added']==12 and loo['omission_diagnostics']==193 and loo['remaining_missing_LOO']==0
assert np.allclose(m[['rho','P','FDR']],before[['rho','P','FDR']],rtol=0,atol=0)
assert not m.FDR.lt(.05).any()
FINAL='G+S+receiver_library+sender_library'
core=m[(m.policy=='aligned_steps')&(m.covariates==FINAL)]
eff=m.edge.eq('ligand_TF')|m.y.eq('State_Effector')
am=m.sender.eq('MacAM'); auc=m.modality.eq('AUC17')
def fmt(x):
 if pd.isna(x): return '—'
 if isinstance(x,(float,np.floating)): return f'{x:.6f}' if abs(x)<.01 and x!=0 else f'{x:.4f}'
 return str(x)
def table(d,cols=None):
 d=d[cols].copy() if cols else d.copy()
 return '<div class="scroll">'+d.map(fmt).to_html(index=False,escape=True,border=0,classes='data')+'</div>'
parts=[]; toc=[]
def sec(title,body):
 i='s'+str(len(toc)+1);toc.append((i,title));parts.append(f'<section id="{i}"><h2>{html.escape(title)}</h2>{body}</section>')
def p(s): return '<p>'+s+'</p>'
def detail(title,body):return '<details><summary>'+html.escape(title)+'</summary>'+body+'</details>'
def stats(d,extra=None):
 return table(d,['spec_id','chain','modality','seed','policy','covariates','edge','y','n','rho','CI95_low','CI95_high','P','FDR','df','LOO_valid','LOO_min','LOO_max']+(extra or []))

sec('1. 本次结论与证据边界',p('<b>两条候选主线保留，证据仍分段且依赖分析条件；本次结果不足以宣布完整因果链成立，也不足以据非显著结果更换 IKZF2 或 STAT1。</b>固定患者巨噬细胞→CD4 表型主线，仅研究 TNF–STAT1、IL18–IKZF2。MacAM 为主，MonoMac 与 MacWide 为并列背景分析；Core 与 CD4-IL7R 分开。')+
 '<div class="cards"><article><b>IL18–IKZF2 / MacAM–IL7R</b><p>同17人主校正：IL18–IKZF2 正相关，IKZF2–效应负相关，IL18–效应方向为负但 P=0.149595。方向相容；第三段尚无名义显著支持。</p></article><article><b>TNF–STAT1 / MacAM–Core</b><p>主校正下 TNF–STAT1、TNF–效应名义显著，STAT1–效应 P=0.067204。加入 Treg 比例的预设敏感性模型中，RNA 与 AUC 三段均名义显著；不能据此选择性替换主模型。</p></article><article><b>共同限制</b><p>当前两组 BH 家族内无 FDR&lt;0.05。受体 RNA 稀疏、来源不平衡、校正依赖及 seed 差异仍在。活性配体和功能性受体未测。</p></article></div>'+
 p('“断点”指证据链中尚未被支持或未被测量的环节，不表示生物通路已经被证伪。P&gt;0.05 不等于基因无效；宽区间表示精度有限，不能反过来证明真实链存在或“只是低功效”。'))

sec('2. 固定设计、复用与统计解释',p('冻结设计先于缺失规格拟合：3发送端×2接收端×2轴＝12链，共1368规格；850个缺失数值签名新拟合，482条历史规格复用，36条规格复用本轮已有数值签名。未重跑历史8264模型网格。最终报告生成只读取结果。')+
 p('TF 表达 RNA 与 SCENIC DropOverlap AUC 为不同测量。AUC 主 seed=17，29/43 仅敏感性；seed 不是独立生物重复。主要表型为冻结 Effector，Activation 平行保留。每个链/模态的主阶梯使用完全一致的供体名单，端点去重叠规则不变。')+
 p('阶梯：无校正→G（CWP/Control）→G+S（研究来源分类）→+receiver_library（接收端深度）→+sender_library（发送端深度）。额外 Treg_fraction 仅 Core，TNK_IFNG 两接收端，IFN_state_NoRegUnion 仅 TNF 的配体–TF/配体–表型，为相互独立的预设敏感性模型，并非全部协变量穷举。≥20细胞仅既定敏感性，不搜索新阈值。')+
 p('rho 为秩变换后、回归去除指定协变量所得残差相关；无校正时为 Spearman。P 采用探索性 t 近似；df=n−设计矩阵秩−1。95% CI 为 Fisher z 近似，标准误1/√(df−1)，非 bootstrap 或因果效应区间。小样本、并列秩和协变量估计下应审慎解释。n≥10且df≥3才做主推断。LOO 逐人遗漏同时重建秩及设计；范围表示影响诊断，不作显著性投票。')+
 p('本次 BH：primary=MacAM aligned_steps，260个唯一签名；sensitivity=其余预设模型，954个唯一签名。家族内按数值签名去重；同一签名可能跨家族，FDR可能不同。旧 FDR 保留为 legacy_FDR，不与本次 FDR 混用。主 family 包含两模态、两表型和近端代理，未仅挑选三条边。来源内可推断结果另作描述性家族。')+
 table(m.groupby(['FDR_family','policy']).agg(specifications=('spec_id','size'),nominal_P_lt05=('nominal_P_lt05','sum'),FDR_lt05=('FDR_lt05','sum')).reset_index()))

rr=roster[roster.chain.eq('MacAM/CD4-IL7R/IL18-IKZF2')]
base=rr[rr.roster.eq('RNA19')][['sample','group','series','sender_cells','receiver_cells']].copy()
for name in ['AUC17','AUC_cells20']:
 q=rr[rr.roster.eq(name)].set_index('sample');base[name]=base['sample'].map(q.included);base[name+'_exclusion_reason']=base['sample'].map(q.reason)
base.to_csv(R/'MacAM_IL7R_19_17_11_trace.tsv',sep='\t',index=False,na_rep='NA')
sec('3. 19→17→11：供体损失逐人可追溯',p('RNA完整链19人（10对照、9病例）；AUC17人（8对照、9病例）。CTL070、CTL103 无冻结 AUC：上游 Core 仅3、18细胞，CD4-IL7R分别3、5细胞，未进行 AUC 插补。')+
 p('MacAM–IL7R 在 AUC17 上要求发送端和接收端均≥20细胞后剩11人：再去除 CTL092（发送4）、CTL106（接收10）、CWP127（发送6）、CWP149（发送3）、CWP89（发送16）、CWP94（发送6）。保留 CTL068、CTL094、CTL098、CTL108、CTL111、CTL117、CWP109、CWP144、CWP152、CWP85、CWP86。不能把这11人与17人结果当作只改变协变量的直接比较。')+
 table(base)+p('Core 的 MacAM ≥20敏感性为12人（多保留 CTL106）；MonoMac两接收端均14人；MacWide Core16人、IL7R15人。完整名单如下，详细原因见 donor_roster_reasons.tsv。')+
 detail('全部12链名单与排除名单',table(rd('chain_roster_summary.tsv')))+
 p('来源与病例并不平衡：RNA G135=2对照；G227=8对照+2病例；G288=7病例。AUC G135仅CTL068，G227=7对照+2病例，G288=7病例。病例/来源联合校正并不能创造缺失的组内对照。'))

overview=[]
for _,r in matrix.iterrows():
 d={'链':r.chain,'n':r.n}
 for e,label in [('LT','配体–TF'),('TFP','TF–效应'),('LP','配体–效应')]:d[label+' rho / P / FDR']=f'{r[e+"_rho"]:+.4f} / {r[e+"_P"]:.6f} / {r[e+"_FDR"]:.6f}'
 d['符号相容']=r.sign_compatible;overview.append(d)
sec('4. 12条同供体链：统一 AUC17 最终校正',p('下表全部 n=17、seed17 DropOverlap AUC、G+S+双方深度，Effector相同定义。配体和表型仍为RNA派生量；“AUC17”指该链采用AUC TF及其固定17人名单，并不表示每条边端点都是AUC。三个相关系数的符号相容只是代数一致性，不是间接效应或中介检验。')+
 table(pd.DataFrame(overview))+p('12链均符号相容，但在该统一主校正规格中，没有一链三段均P&lt;0.05，也没有FDR&lt;0.05。MacAM–Core TNF–效应 P=0.000311、FDR=0.055351，应保留其较强名义关联，同时明确未越过既定FDR门槛。')+
 detail('全部统一最终校正边：RNA/AUC、效应/激活/近端，含区间与LOO',stats(core)))

q=m[m.chain.eq('MacAM/CD4-IL7R/IL18-IKZF2')&auc&eff&m.covariates.eq(FINAL)&m.policy.isin(['aligned_steps','cells20'])]
sec('5. IL18–IKZF2 的方向与实际断点',p('<b>精确对应模型：</b>MacAM→CD4-IL7R，同17位供体，SCENIC DropOverlap seed17，G+S+receiver_library+sender_library。IL18–IKZF2 rho=+0.634403，P=0.026707，FDR=0.247992；IKZF2–Effector rho=−0.737653，P=0.006176，FDR=0.129708；IL18–Effector rho=−0.442630，P=0.149595，FDR=0.421288。')+
 p('因此“上游正、下游负”并不自相矛盾，它对应负向配体–效应方向。这里观察到的第三段点估计确实为负，但其95% CI约[−0.811,+0.176]跨零，不能称为第三段得到统计支持，更不能据此证明低功效解释或链条真实。全链关联也不要求数学上三个边同时显著才可能有复杂机制；本次仅评价既定简单观察链。')+
 stats(q)+p('17人模型的 IL18–IKZF2 LOO 范围约0.532–0.849；IKZF2–效应约−0.871至−0.624；配体–效应约−0.586至−0.110。方向在这些遗漏中保留，但LOO并不解决来源混杂或方向性。≥20细胞的11人最终校正df=4，IL18–IKZF2 rho=0.875490、P=0.022289，而后两段P=0.100638/0.380842，精度更差；不能只突出较大rho。')+
 p('历史 n11 SourceDep（G+S+接收深度）IL18–IKZF2 rho=0.889988、P=0.007260、旧FDR=0.651287；BothEpi进一步加发送深度及sender_epi后rho=0.874173、P=0.052556、旧FDR=0.852854。这些是不同样本/校正的既有结果，不可混记为本次17人主结果。')+
 detail('MacAM 两接收端 IL18 的全部主阶梯与RNA匹配17人结果',stats(m[am&m.axis.eq('IL18')&eff&m.policy.isin(['aligned_steps','RNA_matched_AUC17'])])))

q=m[m.chain.eq('MacAM/Core/TNF-STAT1')&eff&m.policy.isin(['Treg_composition','TNK_competition','IFN_environment','cells20'])]
sec('6. TNF–STAT1：校正依赖与条件性正结果',p('MacAM–Core AUC17 的 TNF–STAT1 主阶梯 rho 约0.036→−0.001→0.073→0.478→0.638，相应最终P=0.025691；关联明显随深度校正改变，不能视为未校正时已稳健存在。最终 STAT1–效应rho=0.544470、P=0.067204，TNF–效应rho=0.861986、P=0.000311。')+
 p('<b>预设 Treg_fraction 敏感性是值得保留的正证据：</b>Core AUC17 的三段rho=0.790406/0.815978/0.872122，P=0.003784/0.002192/0.000464，FDR=0.110743/0.088127/0.070255；RNA19也三段名义显著，P=0.003841/0.005389/0.010441。说明给定组成校正后的条件关联相互一致；不能选择这组规格作为唯一主结论，也不能由加入协变量增强推断其必为正确因果调整。')+
 p('加入TNK_IFNG后，AUC17 TNF–STAT1 P=0.149272、STAT1–效应P=0.520043，而TNF–效应仍P=0.003139。加入冻结IFN_state_NoRegUnion后，TNF–STAT1 P=0.033780、TNF–效应P=0.000777。两种代理所含细胞/信号不同；结果提示环境或组成依赖，但不能确认或排除共享IFN机制，也不能区分混杂、中介和共同反应。IFN代理仅IFIT1/2/3，不能覆盖全部IFN生物学。')+
 detail('TNF–Core敏感性全部端点统计（RNA/AUC）',stats(q))+
 detail('MacAM TNF 两接收端主阶梯与RNA匹配17人',stats(m[am&m.axis.eq('TNF')&eff&m.policy.isin(['aligned_steps','RNA_matched_AUC17'])])))

det=cov[cov.chain.str.startswith('MacAM/')]
sec('7. 受体、配体与近端RNA：能测到什么',p('受体是拓扑/覆盖证据，不要求其RNA随配体RNA正相关。TNF分别看TNFRSF1A、TNFRSF1B；IL18需IL18R1+IL18RAP，同一细胞共检出使用原计数摘要。未按受体检出与否删供体。')+
 p('MacAM–IL7R AUC17：IL18配体RNA在16/17供体检出；IL18R1+IL18RAP同细胞共检出仅5/17供体、11/1371接收细胞，供体中位检出比例为0。Core为9/17供体、32/5684细胞。<b>5/17是稀疏RNA观测与dropout风险，不能作为其余供体不响应或链条不可能的硬排除。</b>细胞数/深度和RNA零值不能代替蛋白检测。')+
 p('MacAM TNF在AUC17中13/17供体检出；IL7R TNFR1为14/17、TNFR2为17/17。RNA零值同时存在于TF：IL7R的STAT1为3/17、IKZF2为4/17；Core为0/17、1/17。AUC并非直接TF蛋白测量。成熟IL18、IL18BP、膜型/可溶TNF、受体蛋白及磷酸化未测，现有RNA不足以判定活性配体供给和真正受体功能。')+
 table(det,['chain','roster','n','ligand_positive_subjects','receptor','receptor_codetected_subjects','receptor_positive_cells','total_receiver_cells','TF_RNA_zero_subjects'])+
 p('统一17人校正下，MacAM–IL7R TNF–NFκB RNA代理rho=−0.348、P=0.268；IL18–NFκB rho=−0.151、P=0.639，IL18–AP1 rho=−0.518、P=0.084。Core近端代理也无名义支持的正相关。RNA19的IL7R IL18–AP1反而为负向名义关联（P≈0.0428）。这些结果不支持“配体RNA越高→近端RNA代理越高”的简单正链，但不能排除瞬时蛋白激活或时间错位。')+
 detail('全部近端代理规格',stats(m[m.edge.eq('ligand_proximal')])) )

modules=rd('reused_module_coverage.tsv'); modules=modules[modules.module.isin(['NFkB_proxy','AP1_components','IFN_state','TregNoTF','Activation'])&modules.variant.eq('NoRegUnion')]
sec('8. 面板覆盖与共享靶基因审计',p('所有表型沿用冻结的去重叠定义：先去除两TF在三个seed的regulon并集及TF自包含，不用后验挑基因增强相关。Effector剩PRF1、IFNG、NKG7，去除GZMA/GZMB；故这是一个窄的效应RNA指标，不能当作全部效应功能。Activation剩CD69、TNFRSF4、IRF4、NR4A1、NR4A2、EGR1，去除ICOS/BATF，IL2RA已预先归入Treg。')+
 p('Treg残余指标只有IKZF4、LRRC32；它与细胞Treg比例不是同一变量。原七基因TregNoTF、残余面板及比例定义均保留在附件，不能笼统称为“已控制全部Treg功能”。近端NFκB/AP1和IFN代理同样规模很小。')+
 table(modules,['module','original_expected','after_overlap','used','used_genes','removed_regulon_genes'])+
 detail('冻结表型逐基因保留/排除',table(rd('phenotype_panel_overlap.tsv')))+detail('regulon重叠原始审计',table(rd('reused_regulon_overlap_audit.tsv'))))

seed=m[am&eff&m.policy.eq('seed_sensitivity')&m.covariates.eq(FINAL)]
sec('9. RNA、AUC与seed一致性',p('固定17人后，RNA与AUC仍不等价：RNA是TF表达，AUC是推断靶集活性；差异可能来自测量、调控层级及regulon不确定性，不能用一个“显著”替换另一个。报告同时保留RNA19、RNA匹配17人与AUC17；全部三段在指定名单和协变量内比较。')+
 p('MacAM–IL7R IL18–IKZF2：seed29/43最终rho=0.415161/0.443591，P=0.179577/0.148609；下游IKZF2–效应仍负向，P=0.003947/0.007114。因此方向比上游显著性更一致。Core IL18–IKZF2仅seed29达到名义显著（P≈0.00918），不能择优汇报。TNF–Core配体–STAT1在seed29 P≈0.0277、seed43 P≈0.0825，下游两seed均未名义显著。')+
 stats(seed)+detail('固定17人RNA与AUC的主要边对照',stats(m[am&eff&m.covariates.eq(FINAL)&((m.policy.eq('RNA_matched_AUC17'))|(auc&m.policy.eq('aligned_steps')))])))

within=rd('within_source_descriptions.tsv')
sec('10. 来源内描述与组成偏倚',p('来源内结果不能称独立验证。AUC G227仅9人、病例2人，按冻结规则仅报告rho/LOO，不给P/FDR。MacAM–IL7R IL18–IKZF2在G227基础/组别校正rho约0.617/0.649，下游IKZF2–效应约−0.417/−0.315；Core IL18同方向。TNF–Core配体–STAT1约0.300/0.292，组别校正LOO可跨零。正向模式可在有限来源内出现，但仍无足够平衡样本支持来源稳健性。')+
 p('G135 AUC只有1人，RNA只有2人，无法作可靠来源内关联；G288为7名病例，仅病例内描述。RNA G227有10人，可按规则给探索性P，但依然只有2病例。所有主模型可拟合；不要把部分来源内不可估计写成全部主模型不可估计。')+
 table(within[within.chain.str.startswith('MacAM/')&within.modality.eq('AUC17')&within.source.eq('G227')&within.edge.isin(['ligand_TF','TF_phenotype','ligand_phenotype'])],['chain','edge','covariates','n','n_case','n_control','rho','P','LOO_min','LOO_max','status'])+
 p('发送端沿用工作注释：MacAM=clusters0/2/11/13（9771细胞），MonoMac=3/7（5135），MacWide=两者合并（14906）。此前审计cluster2约98.8%来自CTL068、cluster11约97.2%来自CTL070，存在供体主导风险；本次没有将已排除的单核/DC/混合簇改标回来。MacAM与MonoMac/MacWide差异用于检查来源和构成，不能挑最高相关者重新定义主发送端。'))

diag=m[am&auc&eff&m.policy.isin(['aligned_steps','cells20','Treg_composition'])&m.covariates.str.startswith(FINAL,na=False)].merge(inf[['spec_id','individual_LOO_available','largest_change_donor','largest_abs_change','sign_flip_donors','note']],on='spec_id',how='left')
sec('11. 精度、共线性与逐供体影响',p('设计诊断分开看：endpoint_VIF衡量端点被协变量解释的程度，max_covariate_VIF衡量协变量间共线性。TNF–Core17人最终模型最大协变量VIF约10.2，加入Treg后约12.3；配体端点VIF约2，二者不能混为同一风险。17人主最终校正df=10，11人同规格df=4，更多协变量会继续减少df。critical_abs_rho是给定df下达到名义P&lt;0.05所需|rho|，不是事后功效估计。')+
 p('MacAM–IL7R17人IL18–IKZF2最大遗漏变化来自CWP127，|Δrho|≈0.215；这说明单供体对幅度有影响，不能因LOO方向全同而称“无单样本影响”。历史已保存的逐人遗漏全部汇入available_per_donor_influence；只有范围、没有逐人标签的历史行明确标记，未编造遗漏供体，也未扩展补算。')+
 table(diag,['spec_id','chain','edge','policy','n','df','condition_number','endpoint_VIF','max_covariate_VIF','nominal_critical_abs_rho','LOO_valid','LOO_min','LOO_max','individual_LOO_available','largest_change_donor','largest_abs_change','sign_flip_donors'])+
 p('缺失的12条历史LOO已完成193个遗漏诊断。写回前后1368行rho、P、FDR逐值完全一致；未重跑主效应。剩余LOO范围缺失数=0，但这不意味着每条历史规格都保留了逐人明细。'))

classes=pd.DataFrame([
 ['未测量的功能环节','全部链：成熟/有效配体、受体蛋白和信号蛋白激活未测','不能由配体或受体RNA正负值判断活性；不是主统计不可拟合'],
 ['测量稀疏与精度有限','IL18受体共检出5/17；11人模型df低、区间宽','不是硬性排除；不能据P值断言真实效应被低功效掩盖'],
 ['方法、来源、组成敏感','TNF–Core随深度/Treg/TNK校正变化；来源病例失衡；seed差异','有条件关联值得保留；需要限定适用模型，不作因果选择'],
 ['可能的候选角色不匹配','IKZF2下游负相关支持抑制性方向，而非预设正向效应增强；STAT1可能反映共享IFN环境','是解释候选，不是已证实错基因；两TF保留'],
 ['既定观察链未完整支持','统一主最终校正中每链至少一段P≥0.05；近端RNA简单正向链未支持','不能称完整统计链已成立；也不能等同生物通路被证伪'],
 ['尚不确定','全部12链FDR未达0.05；无时间/干预方向证据','方向相容、名义显著、条件性一致与因果证明分级表达']
],columns=['断点类别（可重叠）','实际证据','允许的结论'])
sec('12. 断点分类与逐链阅读规则',table(classes)+p('逐链主矩阵附受体检测、近端统计和证据标签；全部1368规格附实际样本ID、端点、校正、效应、区间、P、FDR、旧FDR、复用来源及LOO。对MacAM以外发送端不给等同主链的证据等级；MonoMac/MacWide的名义结果完整公开，不择优换主线。')+detail('12链断点标签原表',table(matrix)))

lit=[
 ['人体疾病组织背景','矽肺患者肺组织观察到巨噬细胞TNF与上皮NFκB活化相关背景','患者肺相关，但上皮结果不能跨接成CD4–STAT1通路；动物系统抑制实验另属不同证据层','https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0005689','原文已查阅'],
 ['人体T细胞干预、亚群限定','TNF在人CD45RA−Treg诱导NFκB程序，并调节抑制功能；抗TNFR2可逆转所述功能变化','支持TNFR2/Treg情境依赖；不是普通肺CD4中TNF直接激活STAT1的证明','https://pubmed.ncbi.nlm.nih.gov/20181891/','原始论文摘要已核对'],
 ['人体CD4干预、条件限定','人外周CD4中IL12+IL18诱导一小部分细胞产生IFNγ；单独IL18不呈相同反应','支持共同刺激与受体情境的重要性；未证明IL18→IKZF2，也未直接验证患者肺','https://pubmed.ncbi.nlm.nih.gov/21572994/','原始论文摘要及图注已核对'],
 ['人体T细胞功能、年龄限定','2026年研究以单细胞/染色质与CRISPR证据识别HELIOS(IKZF2)对婴儿T细胞程序及效应功能的限制作用','为IKZF2负向效应方向提供生物学相容性；不能外推成成人患者肺或IL18上游直接证据','https://pubmed.ncbi.nlm.nih.gov/42230957/','原始论文摘要与图注已核对'],
 ['其他细胞机制背景','TNF可经IRF1依赖的自分泌IFN环路产生STAT1依赖转录反应','原研究以巨噬细胞为背景；提供共享IFN/间接机制候选，不证明受体到CD4 STAT1的直接边','https://pubmed.ncbi.nlm.nih.gov/18345002/','前轮原始摘要已核对，本轮页面读取不稳定'],
 ['人体T细胞干预、待全文复核细节','父任务提供的健康Th1/Th17研究提示TNF/TNFR2作用可随培养条件呈效应抑制','作为方向不应预设单一正向的有限线索；剂量、供体数等细节不用于本报告核心论证','https://pubmed.ncbi.nlm.nih.gov/36012570/','本轮全文受访问限制；依父任务已查摘要/原文线索，降级标注']
]
litdf=pd.DataFrame(lit,columns=['等级/情境','可支持内容','不可跨越的边界','原始来源','核查范围']);litdf.to_csv(R/'literature_evidence_levels.tsv',sep='\t',index=False)
litbody='<div class="scroll"><table><thead><tr>'+''.join('<th>'+c+'</th>' for c in litdf.columns)+'</tr></thead><tbody>'
for row in lit:
 litbody+='<tr>'+''.join('<td>'+(f'<a href="{v}">原始论文</a>' if i==3 else html.escape(v))+'</td>' for i,v in enumerate(row))+'</tr>'
litbody+='</tbody></table></div>'
sec('13. 文献分级：背景支持不能拼接成患者机制',p('这是一轮针对既定两轴的有限原始文献核对，不是系统综述。将人体患者组织、人体细胞干预和跨细胞背景分开；摘要可核对的结论与未能独立核对的全文细节分开。没有在本轮有限核对中找到可直接建立“成人患者肺巨噬细胞IL18→CD4 IKZF2”的匹配干预证据；这不是宣称全世界不存在此证据。')+litbody)

triggers=pd.DataFrame([
 ['先保留现有候选','当前同供体统计有方向相容和条件性正关联，但FDR/来源/测量限制尚在','不因为P=0.150或受体5/17就换基因；先将结论写成有限观察证据'],
 ['若目标是跨来源可复现关联','G135仅对照、G288仅病例、G227只有2病例','在需要进一步证据时优先病例/对照来源平衡、同供体巨噬与CD4匹配的患者样本；不要求现在下载新数据'],
 ['若决定取决于受体/配体是否真正可用','RNA共检出稀疏且活性配体/受体蛋白未测','选能回答成熟IL18/受体蛋白或TNF形式的测量；先明确科学问题，不将任何实验当作完成本轮的前置条件'],
 ['若要求直接方向或中介因果','横断面共同变化不能判定TNF–STAT1共享IFN与直接作用','届时才考虑时间或干预证据；本轮不做高维中介或因果声称'],
 ['何时才考虑更换候选角色','在合适细胞情境与足够测量覆盖下，独立数据方向持续不相容，或有窄区间排除预定有意义效应，或直接功能证据显示角色不匹配','这些条件目前尚未满足；单一非显著P、单seed或小样本筛选不足以触发更换']
],columns=['决策','具体触发条件','对应行动边界'])
sec('14. 新数据或更换候选的具体触发条件',table(triggers)+p('<b>本轮到此停止分析扩展。</b>保留TNF/IL18及STAT1/IKZF2，不增加配体、TF、阈值、聚类或筛选网格；线粒体与焦亡不作为当前主线或交付前置条件。没有下载新的单细胞数据集。'))

sec('15. 交付、核验与复现',p('交付包含：中文HTML报告；完整1368行逐链证据TSV；结果和代码ZIP。ZIP含12链汇总、180个同供体符号组合、全部主规格输入、现存逐供体影响、来源内384行描述、受体覆盖、冻结计划、模型复用/完成记录、独立核验和12条LOO补充证据。')+
 p('850个新拟合已独立秩变换+OLS核验，最大rho差4.44×10⁻¹⁶、最大P差2.95×10⁻¹⁴；482条历史统计复制核验通过；BH和固定名单核验通过；9个已登记模型输入源文件哈希未变（不扩大为“所有项目文件均核验”）。报告生成再次确认rho/P/FDR写回前后逐值一致。')+
 p('代码依赖原项目的冻结前序输入和既有结果路径；ZIP提供本轮逐规格输入供审阅，但不复制全项目或原始单细胞矩阵。原分析脚本有完成保护，不要在完成目录删除标记重跑；如需将来独立复现，应复制到新目录并显式配置输入。frozen_plan.pkl可由冻结计划与原项目输入重建，作为中间缓存不纳入交付。')+
 detail('已完成核验记录', '<pre>'+html.escape(json.dumps({'verification':verify,'LOO_completion':loo},ensure_ascii=False,indent=2))+'</pre>'))

sec('16. 全部1368规格可检索明细',p('下表可按链名、spec_id或策略关键词筛选。完整字段（包括精确sample_ids、legacy_FDR、复用来源、VIF与零值数）见随附TSV；表格数值仅显示精度舍入。')+
 '<input id="search" type="search" placeholder="筛选，例如 MacAM/Core/TNF、BP00001、cells20" aria-label="筛选全部模型"><span id="count"></span><div id="allmodels">'+stats(m)+'</div>')

css='''body{margin:0;background:#f3f6fa;color:#18283b;font:16px/1.75 "Microsoft YaHei",Arial,sans-serif}header{background:#142b44;color:white;padding:48px max(28px,calc((100vw - 1250px)/2));}h1{font-size:32px;margin:8px 0;line-height:1.4}header p{color:#d1e0ed;max-width:1050px}main{max-width:1320px;margin:auto;padding:24px}nav,section{background:white;padding:25px 30px;border-radius:10px;margin-bottom:22px;border:1px solid #dde5ee}h2{font-size:23px;color:#174264;line-height:1.5}h3{font-size:19px}.scroll{overflow:auto;max-height:650px;margin:20px 0;border:1px solid #dae2eb}table{border-collapse:collapse;width:100%;font-size:13px;line-height:1.5}td,th{padding:10px 12px;border-bottom:1px solid #e3e9ef;text-align:left;vertical-align:top;min-width:66px}th{position:sticky;top:0;background:#e8f0f7;color:#173c5c;z-index:1}tr:nth-child(even){background:#f8fafc}td:first-child{min-width:180px;overflow-wrap:anywhere}.cards{display:grid;grid-template-columns:repeat(3,1fr);gap:16px}.cards article{background:#f0f6fa;border-top:4px solid #36738c;padding:18px;border-radius:6px}a{color:#116a94}nav a{display:inline-block;margin:4px 20px 4px 0;font-size:14px}details{margin:18px 0;border:1px solid #d6e2ec;border-radius:6px;padding:12px}summary{cursor:pointer;color:#1c557a;font-weight:bold}pre{white-space:pre-wrap;font-size:13px}input{padding:12px;width:min(600px,80%);font-size:15px;border:1px solid #9eb6ca;border-radius:5px}#count{padding:10px}footer{padding:20px;color:#536a7e;font-size:13px}@media(max-width:800px){main{padding:12px}.cards{grid-template-columns:1fr}section,nav{padding:18px}h1{font-size:26px}}@media print{body{background:white}section{break-inside:auto}.scroll{max-height:none;overflow:visible}details{display:block}nav,input,#allmodels{display:none}table{font-size:9px}header{padding:20px}main{padding:0}}'''
js="""const inp=document.querySelector('#search'),rows=[...document.querySelectorAll('#allmodels tbody tr')];function filter(){let n=0,q=inp.value.toLowerCase();rows.forEach(r=>{let ok=r.textContent.toLowerCase().includes(q);r.hidden=!ok;if(ok)n++});document.querySelector('#count').textContent=n+' / '+rows.length+' 行';}inp.addEventListener('input',filter);filter();"""
doc='<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>TNF/IL18 同供体链断点调查</title><style>'+css+'</style></head><body><header><div>患者肺 · 巨噬细胞 → CD4 T 细胞 · 有限范围详细调查</div><h1>TNF / IL18 同供体链断点调查</h1><p>保留 STAT1 与 IKZF2｜12链 · 1368规格｜主统计已核验，12条历史LOO已补齐<br>版本：2026-10-04 · 固定发送端、接收端、供体、seed与校正条件逐段解释</p></header><main><nav>'+''.join(f'<a href="#{i}">{t}</a>' for i,t in toc)+'</nav>'+''.join(parts)+'</main><footer>结果为观察关联；所有结论以指定模型和患者样本覆盖为界。报告/表格/代码共同构成交付。</footer><script>'+js+'</script></body></html>'
(O/'lung_TNF_IL18_breakpoint_report_zh.html').write_text(doc,encoding='utf-8')
shutil.copyfile(R/'aligned_chain_models.tsv',O/'lung_TNF_IL18_chain_evidence.tsv')
shutil.copyfile(Path(__file__),O/'code/build_breakpoint_delivery.py') if Path(__file__).resolve()!=(O/'code/build_breakpoint_delivery.py').resolve() else None
readme='''TNF/IL18 同供体链断点调查 — 完成交付

先打开 lung_TNF_IL18_breakpoint_report_zh.html。该报告可离线浏览，完整1368规格可筛选。
lung_TNF_IL18_chain_evidence.tsv 是全部1368行证据，不是仅显著结果。
results/chain_breakpoint_matrix.tsv 是12链统一AUC17最终校正汇总。
results/same_donor_chain_sign_audit.tsv 是180组可逐段比较的同供体符号审计。
results/aligned_model_inputs.tsv.gz 是每个规格的实际端点、供体与协变量输入。
results/available_per_donor_influence.tsv.gz 是已有逐供体LOO明细；部分历史行仅有范围，见influence_summary.tsv。
results/MacAM_IL7R_19_17_11_trace.tsv 是19→17→11逐人损失表。
results/literature_evidence_levels.tsv 是分级文献证据。

统计列：rho=偏秩相关；P=探索性t近似；CI95=Fisher z近似；FDR=本轮BH家族；legacy_FDR=历史原家族。
G=病例组别；S=来源；receiver_library/sender_library=两端深度；Treg_fraction与TNK_IFNG/IFN_state各为独立预设敏感性。
sample_ids和covariates共同决定模型条件；不能跨名单拼接链。

原始执行顺序（已完成，勿在原目录重复执行）：
  investigate_breakpoints.py --plan
  investigate_breakpoints.py --run
  audit_breakpoint_results.py
  complete_breakpoint_loo.py
  summarize_breakpoints.py
报告可仅用已有结果重新渲染：python code/build_breakpoint_delivery.py <交付目录>
本报告脚本不拟合模型，仅读取、核对、生成展示与清单。

原项目根：/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA
独立目录：11_reports/breakpoint_TNF_IL18_20261004_v1
原科学环境：/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
统计代码依赖原项目front_half、missing_links、aligned_evidence冻结结果；包不包含全项目原始矩阵。
冻结计划在code/breakpoint_plan.txt；输入文件路径/哈希在results/input_manifest.tsv。
frozen_plan.pkl为可重建缓存，不打包。已完成标记是保护，勿删除以便在原目录重跑。

核验：850新拟合/482历史复制、名单与BH均通过；12历史LOO补193诊断后主rho/P/FDR保持完全一致。
本轮未新增候选/阈值/数据集，线粒体和焦亡未前置，不作因果中介推断。完成后停止分析扩展。
'''
(O/'README.txt').write_text(readme,encoding='utf-8')
receipt={'status':'COMPLETE_VERIFIED','time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'specifications':1368,'chains':12,'new_unique_fits':850,'historical_reused_specifications':482,'new_signature_reused_specifications':36,'LOO_missing_completed':12,'omission_diagnostics':193,'main_rho_P_FDR_unchanged':True,'remaining_analysis':0,'report_and_package_generated':True}
(P/'complete.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8');(O/'STATE').write_text('COMPLETE_VERIFIED\n',encoding='utf-8')
exclude={'frozen_plan.pkl','manifest.tsv','library_identity_receipt.json','report_preview.png','delivery_verification.json','delivery.log','analysis.lock'}
files=[p for p in O.rglob('*') if p.is_file() and p.name not in exclude and p.suffix not in ['.zip','.pyc'] and '__pycache__' not in p.parts]
manifest=pd.DataFrame([{'path':p.relative_to(O).as_posix(),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(files)])
manifest.to_csv(O/'manifest.tsv',sep='\t',index=False)
zpath=O/'lung_TNF_IL18_breakpoint_results_code.zip'
with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for f in files+[O/'manifest.tsv']:z.write(f,f.relative_to(O).as_posix())
with zipfile.ZipFile(zpath) as z:
 assert z.testzip() is None
 for _,r in manifest.iterrows():assert hashlib.sha256(z.read(r.path)).hexdigest()==r.sha256
print(json.dumps({'status':'PASS','report_bytes':len(doc.encode()),'models':len(m),'archive_files':len(files)+1,'archive_bytes':zpath.stat().st_size,'archive_sha256':hashlib.sha256(zpath.read_bytes()).hexdigest()},ensure_ascii=False,indent=2))
