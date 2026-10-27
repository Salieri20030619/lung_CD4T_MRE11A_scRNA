from pathlib import Path
import os,json,html,hashlib,zipfile,datetime
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np,pandas as pd
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/front_half_20261004_v1');R=O/'results';P=O.parent.parent
def rd(n):return pd.read_csv(R/n,sep='\t',low_memory=False)
def wr(d,n):pd.DataFrame(d).to_csv(R/n,sep='\t',index=False,na_rep='NA')
assert json.loads((O/'provenance/verification.json').read_text())['status']=='PASS'
m=rd('new_link_models.tsv');final=rd('fixed_final_models.tsv');rec=rd('receiver_audit_inputs.tsv');sen=rd('sender_state_and_ligands.tsv');comp=rd('sender_composition.tsv');gm=rd('group_context_models.tsv')
review=P/'11_reports/technical_delivery_20261002_ikzf2_stat1_v1/revision_v2/final_release/prior_review_repaired/scenic3_review/results/reviewed_sender_expression.tsv'
oldsen=pd.read_csv(review,sep='\t');oldsen=oldsen[oldsen.sender_filter.eq('All')]
gates=[]
for (scope,sender),z in pd.read_csv(R/'paired_inputs.tsv',sep='\t').groupby(['scope','sender']):
 for ligand in ['IL15','TNF','IL18','CXCL16']:
  rr=['R_TNFR1','R_TNFR2'] if ligand=='TNF' else ['R_'+ligand]
  q=oldsen[oldsen.sender.eq(sender)&oldsen.ligand.eq(ligand)]
  for receptor in rr:
   gates.append(dict(scope=scope,sender=sender,ligand=ligand,receptor_endpoint=receptor,paired_subjects=len(z),sender_positive_subjects=int(z['L_'+ligand].gt(0).sum()),receptor_positive_subjects=int(z[receptor].gt(0).sum()),receptor_positive_pct_median=z[receptor].median(),IL15_IL15RA_joint_subjects=int(q.IL15_IL15RA_pct.gt(0).sum()) if ligand=='IL15' else np.nan,minimum_sender_cells=z.sender_cells.min(),minimum_receiver_cells=z.receiver_cells.min(),interpretation='transcript coverage, not functional ligand/receptor or causal communication'))
wr(gates,'finite_evidence_coverage_matrix.tsv')
wr(rec[['scope','sample','group','series','receiver_cells','Treg_fraction','TregNoTF_full','TregNoTF_NoRegUnion','State_Treg','IKZF2_RNA','STAT1_RNA','IKZF2_AUC_17','STAT1_AUC_17']],'treg_original_fraction_residual_side_by_side.tsv')
old=rd('reused_TF_state_and_TNK_competition.tsv')
q=old[(old.species=='Human')&(old.kind=='receiver')&old.y.isin(['TregNoTF_full','TregNoTF_NoRegUnion'])&old.x.isin(['IKZF2_RNA','IKZF2_AUC_17'])]
wr(q,'treg_original_and_residual_existing_models.tsv')
newstats=m.groupby('FDR_family').agg(specification_rows=('model_id','size'),unique_tests=('test_hash','nunique'),nominal_specifications=('nominal_P_lt_05','sum'),FDR05_specifications=('FDR_lt_05','sum')).reset_index();wr(newstats,'FDR_family_summary.tsv')
nn=rd('reused_step3_nichenet_ranks.tsv');nn=nn[(nn['filter']=='All')&(nn.test=='FDR05')&(nn.direction=='Up')]
matrix=[]
for (scope,sender,lig),z in final[final.kind.eq('ligand')].assign(ligand=lambda x:x.x.str.replace('L_','',regex=False)).groupby(['scope','sender','ligand']):
 cover=pd.DataFrame(gates);c=cover[(cover.scope==scope)&(cover.sender==sender)&(cover.ligand==lig)]
 oldtf=rd('reused_reviewed_sender_TF.tsv');t=oldtf[(oldtf['filter']=='All')&(oldtf.scope==scope)&(oldtf.sender==sender)&(oldtf.ligand==lig)&(oldtf.seed==17)&(oldtf.sample_set=='Matched')&(oldtf.model=='SourceDep')]
 mouse=nn[nn.test_ligand.str.upper().eq(lig)]
 matrix.append(dict(scope=scope,sender=sender,ligand=lig,n=19,receiver_coverage=';'.join(c.receptor_endpoint+':'+c.receptor_positive_subjects.astype(str)+'/19'),new_final_nominal_states=';'.join(z.loc[z.P.lt(.05),'y']),new_final_FDR05_states=';'.join(z.loc[z.FDR.lt(.05),'y']),legacy_sourceTF_rows=len(t),legacy_sourceTF_min_P=t.p.min() if len(t) else np.nan,legacy_sourceTF_min_FDR=t.FDR.min() if len(t) else np.nan,mouse_Core_FDRset_Up_rank=float(mouse['rank'].iloc[0]) if len(mouse) else np.nan,claim='exploratory source/state association only; no demonstrated source-to-TF-to-phenotype chain'))
wr(matrix,'front_half_evidence_matrix.tsv')
def fmt(v):
 if pd.isna(v):return '—'
 if isinstance(v,(float,np.floating)):return f'{v:.4g}'
 return str(v)
def table(d,cls=''):
 return d.to_html(index=False,escape=True,formatters={c:fmt for c in d.columns},classes=cls,border=0,na_rep='—')
parts=[]
def p(t):parts.append('<p>'+t+'</p>')
def h(t):parts.append('<h2>'+t+'</h2>')
parts.append('<div class="eyebrow">LUNG PROJECT · FINITE EVIDENCE AUDIT · 2026-10-04</div><h1>巨噬来源与 T 细胞表型：<br>IKZF2 / STAT1 前半段补强</h1>')
p('<strong>当前可支持的是患者背景下的细胞组成、来源候选与 T 细胞状态关联。</strong>IKZF2 / STAT1 主线保留；尚不能把“巨噬配体→两 TF→T 表型”写成已证实机制。线粒体损伤、OXPHOS 与焦亡未作为本轮成立门槛。')
parts.append('<div class="metrics"><div><b>19</b>人供体 · RNA</div><div><b>17</b>冻结 AUC 供体</div><div><b>185 / 1,675</b>关系 / 规格</div><div><b>70</b>直接复用旧统计</div></div>')
p('新增计算 1,605 条关系规格，另有 125 条组别背景规格；1,675 条关系规格全部可估计。规格数包含校正和敏感性，绝不是新增生物重复。全部逐项输出 n、样本名单、效应、名义 P、FDR 与 LOO；本报告按用户要求将 P&lt;0.05 标为“名义显著”，同时显示 FDR。')
h('1. 既有结果与真实缺口')
p('早期 SCENIC3 marker 复核和 sender–TF 已覆盖 MacAM、MonoMac、MacWide；不能说过去完全遗漏 MonoMac。近期 missing_links 的 8,264 条模型主要聚焦 AM 与 TNK_other。此次补的是巨噬比例、亚群内状态与 T 活化／初始样分化／抑制标记／效应状态的对应，而不是再次计算旧 8,264 条模型。')
p('保留原注释：MacAM=cluster 0/2/11/13；MonoMac=3/7；MacWide=两者并集。它们是 marker 支持的工作注释，MacWide 并非独立来源或独立重复。原单核、DC 与混合群未改名为巨噬；未把连续 SPP1Mac/FOLR2Mac 分数命名成新亚型。')
parts.append(table(rd('cluster_donor_dominance.tsv')[lambda x:x.cluster.isin([0,2,3,7,11,13])][['cluster','scope','cells','subjects','dominant_sample','dominant_fraction']]))
p('映射成功的 MacAM=9,771 细胞，MonoMac=5,135，MacWide=14,906；另有 130 个 Unknown MacAM 只保留注释，不进入供体模型。cluster 2/11 的 CTL068/CTL070 主导风险保留，并有 NoDom 敏感性。比例分母为已选四个原父群的髓系候选范围，不能外推为整肺巨噬绝对丰度。')
h('2. 患者背景与样本设计')
parts.append(table(rd('source_group_design.tsv')))
p('沿用原 CWP / Control 标签（9 / 10），本轮不重命名病例或补充临床诊断。G135 仅对照，G288 仅病例，G227 为 2 病例 / 8 对照；来源校正后的疾病系数主要由有限的来源内重叠支撑，不能称独立队列复制。RNA 包含全部 19 位；CTL070 / CTL103 不在冻结 AUC 名单，未补零。')
gc=gm[gm.y.isin(['MacWide_fraction_candidate','MacAM_fraction_candidate','Fraction_Treg','Fraction_CD4-CCR7'])][['scope','y','covariates','n','effect','P','FDR','LOO_min','LOO_max']]
parts.append(table(gc))
p('未校正时，病例候选髓系范围中的 MacWide / MacAM 比例较低，Core Treg 比例较高，CCR7 亚群比例较低，均仅名义显著。加入来源后这些组别差异不再名义显著；本轮 125 条背景规格没有 FDR&lt;0.05。不能据此宣传患者特异巨噬新亚型或稳定的亚群内疾病效应。')
h('3. 最终充分校正关系：名义 P 与 FDR 并列')
p('配体模型逐步为：基础→组别→组别+来源→加 receiver 深度→加 sender 深度；所有步骤固定同一批供体。Core TF 模型最终再加 Treg 比例；MonoMac 状态最终再加 MacAM / MacWide 比例。偏秩相关效应为 rho；P 为探索性 t 近似。最终 185 个关系没有 FDR&lt;0.05。下表列出全部 17 个最终名义显著关系，不只挑正向结果。')
show=final[final.P.lt(.05)][['kind','scope','sender','x','y','n','rho','P','FDR','LOO_min','LOO_max']].sort_values(['kind','P']);parts.append(table(show))
p('早期步骤并非全阴性：MacAM 候选范围比例与两基因 Treg 残余分数的基础关联通过本轮 FDR，加入组别后仍有一条 FDR 阳性，进一步来源校正后消失。不可把基础阳性和最终校正阴性混称为全程稳健。')
h('4. 哪些方向较稳，哪些依赖条件')
claims=[('IKZF2_AUC_17','State_Effector','CD4-IL7R','receiver','IKZF2 活性与效应状态负向'),('STAT1_AUC_17','State_Effector','Core','receiver','STAT1 活性与效应状态正向'),('L_TNF','State_Effector','Core','MacAM','MacAM TNF 与 Core 效应状态'),('L_TNF','State_Activation','CD4-IL7R','MacWide','MacWide TNF 与 IL7R 活化状态'),('SPP1Mac','State_Activation','Core','MonoMac','MonoMac 内 SPP1 程序与活化状态'),('L_CXCL16','State_Inhibitory','Core','MonoMac','MonoMac CXCL16 与抑制标记状态')]
for x,y,scope,sender,title in claims:
 q=m[(m.x==x)&(m.y==y)&(m.scope==scope)&(m.sender==sender)]
 parts.append('<details open><summary>'+title+'</summary>'+table(q[['covariates','sample_policy','n','rho','P','FDR','LOO_min','LOO_max']])+'</details>')
p('<strong>TF 状态方向：</strong>IL7R 中 IKZF2 AUC→效应分数，n=17，最终 rho=-0.8276、P=0.000480、FDR=0.08846；≥20细胞与 NoDom 的方向和名义显著性保留。Core STAT1 AUC→效应分数，n=17，最终 rho=+0.8344、P=0.000738、FDR=0.08846；方向保留，但明显受深度和 Treg 比例校正影响，来源校正阶段 P=0.1178。两者 LOO 均保持方向；这是方向稳定的探索性状态关系，未达本轮 FDR。')
p('<strong>巨噬来源：</strong>MacAM TNF→Core 效应状态最终 P=0.01868，但加 TNK_IFNG 后 P=0.1217；MacWide TNF→IL7R 活化状态最终 P=0.01288，但≥20细胞 P=0.1849、加 TNK_IFNG 后 P=0.0926。MonoMac SPP1 状态→活化仅在加入比例校正后达到 P=0.04432，≥20细胞 P=0.2402。MonoMac CXCL16→抑制标记分数为负，最终 P=0.001437，但≥20细胞 P=0.5494 且 LOO 可变号。均未获 FDR 支持，不能选取其中最有利规格当定论。')
p('TNK 只用于既有竞争来源解释和预定 IFNG 协变量，不替代巨噬主线。TNK 校正可削弱相关，但不能证明 TNK 因果驱动；也可能反映共享炎症状态及测量相关。NoDom 保持方向同样不能消除来源混杂。')
h('5. IKZF2 / STAT1：RNA、活性、Treg 组成分开')
parts.append(table(rd('reused_human_original_TF_RNA.tsv')[['gene','scope','logFC_original','PValue_original','FDR_original']]))
p('保留历史计数模型：人 Core IKZF2 RNA log2FC=+1.514、P=0.01690、FDR=0.8388，为名义阳性；STAT1 RNA 在对应模型并非疾病差异阳性。新均值／状态关联不取代原 edgeR，也不能把 AUC 与 RNA 合并为同一种证据。完整 RNA 固定到同17位 AUC 供体的敏感性已另列。')
parts.append(table(rd('panel_coverage.tsv')))
p('<strong>Treg 去重后仅 IKZF4、LRRC32 两基因，解释受限。</strong>原 TregNoTF 七基因评分、原注释 Treg 比例、两基因残余分数逐供体并列在 treg_original_fraction_residual_side_by_side.tsv；原面板已有统计保存在 treg_original_and_residual_existing_models.tsv。两基因分数不能代表完整 Treg 身份或抑制功能，亦不能因为名义阳性再增删基因。')
p('五个新使用的状态面板彼此不共享基因，均去除两 TF 三种子原 regulon 靶基因并集。Activation 与 Treg 去重评分经逐供体核实与旧结果完全相同，最大差=0；70 条同供体同校正规格直接复用统计和原 LOO 范围，原 FDR 保留在 legacy_FDR，新家族 FDR 另列。其他状态分数也是预定小 RNA 指标，不是功能实验。')
h('6. 配体／受体与鼠 NicheNet 如何补充')
parts.append(table(pd.DataFrame(gates)[lambda x:x.sender.eq('MonoMac')][['scope','ligand','receptor_endpoint','paired_subjects','sender_positive_subjects','receptor_positive_subjects','IL15_IL15RA_joint_subjects']]))
p('全部 24 个来源×receiver×配体组合见 front_half_evidence_matrix.tsv；受体逐来源覆盖见 finite_evidence_coverage_matrix.tsv。IL15 模型使用单基因 sender RNA，须与 IL15/IL15RA 共检出和 receiver IL2RB/IL2RG 分开解释；TNF 两受体分别审计；IL18 receiver 共检出在人 IL7R 仍偏稀疏；CXCL16/CXCR6 转录不能证明募集速率。新配体关联不能绕过这些覆盖限制。')
parts.append(table(nn[['test_ligand','receiver','test','direction','rank','aupr_corrected','geneset_n','background_n']]))
p('鼠 Core 的既有 FDR 上调靶程序中 IL15 / TNF 的相对 rank 为 3 / 9，IL18 / CXCL16 为 65 / 80，保留为辅助优先级。该86基因疾病程序含 FOXP3、IKZF4、IL2RA、TNFRSF4、PDCD1，与 Treg／活化／抑制标记有描述性交集；它不是这些配体特有的已验证靶集。ActCD4 FDR 集仅1基因，原规则不可打分，不能用其 nominal 集偷偷替换。完整 All/BG1、FDR/Nominal、上下调和预定敏感性排名均复用交付。')
p('没有把鼠 prior 或疾病基因集移植成人类 NicheNet 验证，也没有重新排名以挑选支持。NicheNet 的配体活性用于给候选与靶程序联系排序，不等同实验干预或因果方向；参见 <a href="https://github.com/saeyslab/nichenetr/blob/master/vignettes/ligand_activity_geneset.md">NicheNet 官方方法</a>。')
h('7. 可用于患者主线的表述与停止边界')
parts.append('<blockquote>在现有人肺 CWP / Control 数据中，marker 复核支持肺泡／驻留样与单核–巨噬样工作范围。巨噬比例和亚群内状态与 CD4 表型存在有限、条件敏感的探索关联。IKZF2 与 STAT1 的表达及推断活性应分开解释：IKZF2 活性与部分效应状态负向，STAT1 活性与部分效应状态正向，方向在若干敏感性中保留；尚不足以确定某一巨噬配体经两 TF 导致 T 表型变化。</blockquote>')
p('本次可支持：来源注释与配体转录候选、receiver 转录覆盖、鼠靶程序辅助、两 TF 与人 T 状态的名义关系。仍未覆盖：功能性配体／受体、患者内独立来源效应、直接 TF 靶向与方向性干预。不能把不同链段的相关拼成中介或因果链；不再扩候选、改阈值或增加后验面板。线粒体损伤与焦亡继续保留为探索方向。')
h('8. 核验、完整结果与复现')
verification=json.loads((O/'provenance/verification.json').read_text());parts.append('<pre>'+html.escape(json.dumps(verification,ensure_ascii=False,indent=2))+'</pre>')
p('独立线性回归核验全部 1,675 条关系的偏秩统计，rho 最大误差 4.58e-16，P 最大误差 1.81e-14；BH、固定名单和 LOO 范围一致，19 个输入哈希未变。旧统计复用行只有历史 LOO 摘要，未虚构逐样本旧 LOO；新计算行提供逐一留出结果。')
parts.append(table(newstats))
p('BH 家族跨本次可计算唯一检验：primary（MacAM/MonoMac、比例、状态、TF主规格）、sensitivity（MacWide、≥20细胞、NoDom、TNK、RNA匹配供体）及单独 group-context 家族。相同数值规格按样本/秩/设计去重；家族之间可能含重复数值，不能当独立证据。所有新旧模型字段和解释范围见 README.txt。')
p('服务器目录：<code>'+str(O)+'</code>。代码、冻结方案、所有统计表、逐模型输入／名单／LOO、复用历史表、校验记录均在结果代码 ZIP。未改旧报告、未改网络或凭据，未启动科研任务之外的服务。')
parts.append('<details><summary>展开全部185个最终关系（含阴性）</summary>'+table(final[['model_id','kind','scope','sender','x','y','n','rho','P','FDR','LOO_min','LOO_max']])+'</details>')
css='''body{font:15px/1.65 "Microsoft YaHei","Noto Sans CJK SC",sans-serif;color:#1b2c3e;background:#f4f6f9;margin:0}main{max-width:1250px;margin:auto;background:white;padding:48px}h1{font-size:36px;line-height:1.3;color:#123b52;margin:12px 0 24px}h2{margin-top:42px;border-top:2px solid #d4e1e9;padding-top:22px;color:#174860}.eyebrow{letter-spacing:2px;font-size:12px;color:#627b8b}.metrics{display:flex;gap:12px;flex-wrap:wrap;margin:26px 0}.metrics div{background:#edf5f7;padding:15px 24px;border-radius:8px}.metrics b{display:block;font-size:24px;color:#14566d}table{border-collapse:collapse;display:block;overflow:auto;width:100%;font-size:12px;margin:18px 0}th{background:#e7f0f3;text-align:left;white-space:nowrap}td,th{padding:7px 9px;border-bottom:1px solid #dde5e9;vertical-align:top}tr:nth-child(even){background:#f8fafb}details{padding:8px 0}summary{cursor:pointer;color:#17566a;font-weight:bold}blockquote{border-left:4px solid #2d7182;background:#eef6f7;margin:20px 0;padding:20px}pre{font-size:12px;background:#f4f6f8;white-space:pre-wrap;padding:16px}code{overflow-wrap:anywhere}a{color:#086684}@media print{main{padding:10px}details{display:block}table{font-size:9px}h2{break-after:avoid}}'''
report='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>肺项目前半段有限补强</title><style>'+css+'</style><main>'+''.join(parts)+'</main></html>'
(O/'lung_front_half_report_zh.html').write_text(report,encoding='utf-8')
readme='''肺项目前半段有限补强 / 2026-10-04
先读 lung_front_half_report_zh.html。
results/new_link_models.tsv: 全1675规格，statistic_origin区分新计算与复用。P/FDR是本轮家族；legacy_FDR保留历史值。nominal_P_lt_05 按用户标准标记；不能把它等同FDR。
results/fixed_final_models.tsv: 185关系最终校正规格。
results/fixed_adjustment_transitions.tsv: 同供体基础到最终校正并列。
results/new_link_model_inputs.tsv.gz / new_link_rosters.tsv.gz / new_link_LOO.tsv.gz: 每条规格真实样本值、纳入名单、逐样本LOO。复用行LOO仅历史汇总，见主表；未造旧逐样本LOO。
results/group_context_models.tsv / group_context_inputs.tsv / group_context_LOO.tsv: 125条疾病背景规格，等权样本OLS；不是原edgeR。
results/front_half_evidence_matrix.tsv: 24个来源/receiver/配体证据链单元。
results/finite_evidence_coverage_matrix.tsv: 配体/受体RNA覆盖，不是功能验证。
results/treg_original_fraction_residual_side_by_side.tsv: 原七基因TregNoTF、Treg比例、两基因残余面板及TF RNA/AUC并列。两基因仅受限探索。
results/reused_*: 原历史统计原样子集，不变更旧FDR。reused_mouse_target_program/鼠相关表不作为人类独立验证。
code/finite_plan.txt: 结果前冻结的185关系与统计边界。
code/extract_audit.py: 原计数/冻结注释提取；拒绝覆盖已完成提取。
code/duplicate_audit.py: 新旧Treg/Activation评分相等核对。
code/finite_models.py: 有限样本级模型，70规格复用旧统计；拒绝覆盖已完成运行。
code/verify_results.py: 独立偏秩统计、BH、LOO、固定供体和输入哈希核验。
code/build_delivery.py: 报告和打包。
复现依赖服务器原输入，路径见results/input_manifest.tsv；Python环境为 /mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python。
为保护旧结果，代码有完成标记保护。若需重新执行，应复制代码到新的目录并统一修改O路径；不能删除既有完成标记后直接覆盖。
所有source/phenotype关系都是探索性关联，无实验因果/中介主张。≥20与NoDom为预定敏感性，不挑选阈值。
'''
(O/'README.txt').write_text(readme,encoding='utf-8')
completion=dict(status='COMPLETE_VERIFIED',time=datetime.datetime.now().isoformat(),models=1675,hypotheses=185,reused=70,group_context_models=125,old_models_unchanged=True,report='lung_front_half_report_zh.html')
(O/'provenance/complete.json').write_text(json.dumps(completion,indent=2));(O/'STATE.txt').write_text('COMPLETE_VERIFIED\n')
manifest=[]
for pth in sorted(O.rglob('*')):
 if pth.is_file() and pth.suffix not in ['.zip','.pyc'] and pth.name not in ['manifest.tsv','analysis.lock','delivery.log']:manifest.append(dict(path=str(pth.relative_to(O)),bytes=pth.stat().st_size,sha256=hashlib.sha256(pth.read_bytes()).hexdigest()))
pd.DataFrame(manifest).to_csv(O/'manifest.tsv',sep='\t',index=False)
zipname=O/'lung_front_half_results_code.zip'
with zipfile.ZipFile(zipname,'w',zipfile.ZIP_DEFLATED) as z:
 for pth in sorted(O.rglob('*')):
  if pth.is_file() and pth.suffix not in ['.zip','.pyc'] and pth.name not in ['analysis.lock','delivery.log']:z.write(pth,str(pth.relative_to(O)))
with zipfile.ZipFile(zipname) as z:assert z.testzip() is None
print(json.dumps(dict(completion=completion,report_bytes=(O/'lung_front_half_report_zh.html').stat().st_size,zip_bytes=zipname.stat().st_size,zip_sha256=hashlib.sha256(zipname.read_bytes()).hexdigest()),indent=2))
