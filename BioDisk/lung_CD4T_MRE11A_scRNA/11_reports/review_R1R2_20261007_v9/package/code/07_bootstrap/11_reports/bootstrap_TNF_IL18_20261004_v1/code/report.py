"""Self-contained Chinese HTML, export figures and reproducible archive."""
from pathlib import Path
import json,hashlib,zipfile,html,io,platform
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
B=Path(__file__).resolve().parents[1];R=B/'results';D=B/'deliverables';D.mkdir(exist_ok=True)
def rd(n):return pd.read_csv(R/n,sep='\t',low_memory=False)
s=rd('effect_summary.tsv');ds=rd('paired_difference_summary.tsv');ch=rd('chain_summary.tsv');dg=rd('diagnostic_summary.tsv');alls=rd('all_seeds_chain_summary.tsv')
plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','DejaVu Sans'],'axes.unicode_minus':False,'svg.fonttype':'path','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
MOD=['RNA19','RNA17','AUC17','AUC29','AUC43'];EDGE=['ligand_TF','TF_phenotype','ligand_phenotype'];EN={'ligand_TF':'配体—TF','TF_phenotype':'TF—效应表型','ligand_phenotype':'配体—效应表型'}
COL={'RNA19':'#65778b','RNA17':'#a57149','AUC17':'#176c79','AUC29':'#52868e','AUC43':'#a1bfc2'}
def savefig(fig,name):
 fig.savefig(D/(name+'.png'),dpi=170,bbox_inches='tight',facecolor='white');fig.savefig(D/(name+'.svg'),bbox_inches='tight',facecolor='white');plt.close(fig)
def forest(axis,adjustments,name):
 fig,axs=plt.subplots(1,3,figsize=(13,5.4 if len(adjustments)==1 else 7),sharey=True)
 labels=[f'{m} · {a}' for a in adjustments for m in MOD]
 for ax,e in zip(axs,EDGE):
  for y,(adj,m) in enumerate((a,m) for a in adjustments for m in MOD):
   q=s[(s.scheme=='stratified')&(s.axis==axis)&(s.adjustment==adj)&(s.modality==m)&(s.edge==e)].iloc[0]
   ax.plot([q.CI95_low,q.CI95_high],[y,y],color=COL[m],lw=2.5);ax.plot(q.original_effect,y,'o',color=COL[m],ms=6);ax.plot(q.bootstrap_median,y,'|',color='#172b36',ms=11)
  ax.axvline(0,color='#aab4bb',lw=1,ls='--');ax.set_xlim(-1.03,1.03);ax.set_xticks([-1,-.5,0,.5,1]);ax.set_title(EN[e]);ax.set_xlabel('秩相关 ρ');ax.grid(axis='x',alpha=.15)
 axs[0].set_yticks(range(len(labels)),labels);axs[0].invert_yaxis()
 fig.suptitle(f'{axis}：分层主分析，5000 次供体重抽样',fontsize=16,x=.51,y=.99);fig.text(.51,.01,'圆点：原始效应  ｜  竖线：bootstrap 中位数  ｜  横线：95% 百分位区间',ha='center',color='#526570');fig.tight_layout(rect=[0,.04,1,.95]);savefig(fig,name)
forest('IL18',['B'],'IL18_effects');forest('TNF',['B','B_Treg'],'TNF_effects')
fig,axs=plt.subplots(1,3,figsize=(13,4.5),sharey=True)
for ax,e in zip(axs,EDGE):
 for i,m in enumerate(MOD):
  q=ds[(ds.scheme=='stratified')&(ds.comparison_type=='Treg_enhancement')&(ds.modality==m)&(ds.edge==e)].iloc[0]
  ax.plot([q.CI95_low,q.CI95_high],[i,i],color=COL[m],lw=2.5);ax.plot(q.original_effect,i,'o',color=COL[m]);ax.plot(q.bootstrap_median,i,'|',color='#172b36',ms=11)
 ax.axvline(0,color='#aab4bb',ls='--');ax.set_xlim(-1.1,1.2);ax.set_title(EN[e]);ax.set_xlabel('Δρ = ρ(B + Treg) − ρ(B)');ax.grid(axis='x',alpha=.15)
axs[0].set_yticks(range(5),MOD);axs[0].invert_yaxis();fig.suptitle('TNF：同次供体抽样的 Treg 调整差值',fontsize=16);fig.tight_layout(rect=[0,0,1,.93]);savefig(fig,'TNF_Treg_delta')
fig,axs=plt.subplots(1,2,figsize=(12,4.4),sharey=True)
for ax,e in zip(axs,EDGE[:2]):
 for i,c in enumerate(['AUC29-AUC17','AUC43-AUC17','AUC43-AUC29']):
  q=ds[(ds.scheme=='stratified')&(ds.axis=='IL18')&(ds.adjustment=='B')&(ds.edge==e)&(ds.comparison==c)].iloc[0]
  ax.plot([q.CI95_low,q.CI95_high],[i,i],lw=3,color='#176c79');ax.plot(q.original_effect,i,'o',color='#176c79');ax.plot(q.bootstrap_median,i,'|',color='#172b36',ms=11)
 ax.axvline(0,color='#aab4bb',ls='--');ax.set_xlim(-.9,.6);ax.set_title(EN[e]);ax.set_xlabel('同次配对 Δρ');ax.grid(axis='x',alpha=.15)
axs[0].set_yticks(range(3),['seed29 − seed17','seed43 − seed17','seed43 − seed29']);axs[0].invert_yaxis();fig.suptitle('IL18：基础调整后的 seed 配对差值',fontsize=16);fig.tight_layout(rect=[0,0,1,.92]);savefig(fig,'IL18_seed_delta')
rename={'axis':'轴','modality':'模态','adjustment':'调整','edge':'边','original_effect':'原效应','bootstrap_median':'重抽中位数','CI95_low':'95%下限','CI95_high':'95%上限','sign_preservation':'符号保持比例','valid':'有效次数','failure_rate':'失败率','positive_fraction':'正值比例','scheme':'抽样','comparison':'配对差','direction_compatible_fraction':'三边方向相容','all_three_original_sign_fraction':'三边同保原符号','failed':'失败次数','roster':'名单','unique_donors_median':'独特供体中位数','unique_donors_min':'独特供体最少','multiplicity_effective_n_median':'重复权重有效数中位数','low_unique_support_fraction':'支持df<3比例','rank_deficient':'秩亏次数','unique_support_df_min':'支持df最少','nominal_df_min':'名义df最少','compatible_all_attempts_lower':'全尝试相容下界','compatible_all_attempts_upper':'全尝试相容上界','all_seeds_compatible_fraction':'三seed均相容比例','all_seeds_original_sign_fraction':'三seed均保原符号'}
def tbl(df,cols=None):
 if cols is not None:df=df[cols]
 df=df.copy()
 if 'edge' in df:df['edge']=df.edge.map(EN)
 return '<div class="scroll">'+df.rename(columns=rename).to_html(index=False,border=0,na_rep='NA',float_format=lambda x:f'{x:.4f}')+'</div>'
def svg(name):
 text=(D/(name+'.svg')).read_text(encoding='utf-8');return '<figure>'+text[text.index('<svg'):]+'</figure>'
def p(t):return '<p>'+t+'</p>'
def get(axis,mod,edge,adj='B',scheme='stratified'):
 return s[(s.axis==axis)&(s.modality==mod)&(s.edge==edge)&(s.adjustment==adj)&(s.scheme==scheme)].iloc[0]
def fmt(q):return f'{q.original_effect:.3f}，95%区间 [{q.CI95_low:.3f}, {q.CI95_high:.3f}]，符号保持 {q.sign_preservation:.1%}'
ilfirst=[get('IL18',m,'ligand_TF') for m in MOD[2:]];ilsecond=[get('IL18',m,'TF_phenotype') for m in MOD[2:]];illast=get('IL18','AUC17','ligand_phenotype')
parts=['''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>TNF / IL18 供体配对 bootstrap</title><style>
:root{--ink:#243947;--teal:#176c79;--muted:#62717b}*{box-sizing:border-box}body{margin:0;background:#eef2f4;color:var(--ink);font:16px/1.8 "Microsoft YaHei",sans-serif}main{max-width:1280px;margin:28px auto;background:white;padding:38px 44px;border-top:7px solid var(--teal)}h1{font-size:32px;line-height:1.4}h2{font-size:24px;margin-top:38px;border-bottom:1px solid #d4dfe3;padding-bottom:8px}h3{font-size:19px}small,.note{color:var(--muted)}.lead{background:#edf6f6;border-left:5px solid var(--teal);padding:17px 24px}.warn{background:#fcf6e9;border-left:4px solid #bc8f39;padding:14px 20px}.grid{display:grid;grid-template-columns:repeat(3,1fr);gap:12px}.card{background:#f1f5f7;padding:13px 18px}.card strong{font-size:26px;color:var(--teal)}figure{margin:20px 0}figure svg{width:100%;height:auto}.scroll{overflow-x:auto;margin:16px 0}table{border-collapse:collapse;font-size:13px;width:100%;font-variant-numeric:tabular-nums}th,td{padding:8px 10px;border-bottom:1px solid #dde5e9;text-align:left;white-space:nowrap}th{background:#eaf1f4}tr:nth-child(even){background:#f7f9fa}details{margin:20px 0}summary{cursor:pointer;color:var(--teal);font-weight:bold}pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:13px;background:#f5f7f8;padding:16px}code{overflow-wrap:anywhere}a{color:var(--teal)}@media(max-width:700px){main{padding:18px;margin:0}h1{font-size:25px}.grid{grid-template-columns:1fr}}@media print{body{background:white}main{margin:0;padding:12px}details{display:block}.scroll{overflow:visible}th,td{white-space:normal;font-size:10px;padding:4px}figure{break-inside:avoid}}</style></head><body><main>
<small>2026-10-04 · 冻结输入 · 供体整行配对 · 探索性抽样稳定性</small><h1>TNF–STAT1 与 IL18–IKZF2<br>供体级配对 bootstrap</h1>
<div class="grid"><div class="card"><strong>19 / 17</strong><br>RNA / AUC 独立供体</div><div class="card"><strong>5000 次</strong><br>每方案、每固定名单</div><div class="card"><strong>3 个 seed</strong><br>同一次供体索引配对</div></div>''']
parts+=['<div class="lead"><b>结论先读</b>'+p('两条线均保留为候选，但稳定性集中在不同位置。IL18 的 IKZF2—效应表型负向关系比首边稳定；首边对 seed 和供体抽样更敏感，末边精度仍有限。TNF 的正向结构较一致，加 Treg 后的相关与“增加了多少”必须分开判断；增强差值的区间比原始点差更能显示不确定性。')+p('这次没有选择最佳 seed、换基因或重建 regulon。区间和方向比例只描述冻结分数在当前供体重抽样下的行为，不是 P 值、真实方向概率或因果概率。')+'</div>']
parts+=['<h2>1. 冻结范围与可估计性</h2>',p('两条线分别为 MacAM→CD4-IL7R：IL18–IKZF2，及 MacAM→Core：TNF–STAT1。唯一效应终点为冻结 State_Effector。每条线只计算配体—TF、TF—效应表型、配体—效应表型三边。RNA19、RNA17 与 AUC17/29/43 并列；RNA17 是明确匹配 AUC 的同一 17 人，AUC17 指 seed17，并非额外的名单。'),p('B=疾病组 G + 来源 Source + 两端 library 深度。连续端点、深度和 Treg 比例每次重新计算平均秩；对固定设计求残差相关。未调整、B 均计算；TNF 另有 B+Treg。共 75 个展示规格，主分析与敏感性各一套；重复末边保留用于三边配对展示，不能当独立证据。'),tbl(rd('design_cells.tsv'),['roster','series','group','n']),p('原样本所有 75 项均可估计，69 项与旧表直接对应，差值小于 10⁻¹⁰；6 项为预定 TNF 加 Treg 的 RNA17/AUC29/AUC43 前两边。未调整秩=1；B 秩=6，RNA19/17 名义 df=12/10；B+Treg 秩=7，df=11/9。RNA19→17 排除 CTL070、CTL103 是冻结上游 AUC 不可用所致，没有按本次效应筛人。')]
parts+=['<h2>2. 抽样与区间为什么这样定</h2>',p('主分析按“来源×疾病组”保留原人数、层内有放回抽整名供体。两条线、三个 seed、同名单 RNA/AUC、所有调整和三边共享每次索引；RNA19 单独抽样，不与 17 人结果构造伪配对差。预定唯一抽样敏感性为不分层整供体抽样，同样 5000 次。'),'<div class="warn">AUC 名单的 G135 仅 CTL068 一名对照，分层抽样中此人固定。G227 仅 2 名病例，G288 全为病例。因此主区间条件于当前来源/疾病构成，可能低估稀疏单元的总体不确定性；不分层敏感性用于显示这种限制，不能以较窄区间取胜。</div>',p('预定百分位区间：有效重复效应的 2.5%、97.5% 分位数，线性插值；不试选 BCa、Fisher 或其他区间。失败不补抽、不删协变量挽救。每项列出有效分母和失败率；区间以该项可估计为条件。小样本秩统计、重复权重与近退化设计可使区间覆盖失准，尤其不能把 5000 次当成 5000 名患者。'),p('保存每次独特供体数 u、重复权重有效数 n²/Σm²、名义 df=n−设计秩−1，以及支持诊断 u−设计秩−1。后者不是校准推断 df；低于 3 时明确标记但不事后剔除。数学上仍可计算的边界相关保留。秩亏、端点/残差零方差、名义 df<3 或非有限输入均记录为失败。')]
parts+=['<h2>3. IL18：跨 seed 的哪一段更稳？</h2>',svg('IL18_effects'),p('B 调整首边：'+ '；'.join(f'seed{m[3:]}：{fmt(q)}' for m,q in zip(MOD[2:],ilfirst))+'。'),p('B 调整第二边：'+ '；'.join(f'seed{m[3:]}：{fmt(q)}' for m,q in zip(MOD[2:],ilsecond))+'。'),p('末边不含 TF/AUC，三个 seed 在同供体名单上完全相同：'+fmt(illast)+'。它不是三个独立 seed 的重复验证。'),svg('IL18_seed_delta'),p('seed 差必须由同一次抽样的两个相关相减；并列区间是否重叠不能代替该差值区间。以下配对差固定用后者减前者，不选择最大效应 seed。'),tbl(ds[(ds.scheme=='stratified')&(ds.axis=='IL18')&(ds.adjustment=='B')&(ds.comparison_type=='seed_paired')],['comparison','edge','original_effect','bootstrap_median','CI95_low','CI95_high','sign_preservation','valid','failure_rate']),p('这里比较的是三个已冻结的计算产物在相同供体扰动下的差异。靶集合和 AUCell 排名生成的随机性没有被重新抽取，不能由本次 bootstrap 单独量化或归因。')]
parts+=['<h2>4. TNF：正向关联与 Treg 增强分开看</h2>',svg('TNF_effects'),svg('TNF_Treg_delta'),p('Δρ=ρ(B+Treg)−ρ(B)。正差只表示相关数值向正方向增加；不是 Treg 的生物学抑制/促进或中介效应。每次两种调整使用同一供体索引，差值仅在双方都可估计时计算。'),tbl(ds[(ds.scheme=='stratified')&(ds.comparison_type=='Treg_enhancement')],['modality','edge','original_effect','bootstrap_median','CI95_low','CI95_high','positive_fraction','valid','failure_rate']),p('此前冻结诊断中的负残差协方差项和方差归一化，可以解释原始点估计为何增加；本次则检查这种增加对供体抽样有多稳定。即使调整后相关区间偏正，也不能据此推断增强差值的区间一定大于零。Core 细胞组成与细胞内状态仍未分离。')]
parts+=['<h2>5. 三边同次方向相容，不是三次独立支持</h2>',p('方向相容定义：sign(配体—TF)×sign(TF—效应)=sign(配体—效应)，三边非零且均有效；它允许相对于原图的多条边同时翻转。另列“三边同保原符号”，以区分两者。以下分母为同次三边都有效的次数。'),tbl(ch[(ch.scheme=='stratified')&ch.adjustment.isin(['B','B_Treg'])],['axis','modality','adjustment','valid','direction_compatible_fraction','all_three_original_sign_fraction','failure_rate']),p('更严格地要求三个 seed 在同一次抽样都满足三边相容：'),tbl(alls,['scheme','axis','adjustment','valid','all_seeds_compatible_fraction','all_seeds_original_sign_fraction','compatible_all_attempts_lower','compatible_all_attempts_upper']),p('全尝试上下界把失败重复分别视作全不相容/全相容，展示失败造成的信息缺口；不把失败算作已证实方向。代数相容不估计间接效应，也不建立配体→TF→表型因果链。')]
parts+=['<h2>6. 秩亏、供体多样性与不分层敏感性</h2>',tbl(dg,['scheme','axis','roster','adjustment','rank_deficient','unique_donors_min','unique_donors_median','multiplicity_effective_n_median','unique_support_df_min','low_unique_support_fraction']),p('不分层抽样允许稀疏单元消失；固定 G 与 Source 列不会随之被删掉。其可估计重复因此是被选择后的子集，高失败率使条件区间不能被视作充分解决来源问题。全部尝试和原因组合保留在逐次结果中，零方差、残差退化及条件数见完整诊断表。'),tbl(s[(s.scheme=='unstratified')&s.adjustment.isin(['B','B_Treg'])&s.modality.isin(['RNA19','AUC17'])],['axis','modality','adjustment','edge','original_effect','bootstrap_median','CI95_low','CI95_high','sign_preservation','valid','failure_rate'])]
parts+=['<h2>7. 明确同供体的 RNA/AUC 比较</h2>',p('下面只比较匹配 17 人、同一次索引的 AUC−RNA17 相关差。RNA19 保留为全名单分析，但其与 AUC17 的差同时涉及名单变化，不能解读为纯模态差异。RNA 和 regulon AUC 测量对象不同，配对相关差本身也不是测量优劣。'),tbl(ds[(ds.scheme=='stratified')&(ds.comparison_type=='AUC_minus_matched_RNA')&(ds.adjustment=='B')],['axis','comparison','edge','original_effect','bootstrap_median','CI95_low','CI95_high','valid','failure_rate'])]
parts+=['<h2>8. 结论边界与可复现交付</h2>',p('本次可支持的是有限两条线在冻结供体表上的抽样稳定性描述。IL18 保留 seed 相关的首边不确定性；TNF 保留组成/调整敏感性。没有引入损伤焦亡前置、候选扩展或新数据。若要区分细胞内机制、来源混杂和 regulon 生成随机性，需要独立于本次重抽样的新测量或专门设计；本报告不据此追加分析。'),p('主估计代码和方案在重抽样前已记录哈希；四组索引先保存再估计。原始效应一致性、每次分层人数、三个 seed 的配对索引、非 TF 末边恒等性，以及选定重复的独立 Spearman / QR 残差实现均已核验。未计算新 P/FDR 或“显著比例”。'),'<details><summary>查看完整 150 行单边汇总（含未调整及两种抽样）</summary>'+tbl(s,['scheme','axis','modality','adjustment','edge','original_effect','bootstrap_median','CI95_low','CI95_high','sign_preservation','valid','failure_rate','abs_rho_boundary_fraction'])+'</details>','<details><summary>完整配对差汇总（含 TNF 三 seed）</summary>'+tbl(ds,['scheme','axis','modality','adjustment','comparison','edge','original_effect','bootstrap_median','CI95_low','CI95_high','positive_fraction','valid','failure_rate'])+'</details>','<details><summary>完整可估计性诊断汇总</summary>'+tbl(dg)+'</details>','<details><summary>运行前冻结方案原文</summary><pre>'+html.escape((B/'FROZEN_PLAN.txt').read_text(encoding='utf-8'))+'</pre></details>','<details><summary>核验记录</summary><pre>'+html.escape((B/'provenance/VERIFIED.json').read_text())+'</pre></details>',p('ZIP 含 inputs/、code/、results/、provenance/ 和图表。完整逐次效应 750,000 行、设计诊断 100,000 行；重复展示末边不算独立生物学证据。随机数为 NumPy PCG64，4 个种子为 2026100401–2026100404。复现请先阅读 README.txt；不要运行旧目录里的历史模型脚本。'),'<p class="note">数据和所有结论来源：用户授权项目的两套冻结供体输入、旧结果与方法代码，以及本次预定 bootstrap。报告不新增临床或因果推断。</p></main></body></html>']
parts.insert(2,'<div class="lead">主分析的关键实数：IL18 三 seed 的首边区间均跨零，第二边均为负区间，末边跨零；首边 seed 配对差区间也均跨零。TNF 五种模态展示的三边 Treg 增强差值，15 个区间全部跨零。它们保留了方向倾向，但不能据此断言稳定增强或 seed17 优于其他 seed。</div>')
parts.insert(-8,p('失败实数：分层主分析无秩亏；TNF 的 17 人 B+Treg 有 2/5000 次零残差失败（0.04%），其余主规格均有效。不分层 B 的 RNA19 秩亏 1207/5000（24.14%），17 人秩亏 2192/5000（43.84%）；B+Treg 的 17 人秩亏 2193 次，另需结合零残差失败。分层独特供体中位数为 RNA19 的 13 人、17 人名单的 12 人，重复权重有效数中位数约 10.94 / 9.97。17 人 B+Treg 中 16.44% 的尝试支持诊断 df<3，均已显式标记。'))
checks=rd('original_legacy_checks.tsv')
legacy=pd.read_csv(B/'inputs/aligned_chain_models.tsv',sep='\t')
legacy_ids=set(';'.join(checks.legacy_ids.dropna()).split(';'))
reused=legacy[legacy.spec_id.isin(legacy_ids)].copy()
reused.to_csv(R/'legacy_statistics_REUSED.tsv',sep='\t',index=False,na_rep='NA')
insert_at=parts.index('<h2>8. 结论边界与可复现交付</h2>')+1
parts.insert(insert_at,p('旧单次 P/FDR 完整保留，bootstrap 不替代或改写既有单次推断：inputs/aligned_chain_models.tsv 原文件（含全部旧规格的 P、FDR、历史 FDR、family 与计算方法）原样保存。本次对应原始效应的旧结果另存 legacy_statistics_REUSED.tsv 并在下表展示；6 个无直接旧规格的 Treg 补充项不生成新 P/FDR。旧 P/FDR 保持各自原有校正 family，不能改当这75个展示规格的新校正结果。'))
parts.insert(insert_at+1,'<details><summary>对应规格的旧单次 P/FDR（原样复用，非 bootstrap 推断）</summary>'+tbl(reused,['spec_id','chain','modality','n','edge','x','y','covariates','rho','P','FDR','legacy_FDR','FDR_family','P_method'])+'</details>')
insert_at=parts.index('<h2>2. 抽样与区间为什么这样定</h2>')+1
parts.insert(insert_at,p('分层单元人数按 G135-Control、G227-Control、G227-CWP、G288-CWP 排列：RNA19 为 2、8、2、7；匹配17人为 1、7、2、7。其余来源×疾病组合在原数据中不存在。AUC 单人层 CTL068 每次恰出现1次；两人层只在这两名已观测供体间重抽。分层不能制造独立来源验证，不分层可估计子集也不能支持无条件稳健结论。'))
parts.insert(insert_at+1,p('沿用旧定义，seed17 为主 seed，seed29/43 为预定并列敏感性；主 seed 的保留不是依据本次结果选择。全部 seed 的原效应、区间、配对差和失败率都报告。'))
brief=[]
for axis,adj in [('IL18','B'),('TNF','B'),('TNF','B_Treg')]:
 for m in MOD:
  q=s[(s.scheme=='stratified')&(s.axis==axis)&(s.adjustment==adj)&(s.modality==m)].set_index('edge')
  c=ch[(ch.scheme=='stratified')&(ch.axis==axis)&(ch.adjustment==adj)&(ch.modality==m)].iloc[0]
  rec={'轴':axis,'模态':m,'调整':adj}
  for e in EDGE:
   r=q.loc[e];rec[EN[e]+' 95%CI']=f'[{r.CI95_low:.3f}, {r.CI95_high:.3f}]'
  rec['三边均有效率']=f'{c.valid/5000:.2%}';rec['同次三边相容']=f'{c.direction_compatible_fraction:.2%}'
  brief.append(rec)
parts.insert(3,'<details open><summary>调整结果摘要：逐边区间、有效率与联合方向（分层主分析）</summary>'+tbl(pd.DataFrame(brief))+'</details>')
(D/'TNF_IL18_donor_bootstrap_zh.html').write_text('\n'.join(parts),encoding='utf-8')
(B/'README.txt').write_text('''TNF/IL18 供体配对 bootstrap；先阅读 FROZEN_PLAN.txt。
结果：deliverables/TNF_IL18_donor_bootstrap_zh.html 为自包含中文图文报告。
代码：bootstrap.py --preflight 冻结原始名单/效应/输入哈希；bootstrap.py --run 使用既定索引规则执行；summarize.py 汇总并核验；report.py 生成报告和图表。
原估计环境：Python 3.11.15, NumPy 2.4.6, pandas 2.3.3, SciPy 1.17.1；绘图可用 matplotlib。
现有目录含完成标记与锁，阻止误重跑。要独立复现，在另一个全新目录仅复制 inputs/、code/、FROZEN_PLAN.txt，再按上述顺序运行。不要删当前目录的完成标记或覆盖结果。
所有抽样单位为供体。RNA19独立索引；RNA17和AUC17/29/43及两条线、三边、调整共享17人索引。indices*.npy/tsv.gz为0-based索引，roster*.tsv给定排序。逐次结果失败行保留NA和原因，不补抽、不删协变量。
percentile CI仅在共同有效重复上汇总，显式报告分母；不分层高失败率需结合诊断理解。方向比例非P、非因果概率。重复末边与三个seed不是独立生物重复。
输入来源：/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports 下 breakpoint_TNF_IL18_20261004_v1 与 dual_line_TNF_IL18_20261004_v1；旧目录只读。
''',encoding='utf-8')
print(json.dumps({'report':str(D/'TNF_IL18_donor_bootstrap_zh.html'),'figures':4,'bytes':(D/'TNF_IL18_donor_bootstrap_zh.html').stat().st_size}))
