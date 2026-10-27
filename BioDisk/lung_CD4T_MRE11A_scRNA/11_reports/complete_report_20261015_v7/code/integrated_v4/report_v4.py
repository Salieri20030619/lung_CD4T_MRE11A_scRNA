from pathlib import Path
import sys,json,hashlib,html
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(sys.argv[1]);rd=lambda p:pd.read_csv(O/p,sep='\t')
f=rd('03_receiver/IKZF2_STAT1_primary_and_sensitivity.tsv'); print('FOCAL_COLUMNS',list(f.columns))
p=rd('05_human/reviewed_sample_scope_roster.tsv');print('HUMAN_COLUMNS',list(p.columns))
fig,axes=plt.subplots(2,1,figsize=(11,6),layout='constrained')
mods=['PyroExec','Inflam','ApoExec','OXPHOS','MitoISR','MitoROS','Mitophagy'];states=['Core','Naive','ActCD4','Treg']
for ax,tf in zip(axes,['Ikzf2(+)','Stat1(+)']):
 b=f[f.x==tf].drop_duplicates(['scope','y']); mat=b.pivot(index='scope',columns='y',values='rho').reindex(index=states,columns=mods);im=ax.imshow(mat,vmin=-1,vmax=1,cmap='RdBu_r',aspect='auto');ax.set_xticks(range(7),mods);ax.set_yticks(range(4),states);ax.set_title(tf+' activity vs program (6 mice; seed 17)')
 for i,s in enumerate(states):
  for j,m in enumerate(mods):
   z=b[(b.scope==s)&(b.y==m)]
   if len(z):ax.text(j,i,f'{z.iloc[0].rho:.2f}'+('*' if z.iloc[0].P<.05 else ''),ha='center',va='center',color='black',fontsize=9)
fig.colorbar(im,ax=axes,label='Spearman rho',shrink=.8);fig.suptitle('* nominal exact P < 0.05; focal hits do not pass BH < 0.05',fontsize=11);fig.savefig(O/'figures/TF_programs.png',dpi=170);plt.close(fig)
c=rd('06_support/Core_fixed_composition_focal.tsv');fig,ax=plt.subplots(figsize=(8,4.3),layout='constrained')
for i,(tf,m) in enumerate([('Ikzf2(+)','ApoExec'),('Stat1(+)','MitoROS')]):
 orig=f[(f.scope=='Core')&(f.x==tf)&(f.y==m)].iloc[0];fix=c[(c.regulon==tf)&(c.module==m)&(c.model=='Pooled')].iloc[0]
 ax.plot([0,1],[orig.rho,fix.rho],'o-',label=tf+' / '+m,lw=2)
 for x,z in enumerate([orig,fix]):ax.annotate(f'rho={z.rho:.3f}; P={z.P:.4f}; BH={z.BH:.3f}',(x,z.rho),xytext=(5,10 if i==0 else -19),textcoords='offset points',fontsize=8)
ax.axhline(0,color='grey',lw=.8);ax.set_xlim(-.12,1.65);ax.set_ylim(-1.15,1.25);ax.set_xticks([0,1],['Observed Core composition','Fixed subtype composition']);ax.set_ylabel('Spearman rho');ax.set_title('IKZF2 association changes sign after composition adjustment');ax.legend(loc='lower right');fig.savefig(O/'figures/composition_sensitivity.png',dpi=170);plt.close(fig)
# Internal working report: no customer translation or causal relabeling.
def table(path,cols=None,limit=None):
 d=rd(path)
 if cols:d=d[[x for x in cols if x in d.columns]]
 if limit:d=d.head(limit)
 return d.to_html(index=False,escape=True,float_format=lambda x:f'{x:.5g}',na_rep='NA')
def link(p,label=None):return '<a href="'+p+'">'+html.escape(label or p)+'</a>'
body='''<h1>矽肺巨噬细胞—CD4 整合路径：v4 内部分析记录</h1><p>2026-10-15 · 独立修订目录 · 保留既有结果 · 尚非客户交付稿</p>
<h2>主框架可以保留，但完整机制链尚未得到支持</h2><p>本轮围绕“巨噬细胞信号 → CD4 调控/响应 → 线粒体及死亡相关转录程序”重新组织既有数据。状态内拟时序和小鼠真实时间均作为辅助证据，不作为候选准入硬门槛，也不要求历史六基因来自自动交集。尚未形成原三配体→IKZF2/STAT1→焦亡程序的双段名义显著链。</p>
<p><strong>STAT1 可保留为优先的 CD4 接收端候选；IKZF2 保留为细胞状态/组成相关调控线索。</strong>这不等同于已经证明线粒体损伤或焦亡，更未证明巨噬细胞对其因果作用。</p>
<h2>STAT1 的线粒体程序关联较稳，IKZF2 的 Core 关联受组成影响</h2><ul><li>Core IKZF2 活性与 ApoExec：rho=1，P=0.002778，BH=0.2743；固定亚群组成后 rho=-0.7714，P=0.1028。不能将原正相关解释为状态内促凋亡作用。</li><li>Core STAT1 活性与 MitoROS：rho=0.9429，P=0.01667，BH=0.3678；固定组成后 rho=0.8857，P=0.03333，BH=0.4467。方向保留；控制疾病组别后原关联 P=0.2778，固定组成后 P=0.1667。</li><li>Treg STAT1–MitoROS：rho=0.8857，P=0.03333，BH=0.3790；种子43不再名义显著。上述关联均为探索性证据。</li><li>IKZF2/STAT1 活性与 PyroExec 未见名义显著关联。MitoROS 是相关转录程序，不是直接的线粒体功能损伤测量。</li></ul>
<img src="figures/TF_programs.png"><img src="figures/composition_sensitivity.png">
<h2>TNF、IL15 的预测排序可复用，上游配对证据仍不足</h2><p>AMlike 来源限定后，Core 上调靶基因预测中 TNF 排第3、IL15 第5；ActCD4 中分别第4、第6。有表达/受体支持，但 AMlike TNF 和 IL15 疾病差异原始 P 分别为0.3230和0.3588，不能称为疾病上调信号；也没有得到它们与 STAT1 活性的名义显著配对关联。预测排序无生物学效应 P 值。</p>
<p>AMlike IL18 下调：logFC=-0.7132，P=0.004049，BH=0.1361；AllMac 下调更强，logFC=-2.0459，P=0.000106，BH=0.007101。AMlike IL18 与 ActCD4 IKZF2 活性负相关：rho=-0.9429，P=0.01667，BH=0.65；其他种子方向一致。这是一段负向关联，并未接上相同状态/相同活性读数的焦亡或线粒体表型链。</p>
<p>IFNG 预测排名靠前，但病例巨噬细胞中的表达极稀疏，来源审计不足以支持其作为主要巨噬细胞信号。不能凭配体排名跳过来源核查。表达共现、配体靶基因先验、TF 活性和基因 RNA 是不同证据，不跨读数拼接机制链。</p>
<h2>原基因保留各自角色，不为保留目标调整阈值</h2><p>原六候选在 ActCD4 中仍全部达到原始 P&lt;0.05；其中 IL1R2 的 BH=0.04044。其余候选的完整效应、P 与 BH 见下表。它们保留为差异响应/验证指标，不能因为有差异就自动升级成上游调控靶点。STAT1 在 Treg 的 RNA 差异为 logFC=2.1506，P=2.5823×10⁻⁵，BH=0.005589；IKZF2 在 Core 名义显著，但 ActCD4/Treg RNA 差异不显著。RNA 差异和 TF 活性证据分别解释。</p>
<p>AHNAK、MDFIC、TTC39C、ST6GALNAC3 保留为补充候选，之前交集结果不再强制替代历史候选。创新性只能用于同等证据下的排序，不能修复统计或方向缺口；既有创新性与方向矛盾记录继续有效。</p>'''
body+=table('03_receiver/candidate_disease_DE.tsv')
body+='''<h2>患者证据定位为有限支持和细胞定位</h2><p>GSE174725 复核注释后，病例中仅1个样本具有 BALF_AM。未设置20或10细胞门槛，仍无法进行有效的人类巨噬细胞组间或配对推断。不以单核细胞代替肺泡巨噬细胞。CD4 程序组间比较未见原始 P&lt;0.05；人类 IKZF2/STAT1 RNA 亦未达到名义显著。不能把这些结果写成患者层面的机制验证。</p><p>GSE296682 仅提供患者细胞定位，保留原作者标签；不作病例对照推断。GSE288887 煤工尘肺退出主发现及选靶，不与矽肺患者混合推断。本轮未下载 PRJNA1289989。</p>'''
body+=table('05_human/reviewed_sample_scope_roster.tsv')
body+='''<h2>建议的基础实验顺序</h2><ol><li>先验证表型：在明确 CD4 亚群中分别测线粒体功能/应激和死亡类型；线粒体膜电位、线粒体 ROS 与细胞存活读数配套，焦亡需裂解型 gasdermin、相应蛋白酶及膜通透性等独立证据。转录评分不足以确定死亡类型。</li><li>再验证细胞间影响：设置巨噬细胞条件培养/共培养、对照和硅尘直接作用对照，区分巨噬细胞介导与残留颗粒直接效应；记录 T 细胞亚群比例和亚群内表型。</li><li>表型建立后验证接收端：STAT1 优先作为受体侧候选，配合活化及扰动读数；IKZF2 需在相同状态内验证，避免把亚群比例改变当作细胞内作用。</li><li>最后定位上游：TNF/IL15 作为候选分别阻断或补充，并设置适当对照及救援。IL18 按负向/补偿假设独立考察。不得预设阻断方向或必然改善焦亡。</li></ol><p>这是实验优先级建议，不是已经确认的药理靶点。若目标是较少改动旧报告，可保留背景、差异结果和通讯框架；必须改正组成敏感的 IKZF2 解释、未闭合链条及“已发生焦亡”的措辞。</p>
<h2>统计方法与复核范围</h2><p>发现队列 GSE314943：3 saline + 3 silica，20周。以独立小鼠为单位，不将细胞或随机种子当重复。所有非空可计算样本/状态保留；冻结差异表及其原 BH 复用。TF 主读数为种子17的 DropOverlap 活性，29/43为敏感性；先验网络不代表实测因果。主要关联为6样本精确 Spearman 置换（720），疾病组别调整采用组内置换（36）。BH 在完整同类假设家族中计算，包括全部评估 TF，不仅对选中的两个 TF 校正。原始 P&lt;0.05 为名义筛选，完整 P、BH、效应和检验家族均保留。</p><p>固定组成分析以全部6只小鼠的平均亚群比例为统一权重，重新组合各状态活性和程序。真实时间 GSE274372 D0/D7/D28/D56 与细胞拟时序分开；GSE183682 D7/D56 每条件/时间仅1只小鼠，只作方向支持，且不能与20周发现终点等同。</p><p>独立数值核查通过：1576项焦点 P 值和100个 BH 家族；最大 BH 数值误差约1.58×10⁻¹⁴。复核不等于新增独立生物学验证。小样本、疾病组间分离、程序与活性签名剩余重叠、来源表达稀疏和人类样本覆盖限制共同约束解释。DropOverlap 去除了预定面板/候选重叠，匹配对照基因仍可能重叠，不能视作完全独立测量。</p>
<h2>可复核文件</h2>'''
for path in ['provenance/prespecified_analysis.json','provenance/independent_validation.json','03_receiver/IKZF2_STAT1_primary_and_sensitivity.tsv','03_receiver/candidate_disease_DE.tsv','04_chains/role_specific_candidate_matrix.tsv','04_chains/old_TF_chain_cards.tsv','02_sender/source_specific_activity_ranks.tsv','02_sender/lead_ligand_source_identity_audit.tsv','06_support/Core_fixed_composition_focal.tsv','05_human/human_program_group_tests.tsv','provenance/BH_families.tsv']:
 body+='<p>'+link(path)+'</p>'
body+='<p>先前阶段目录：../customer_sequence_20261015_v2；../candidate_novelty_direction_20261015_v3。旧结果未覆盖。完整输出、脚本、来源和退出记录随本目录保存。</p>'
(O/'REPORT_INTERNAL.html').write_text('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>矽肺整合分析 v4</title><style>body{font:16px/1.7 "Microsoft YaHei",sans-serif;max-width:1200px;margin:35px auto;padding:0 24px;color:#203047}h1,h2{color:#163c60}h2{margin-top:36px}table{border-collapse:collapse;font-size:12px;display:block;overflow:auto}td,th{padding:5px 9px;border:1px solid #ccd6df}img{width:100%;max-width:1100px}a{color:#086b99}</style>'+body+'</html>',encoding='utf-8')
(O/'README.txt').write_text('v4 内部分析记录，2026-10-15。先打开 REPORT_INTERNAL.html。所有结论以样本级统计表为准；并未确立完整巨噬细胞→TF→焦亡机制。状态内拟时序不作必须条件。旧结果保留。\n',encoding='utf-8')
print('REPORT_AND_FIGURES_COMPLETE')