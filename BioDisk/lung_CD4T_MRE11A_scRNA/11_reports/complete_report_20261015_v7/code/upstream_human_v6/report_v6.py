from pathlib import Path
import sys,json,hashlib,itertools
import numpy as np,pandas as pd
from scipy.stats import false_discovery_control
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(sys.argv[1]);r=lambda n:pd.read_csv(O/'tables'/n,sep='\t');stats=r('human_donor_group_effects.tsv');expr=r('human_donor_expression.tsv');checks=0
for z in stats.itertuples():
 b=expr[(expr.variant==z.variant)&(expr.scope==z.scope)&(expr.gene==z.gene)];c=b[b.group=='Exposed'].CP10k.dropna().to_numpy();s=b[b.group=='Silicosis'].CP10k.dropna().to_numpy()
 if np.isfinite(z.P):
  v=np.r_[c,s];vals=[]
  for ix in itertools.combinations(range(len(v)),len(c)):
   m=np.zeros(len(v),bool);m[list(ix)]=True;vals.append(abs(v[~m].mean()-v[m].mean()))
  p=np.mean(np.array(vals)>=abs(s.mean()-c.mean())-1e-12);assert abs(p-z.P)<1e-12;checks+=1
for _,b in stats.groupby('family'):
 ok=b.P.notna()
 if ok.any():assert np.allclose(false_discovery_control(b.loc[ok,'P']),b.loc[ok,'BH'])
for x in json.loads((O/'provenance/input_hashes.json').read_text()):assert hashlib.sha256(Path(x['path']).read_bytes()).hexdigest()==x['sha256']
pat=r('patient2_localization.tsv');joint=r('patient2_complex_codetection.tsv');ident=r('patient2_CD4_identity_codetection.tsv');dec=r('human_Core_descriptive_decomposition.tsv');loo=r('human_leave_one_donor_effects.tsv');up=r('mouse_upstream_integrated.tsv');ids=['B2_RZ725','B6_RZ725']
fig,axes=plt.subplots(1,3,figsize=(12,4),layout='constrained')
for ax,(scope,g) in zip(axes,[('Alveolar macrophage','CXCL16'),('CD4 effector T cell','CXCR6'),('CD4 effector T cell','GSDMD')]):
 b=pat[(pat.scope==scope)&(pat.gene==g)].set_index('sample').loc[ids];ax.bar(ids,b.pct,color=['#276d96','#cc653d']);ax.set_ylim(0,max(25,b.pct.max()*1.25));ax.set_ylabel('Detected cells (%)');ax.set_title(scope+'\n'+g,fontsize=10)
 for i,(_,z) in enumerate(b.iterrows()):ax.text(i,z.pct+.7,f'{int(z.positive_cells)}/{int(z.cells)}',ha='center',fontsize=10)
fig.suptitle('Two silicosis patients: localization only, no disease-control inference');fig.savefig(O/'figures/patient_localization.png',dpi=170);plt.close(fig)
fig,ax=plt.subplots(figsize=(8,4),layout='constrained')
for i,g in enumerate(['STAT1','GSDMD']):
 b=dec[dec.gene==g];ax.bar(i-.16,b.composition_component.sum(),width=.3,label='Composition' if i==0 else None,color='#6795b2');ax.bar(i+.16,b.within_state_component.sum(),width=.3,label='Within-state component' if i==0 else None,color='#ce7552')
ax.set_xticks([0,1],['STAT1','GSDMD']);ax.axhline(0,color='black',lw=.7);ax.set_ylabel('Contribution to case-control mean CP10k difference');ax.set_title('Human Core lower expression is mainly within-state');ax.legend();fig.savefig(O/'figures/human_composition_decomposition.png',dpi=170);plt.close(fig)
def table(d):return d.to_html(index=False,float_format=lambda x:f'{x:.5g}',na_rep='NA')
b='''<h1>上游配体与矽肺患者支持：v6 补充结果</h1><p>2026-10-15 · 内部分析 · 新目录保存，旧结果不覆盖。</p><h2>CXCL16–CXCR6 的来源与患者定位得到补强，直接连接 STAT1 的证据仍缺失</h2><p>本轮优先保留 CXCL16–CXCR6 为巨噬细胞与 CD4 募集/滞留相关候选；STAT1–GSDMD 作为接收端探索线索继续独立评价。二者之间画待验证连接，不写成已成立的连续机制。C3 的小鼠差异和相关性更强，但特异受体支持弱，不宜仅按 P 值将其升为首选。</p><h2>小鼠来源差异与配对关联不能互相替代</h2><p>AMlike CXCL16：logFC=1.8758，P=0.002370，BH=0.1065，疾病升高达到名义显著；但 AMlike CXCL16 与 ActCD4 STAT1 活性 rho=0.02857、P=1。AllMac CXCL16 与该活性 rho=1、P=0.002778、BH=0.4811，但 AllMac 的 CXCL16 疾病差异 P=0.1487。来源定义改变了结论，不能将 AMlike 的差异与 AllMac 的相关性拼成一条已验证链。</p><p>病例 ActCD4 的 CXCR6 检出率为28.17%、30.00%、53.78%，提供接收侧覆盖。C3 在 AllMac 中上调：logFC=1.2862，P=0.000264，BH=0.01129；然而小鼠 ActCD4/Treg 的 C3AR1 仅1/3病例有检出，另两例为0。C3 的广义先验受体名单不能替代具体 C3a–C3AR1 支持；也未测 C3 蛋白裂解产物。</p><p>IL27 虽有 CD4–STAT1 的机制文献先验，但来源亚基和完整受体覆盖不稳。IL12A 在 AllMac 的疾病差异通过 BH，不等于完整 IL-12 蛋白存在。TNF/IL15 保留旧候选，当前不能升级为直接驱动 STAT1 的证据。</p>'''
b+=table(up[(up.model=='Pooled')&(up.scope=='ActCD4')][['sender','x','rho','P','BH','ligand_logFC','ligand_P','ligand_BH']])
b+='''<h2>两位患者提供了新增定位支持</h2><p>从与旧核查一致的官方 GSE296682 版本提取31个配体、受体和执行相关基因；复用官方范围缓存，本轮新增网络块请求为0。全部4579个患者细胞按作者标签保留，没有最低细胞数门槛。</p><p>两位患者肺泡巨噬细胞 CXCL16 阳性分别264/323（81.73%）和601/782（76.85%）；CD4 effector T cell 的 CXCR6 阳性分别5/33（15.15%）和2/27（7.41%）。这补强了配体来源与受体定位，同时不证明分泌、物理接触、疾病上调或 STAT1 激活。</p><img src="figures/patient_localization.png"><p>有限身份核查中，两位患者各有1个 CXCR6 阳性 CD4 effector 细胞同时检出至少2个 T 细胞标记；CXCR6 阳性细胞没有达到预设的至少2个髓系标记共表达标记。其余阳性细胞不能因 T 标记掉落直接删除，但需进一步定位确认。不能把作者 CD4 effector 标签等同于小鼠 ActCD4。</p><p>GSDMD 在两位患者 CD4 effector 中分别检出4/33、2/27；STAT1+GSDMD 同细胞共检出为3/33和0/27，前者3个均带髓系混合标记。CXCR6+STAT1+GSDMD 三重共检出均为0。因此，新增数据支持候选分子的定位，尚不能补强完整受体→STAT1→GSDMD 的同细胞证据。</p><p>C3AR1 在两位患者所有作者 CD4 effector、helper、Regulatory T 标签中均未检出；IL27 完整受体 IL27RA+IL6ST 在 CD4 effector 为0/33和4/27，供者覆盖不足。上述未检出受单细胞检测灵敏度约束，不作绝对不存在的结论。</p>'''
b+=table(ident)
b+='''<h2>患者疾病方向不能靠组成或单个供者解释掉</h2><p>GSE174725 保持原复核注释与独立供者统计：Core STAT1 差值=-0.9190 CP10k，P=0.10；GSDMD=-0.5902，P=0.30。本轮31基因家族 BH 分别0.8667和1，与旧较小家族 BH 不混用。排除既有 doublet/混合标记后方向及效应近似不变。</p><p>逐个剔除任一供者，Core STAT1 差值始终为负（-1.1431至-0.6949），GSDMD 也始终为负（-1.0287至-0.1517）。这些留一结果仅报告描述性效应，只有1个对照的情况不作推断检验。</p><p>等供者 Core 均值作对称组成分解：STAT1 总差值-0.9190中，组成部分约-0.0476，状态内部分约-0.8714；GSDMD 总差值-0.5902中，组成部分约-0.0862，状态内部分约-0.5040。降低主要不来自亚群比例变化。分解属于描述性恒等分解，不是因果分解，也没有为缺失状态填入虚构表达值。</p><img src="figures/human_composition_decomposition.png"><p>MemCD4 中 CXCR6 和 GSDMD 的差值分别+0.2352和+0.7211 CP10k，原始 P 均0.3333；STAT1 差值+0.0234，P=1。方向可保留为特定状态线索，但该状态仅2对照+2病例有细胞，不能用它替代 Core 的关键反向结果。</p>'''
b+=table(stats[(stats.variant=='Reviewed')&stats.scope.isin(['Core','MemCD4','NaiveCD4','Treg'])&stats.gene.isin(['STAT1','GSDMD','CXCR6','IL27RA'])])
b+='''<h2>本轮后的候选优先级</h2><table><tr><th>候选</th><th>可保留的作用假设</th><th>尚缺证据</th></tr><tr><td>CXCL16–CXCR6</td><td>优先的募集/滞留相关上游候选；小鼠来源差异及两位患者定位</td><td>同一来源的稳定配对关联、STAT1连接、功能干预</td></tr><tr><td>C3</td><td>小鼠巨噬细胞疾病响应/通讯探索线索</td><td>具体活化产物和CD4特异受体；患者C3AR1未检出</td></tr><tr><td>IL27</td><td>有直接CD4–STAT1文献先验的机制对照候选</td><td>来源完整亚基、跨供者完整受体及项目内统计关联</td></tr><tr><td>TNF/IL15</td><td>保留原有候选供比较</td><td>当前仍无证据将其确立为STAT1直接上游</td></tr></table><p>实验上可先验证巨噬细胞 CXCL16 蛋白和 CD4 表面 CXCR6，再区分迁移/滞留效应与细胞内 p-STAT1、GSDMD 活化效应。若只影响迁移，不应强行表述为直接诱导焦亡。v5 已发现的环境 RNA 问题尚未解决，相关小鼠基因关联继续作为探索性结果。</p><h2>文献与统计边界</h2><p>CXCL16/CXCR6 与肺部 T 细胞定位有既有研究：<a href="https://rupress.org/jem/article/216/12/2748/132521/CXCR6-regulates-localization-of-tissue-resident">肺部 CXCR6 定位研究（CD8背景）</a>；<a href="https://www.nature.com/articles/s42003-024-06662-z">硅尘 CD4 组织驻留 T 细胞研究</a>。因此，本轮不把这条轴本身称为全新发现。<a href="https://pubmed.ncbi.nlm.nih.gov/15356135/">IL27–STAT1 的 CD4 实验研究</a>提供先验机制，不代表本项目已验证该通路。</p><p>旧 Human SCENIC 活性文件来自煤工尘肺，已明确排除。本轮没有给 GSE174725 伪造 STAT1 活性，只分析其 RNA、表达覆盖及既有注释。GSE296682 仅患者定位，没有对照、不计算疾病差异或两供者相关 P。GSE174725 的 AM 病例有效覆盖仍仅1例，不计算巨噬细胞组间检验；单核细胞对照保持独立标签。</p><p>以供者为统计单位，缺失保持 NA，全零不可检验，不设至少20细胞门槛；原始P筛选与BH证据分开报告。小鼠上游BH保持原完整家族，人类本轮以每版本/状态的31基因为同类家族。新增所有可检验P已独立枚举验证、BH回算、组成分解恒等关系和数据来源哈希均通过。</p>'''
for n in ['mouse_upstream_integrated.tsv','mouse_receptor_specificity.tsv','human_donor_group_effects.tsv','human_leave_one_donor_effects.tsv','human_Core_descriptive_decomposition.tsv','human_complex_codetection.tsv','patient2_localization.tsv','patient2_complex_codetection.tsv','patient2_CD4_identity_codetection.tsv']:
 b+=f'<p><a href="tables/{n}">{n}</a></p>'
(O/'REPORT_INTERNAL.html').write_text('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>上游与患者支持v6</title><style>body{font:16px/1.75 "Microsoft YaHei",sans-serif;max-width:1180px;margin:35px auto;padding:0 24px;color:#21354a}h1,h2{color:#174363}h2{margin-top:32px}table{border-collapse:collapse;display:block;overflow:auto;font-size:12px}td,th{border:1px solid #ccd5df;padding:6px 8px}img{width:100%}a{color:#086b99}</style>'+b+'</html>',encoding='utf-8')
(O/'README.txt').write_text('v6内部分析；打开REPORT_INTERNAL.html。CXCL16-CXCR6来源与患者定位得到补强，但STAT1-GSDMD完整机制未确立。人类Core反向结果保留。\n',encoding='utf-8')
(O/'provenance/validation.json').write_text(json.dumps({'status':'PASS','independent_exact_P_checks':checks,'BH_families':stats.family.nunique(),'source_hash_checks':5,'composition_identity':'PASS','patient2_STAT1_parity_with_verified_official_cache':'PASS','causal_chain':'not established'},indent=2),encoding='utf-8');print('COMPLETE',checks,stats.family.nunique())