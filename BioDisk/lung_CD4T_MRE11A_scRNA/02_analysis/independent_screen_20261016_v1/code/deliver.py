from pathlib import Path
import pandas as pd,html,json
O=Path('C:/codex_test/output/silicosis_independent_v1/independent_screen_20261016_v1')
read=lambda p:pd.read_csv(O/p,sep='\t')
h=read('intersection/nominal_intersection.tsv');c=read('intersection/screen_counts.tsv');s=read('scores/subtype_rankings.tsv');m=read('scores/module_members.tsv');t=read('trajectory/all_gene_and_program_tests.tsv.gz')
h['program_membership']=h.feature.map(m.groupby('gene').module.apply(lambda v:';'.join(sorted(set(v))))).fillna('')
h.to_csv(O/'intersection/nominal_intersection.tsv',sep='\t',index=False,encoding='utf-8')
print('UNIQUE',h.feature.nunique(),'subtype_unique',h[h.scope.ne('Core')].feature.nunique())
print('AUTOPHAGY MEMBERS',h[h.program_membership.str.contains('Autophagy|Mitophagy')][['scope','feature','program_membership']].to_string(index=False))
print('STABLE ACT',','.join(h[(h.scope=='ActCD4')&(h.nominal_configurations==4)].feature))
print('STABLE TREG',','.join(h[(h.scope=='Treg')&(h.nominal_configurations==4)].feature))
print(read('trajectory/common_support.tsv').to_string(index=False));print(read('time/sample_inclusion.tsv').to_string(index=False))
parts=['''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>矽肺 CD4 独立筛选</title><style>body{font-family:Microsoft YaHei,Arial,sans-serif;max-width:1200px;margin:40px auto;line-height:1.65;color:#193148}table{border-collapse:collapse;font-size:12px;width:100%;margin:20px 0}th,td{border:1px solid #cad5df;padding:6px;text-align:left}th{background:#edf3f7}h1,h2{color:#125b70}p{max-width:1050px}</style><h1>矽肺 CD4 独立筛选结果</h1><p>2026-10-16 · 独立探索版 v1 · 未使用既有候选名单或疾病差异表达作为准入条件。</p>''']
def table(d):return d.to_html(index=False,border=0,float_format=lambda x:f'{x:.5g}',escape=True)
parts+=['<h2>亚群筛选产生 18、12、9 个交集候选</h2><p>ActCD4、Naive、Treg 分别获得 18、12、9 个候选，四种图参数／起点配置均保留的分别为 16、0、7 个。Core 是全部 CD4 的合并分析，与亚群结果有重叠，不应相加。候选均满足主分析拟时序原始 P&lt;0.05 和真实时间总体原始 P&lt;0.05；没有候选同时通过两侧 FDR&lt;0.05。</p>',table(c)]
parts+=['<h2>焦亡评分最高的亚群为 ActCD4，凋亡评分最高的为 Treg</h2><p>以下是硅尘组样本等权均值的亚群内排序。不同程序的分值不用于相互比较。Treg 的凋亡分值虽最高，相对对照却下降；不能将“最高”解释为“诱导最强”。所有疾病组比较的双侧精确置换 P 均不低于 0.10。</p>',table(s[['module','scope','case_mean','difference','case_score_rank','increase_rank','P','BH']])]
parts+=['<h2>自噬轨迹未达到名义显著，不能把交集直接解释为自噬驱动</h2><p>四个分析范围的广义自噬、线粒体自噬轨迹均未达到原始 P&lt;0.05。ActCD4 的 MitoISR 轨迹为名义显著（P=0.03125，FDR=0.25），但不是线粒体损伤发生或自噬通量改变的直接测量。以下完整保留八种程序轨迹检验。</p>',table(t[(t.run=='K15_R1')&(t.type=='program')][['scope','feature','endpoint_difference','shape','P','BH']])]
parts+=['<h2>全部交集候选：按参数稳定性及轨迹幅度排列</h2><p>排序未使用文献新颖性或既有目标。curve_range 是共同支持区间内拟合平均轨迹的 log1p(CP10k) 幅度；不是疾病组差异效应。nominal_configurations 表示四种配置中达到轨迹 P&lt;0.05 的次数，不等于方向一致性。D28_D7_logFC 是真实时间差异；时间与状态轴方向不同不作为剔除条件。program_membership 仅表示预定义基因集成员关系，不证明调控作用。</p>',table(h[['scope','feature','curve_range','shape','P','BH','time_P','time_FDR','D28_D7_logFC','nominal_configurations','program_membership']])]
parts+=['''<h2>方法与可复核口径</h2><p>GSE314943：使用固定 QC 后的 2,311 个 CD4 细胞和既有 Naive／ActCD4／Treg 标签，6 个独立生物学样本（对照、硅尘各 3）。原始计数按原始细胞文库量归一化为 CP10k 后 log1p。八种程序用基因均值减匹配表达背景均值重新计算，背景不排除任何预先指定的候选基因；每个程序基因选择最多 20 个最近表达背景基因，控制集合取并集。组间分析先按生物学样本平均，遍历全部 20 种 3 对 3 标签分配作双侧检验。</p><p>用原有 30 维 Harmony 坐标新建 15／30 邻居图。每张图使用包含全部 6 个样本及 Naive 的最大连通分量。根从 TCF7／LEF1／CCR7 均值最高四分位的 Naive 细胞中选取最接近该集合中心者，第二根来自不同样本。四种配置均进行 DPT；K15_R1 为预定主分析。每个分析范围仅使用六个样本共有的拟时序范围；基因表达门槛为至少 3 个样本伪 bulk CPM≥1。逐样本拟合 log1p(CP10k) 对拟时序线性项、二次项及中心化 log1p 文库量的回归。对六个样本的两项轨迹系数作全部 64 种联合符号翻转，检验平均曲线变化，BH 在每个范围／配置的完整基因家族内计算；程序另行校正。</p><p>GSE274372：12 个 TLike 样本，D0、D7、D28、D56 各 3 个。使用原始伪 bulk 计数、edgeR filterByExpr、保留原始文库量、TMM、robust 离散度及 quasi-likelihood 模型，以分类时间建模。总体检验为三个相对 D0 对比的联合 F 检验；另作 D7−D0、D28−D0、D56−D0、D28−D7、D56−D28。各检验在完整基因家族中分别 BH。共有 732 个基因通过时间数据表达门槛，其中 55 个总体原始 P&lt;0.05。规范化基因符号用于交集，保留原始 ID。</p><h2>解释边界</h2><p>这轮结果定位为探索性筛选。六样本轨迹检验的双侧最小 P 为 0.03125，疾病组精确比较最小 P 为 0.10；细胞数不作为统计样本量。符号翻转依赖独立样本轨迹系数在零假设下的对称性。拟时序反映转录状态排序，不能等同真实病程，也不能证明 Naive 经某个分支分化为 Treg。共同区间、根与邻居数敏感性已单独保存，Naive 候选稳定性较弱。</p><p>真实时间数据为 TLike，不是已证实纯 CD4 或对应亚群；其中 D0 一个文库仅 8 个细胞，表达门槛后仅 732 个基因可检验，时间侧覆盖是交集的重要限制。后续时间点缺少同期 saline 对照，时间变化不能完全分离年龄与暴露时长。程序评分是转录水平代理指标，不能直接确认焦亡、凋亡、线粒体损伤或自噬通量。此次保留所有名义显著候选，不因 FDR 未通过全部排除，也不将其称为校正后显著。</p></html>''']
out=O.parent/'silicosis_independent_screen_v1.html';out.write_text('\n'.join(parts),encoding='utf-8')
assert len(h)==82 and h[['scope','feature']].duplicated().sum()==0
assert (h.P<.05).all() and (h.time_P<.05).all()
assert read('scores/disease_score_tests.tsv').P.min()>=.1
print('VERIFIED',out)
