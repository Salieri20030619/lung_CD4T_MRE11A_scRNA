from pathlib import Path
import sys,json,hashlib,html
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy.stats import false_discovery_control
O=Path(sys.argv[1]);V=O.parent/'integrated_macrophage_CD4_20261015_v4';P=O.parents[1];rd=lambda n:pd.read_csv(O/'tables'/n,sep='\t');a=rd('STAT1_gene_associations.tsv');e=rd('gene_sample_expression.tsv');act=pd.read_csv(V/'03_receiver/TF_all_regulon_sample_activity.tsv.gz',sep='\t');S=['Saline_1','Saline_2','Saline_3','Silica_1','Silica_2','Silica_3']
fig,axes=plt.subplots(1,3,figsize=(13,4),layout='constrained')
for ax,scope in zip(axes,['Core','ActCD4','Treg']):
 x=act[(act.scope==scope)&(act.seed==17)&(act.regulon=='Stat1(+)')].set_index('sample').loc[S,'mean_auc'];y=e[(e.scope==scope)&(e.gene=='GSDMD')].set_index('sample').loc[S,'CP10k'];r=a[(a.scope==scope)&(a.seed==17)&(a.model=='Pooled')&(a.y=='GSDMD')].iloc[0]
 for i,s in enumerate(S):ax.scatter(x.iloc[i],y.iloc[i],color='#cd6038' if i>=3 else '#296f9b',s=65);ax.annotate(s,(x.iloc[i],y.iloc[i]),fontsize=7,xytext=(3,-13 if scope=='Core' and s=='Silica_1' else 5),textcoords='offset points')
 ax.xaxis.set_major_locator(MaxNLocator(4));ax.margins(.22);ax.set_title(f'{scope}: rho={r.rho:.3f}\nP={r.P:.4f}; BH={r.BH:.4f}',fontsize=10);ax.set_xlabel('STAT1 frozen activity (seed 17)');ax.set_ylabel('GSDMD mean CP10k')
fig.suptitle('Mouse STAT1-GSDMD transcript association; no protein cleavage measurement');fig.savefig(O/'figures/STAT1_GSDMD_samples.png',dpi=170);plt.close(fig)
h=rd('human_gene_effects.tsv');fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
for ax,g in zip(axes,['STAT1','GSDMD']):
 b=h[h.gene==g].set_index('scope').loc[['Core','MemCD4','NaiveCD4','Treg']];ax.bar(b.index,b.difference,color=['#296f9b' if z<0 else '#cd6038' for z in b.difference]);ax.axhline(0,color='black',lw=.7);ax.set_title(g);ax.set_ylabel('Silicosis - exposed mean CP10k');ax.tick_params(axis='x',rotation=15)
fig.suptitle('Human reviewed CD4 expression: no nominal replication');fig.savefig(O/'figures/human_direction.png',dpi=170);plt.close(fig)
# Output/source integrity and all newly added BH families.
bhchecks=0
for name in ['STAT1_gene_associations.tsv','fixed_composition_STAT1_gene.tsv','human_gene_effects.tsv']:
 d=rd(name)
 for _,b in d.groupby('family'):
  ok=b.P.notna();assert np.allclose(false_discovery_control(b.loc[ok,'P']),b.loc[ok,'BH']);bhchecks+=1
sources=rd('input_manifest.tsv').to_dict('records')+json.loads((O/'provenance/human_sources.json').read_text())
for r in sources:assert hashlib.sha256(Path(r['path']).read_bytes()).hexdigest()==r['sha256']
qc=rd('CD4_executor_identity_audit.tsv');roster=qc.drop_duplicates(['scope','sample'])[['scope','sample','cells','kept_cells']];roster.to_csv(O/'tables/strict_identity_sensitivity_roster.tsv',sep='\t',index=False)
# Read-only inventory: identify available ambient-correction inputs, no new cleanup run.
paths=[]
for root in [P/'00_raw/GSE314943',P/'02_analysis/GSE314943']:
 if root.exists():
  for f in root.rglob('*'):
   if f.is_file() and any(s in f.name.lower() for s in ['raw_feature','filtered_feature','ambient','soupx','decont','emptydrop']):paths.append(str(f))
(O/'provenance/ambient_input_inventory.json').write_text(json.dumps({'searched_roots':[str(P/'00_raw/GSE314943'),str(P/'02_analysis/GSE314943')],'matching_files':paths,'interpretation':'filename inventory only; not evidence that every contamination-control input is absent or present'},indent=2),encoding='utf-8')
def tbl(d):return d.to_html(index=False,float_format=lambda x:f'{x:.5g}',na_rep='NA')
b='''<h1>STAT1—焦亡链条补充分析 v5</h1><p>2026-10-15 · 内部工作报告 · 新目录保存，旧版本未修改。</p>
<h2>新增证据支持优先验证 STAT1–GSDMD；尚不能接通巨噬细胞来源的完整机制</h2><p>将旧 PyroExec 分数拆开后，发现 STAT1 活性与 GSDMD 转录表达在 Core、ActCD4、Treg 中正向关联。原 PyroExec 五基因总评分仍无名义显著关联，没有替换评分或修改旧结果。新结果是明确标注的事后探索，不作为独立外部验证。</p><p><strong>更关键的新限制是细胞来源可靠性：</strong>当前冻结 CD4 矩阵有广泛肺上皮标记共表达。严格排除异源标记阳性细胞后 Saline_2 全部排空，六样本敏感性无法评估。不能据此断定全部为双细胞，也不能宣称已排除环境 RNA 对 GSDMD 的贡献。应先核实背景 RNA，再开展机制主张。</p>
<h2>执行基因分解补出了 GSDMD 分支</h2><p>旧 PyroExec 包含 PYCARD、CASP1、GSDMD、IL1B、IL18；其中有执行成分，也有炎症因子。GSDME 不在旧面板中。本次固定检查21个名称，鼠 Casp4 是历史 caspase-11 的官方基因名；CASP11 不另作独立基因，鼠可测的独立特征为20个。GSDMD 等旧面板基因已从冻结 TF regulon 中剔除，关联不是直接将 GSDMD 同时计入活性和终点。新加的 IRF/GBP/ZBP1 等未保证从旧 regulon 中剔除，其关联标记为非独立探索。</p>'''
key=a[(a.seed==17)&a.y.isin(['GSDMD','CASP1'])][['scope','model','y','rho','P','BH','loo_same_sign']];b+=tbl(key)
b+='''<img src="figures/STAT1_GSDMD_samples.png"><p>固定亚群组成后 Core STAT1–GSDMD：rho=1，P=0.002778，BH=0.02778，达到这次定向基因家族的 BH&lt;0.05；疾病组别调整后 P=0.05556、BH=0.5556。这个校正后显著属于组成敏感性分析的限定检验家族，不能写成原始主分析全面通过 FDR，也不能把控制疾病组别后的结果写成独立关联已确立。</p><p>三个随机种子中 GSDMD 的关联方向保持一致；种子只是算法敏感性，不增加生物学样本量。Core CASP1 与 STAT1 的关联在固定组成后不再名义显著。GSDME 主种子未达名义显著，不能换种子43来确立另一分支。IRF1 未出现一致上调或与 STAT1 的关联，不将文献中的 STAT1→IRF1 直接套到本数据。</p>'''
b+=tbl(rd('gene_disease_DE_original.tsv').loc[lambda z:z.canonical.isin(['STAT1','IRF1','CASP1','GSDMD','GSDME']),['scope','canonical','logFC','PValue','FDR']])
b+='''<h2>上游有探索性替代线索，但尚无合格的特异性配体</h2><p>TNF、IL15、IL18→STAT1 活性原分析仍未达到配对关联筛选标准。IL27、IFNG 也未形成充分统计支持。小鼠 IFNAB 对应原始 Ifnab，不能当作所有 I 型干扰素的合并读数；其部分负相关不支持“上调干扰素驱动 STAT1”的简单模型。</p><p>C3 与 CXCL16 可作为来源侧待核实线索：AMlike C3–Core STAT1 rho=1、P=0.002778、BH=0.4353；AllMac CXCL16–ActCD4 STAT1 rho=1、P=0.002778、BH=0.4811。它们不能因相关就被称作直接激活 STAT1 的配体，需核实兼容受体、信号机制和疾病混杂。所有上游 BH 保留原完整配体检验家族，未对筛出的几项重新缩小校正池。</p><p>病例 AMlike 中 IFNG 阳性细胞合计1个，IFNB1/IFNA1均为0；AllMac 中分别2、3、0个。IL12A 虽在 AllMac–ActCD4/Treg 有名义正相关，但病例 AMlike 中 IL12A+IL12B 同细胞共检出为0，AllMac 为1；IL27+EBI3 对应为0和5。单亚基、掉落与低表达限制必须分别解释，未检出不等于证明没有蛋白。</p><p>程序输出279条“上游关联+执行基因关联”组合，含正负方向及重叠状态，只是全量事后枚举；不是279个验证通路，不生成链条 P 值，也不以任意两段名义显著作为因果成立。未在6只小鼠上拟合中介效应模型。</p>
<h2>身份审计暴露了需优先解决的环境 RNA 问题</h2><p>采用固定诊断规则：至少2个髓系标记或至少2个上皮标记检出，则标记为需检查；另记录至少2个 T 细胞标记支持。规则用于压力测试，不用于重新注释。Saline_2 的424个 CD4 细胞均检出 SFTPC 和 SCGB1A1，其中不少也有 CD3/TRAC 支持，因此不能简单删除为非 T 细胞。该样本 GSDMD 阳性31个细胞，每个仅1 UMI；存在环境来源解释空间。</p>'''
b+=tbl(qc[(qc.scope=='Core')&(qc.gene=='GSDMD')]);b+='<p>严格排除后的实际样本覆盖：</p>'+tbl(roster)
b+='''<p>全部四个状态均因 Saline_2 为0而无法保留六个配对样本。未把缺失补零，未悄悄改用5样本，未放松规则以追求结果。单细胞转录质控不能替代分选 CD4 的蛋白验证；本次不覆盖冻结注释或原矩阵。下一轮计算应从原始/未过滤液滴或已校正对象识别环境背景，比较校正前后 GSDMD、CASP1、STAT1 与原始疾病差异的稳定性。</p>
<h2>患者数据尚不支持同方向的 STAT1–GSDMD 升高</h2><p>GSE174725 复核注释下，Core GSDMD 病例均值比暴露对照低0.5902 CP10k，P=0.30、BH=1；STAT1 同样较低。MemCD4 的 GSDMD 方向向上但仅2对照+2病例，P=0.3333。不能用 MemCD4 替代整体结果，也不能把跨物种阴性静默删除。人类患者 AM 覆盖缺口仍在，没有新增有效的人类配体—STAT1 配对验证。</p><img src="figures/human_direction.png">'''
b+=tbl(h[h.gene.isin(['STAT1','GSDMD','CASP1','GSDME'])])
b+='''<h2>真正接通机制需补齐的证据与顺序</h2><table><tr><th>优先级</th><th>缺口</th><th>补充与通过条件</th></tr><tr><td>1</td><td>信号确实位于 CD4</td><td>先核查/校正环境 RNA；独立分选或定位 CD3+CD4+ 亚群，检测 STAT1、GSDMD，并排除巨噬细胞混入。校正后方向稳定才推进主要计算机制解释。</td></tr><tr><td>2</td><td>转录状态对应真实焦亡</td><td>同步测 p-STAT1、GSDMD 裂解/膜定位、炎性 caspase 活化及膜通透性/裂解性死亡。LDH、IL1B 或 RNA 单独阳性不充分；应证明 CD4 内的 GSDMD 依赖性。线粒体 ROS/膜电位独立测定，不能用“线粒体损伤（焦亡）”视为同义。</td></tr><tr><td>3</td><td>巨噬细胞介导作用</td><td>同来源巨噬细胞条件培养或共培养，匹配对照、颗粒残留/直接硅尘作用对照；先确认能诱导 CD4 表型，再测上清蛋白和接收端 p-STAT1 的时间顺序。</td></tr><tr><td>4</td><td>锁定配体与受体</td><td>依据上清实测及受体选择候选，用中和/受体干预和适当补充/救援验证。TNF/IL15 保留旧候选；IFN/IL27 检查但不预设来自巨噬细胞；C3/CXCL16 为低等级替代线索。</td></tr><tr><td>5</td><td>STAT1 与死亡执行的因果次序</td><td>接收端 STAT1 干预应改变 GSDMD 活化及死亡，GSDMD 干预应影响死亡但不必影响早期 STAT1 活化。配合适当救援和重复，区分一般存活/激活变化。具体方案按实验体系确定。</td></tr></table>
<p>状态内拟时序不是这些步骤的必要条件，也不能替代上述干预证据。现阶段实验主线宜写作“巨噬细胞相关信号与 CD4 STAT1–GSDMD 相关转录状态及焦亡易感性的关系”，仍须先验证细胞身份和执行表型。</p>
<h2>相关文献与创新性边界</h2><p>矽肺中的 GSDMD/GSDME 机制已有以巨噬细胞为重点的研究，不能把“矽肺+焦亡”本身当成新发现：<a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC9069405/">GSDMD 与硅尘炎症/纤维化</a>；<a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC9718385/">Caspase-1/Gsdmd 与 Caspase-3/-8/Gsdme 分支</a>。STAT1 的矽肺蛋白修饰亦已有研究：<a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC11107463/">2024 肺组织蛋白组研究</a>。可探索的新问题是特定 CD4 状态中的执行事件及巨噬细胞到该状态的功能联系，不能据本轮检索宣称首次。</p><p><a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC7674074/">TNF 与 IFNG 联合驱动 STAT1/IRF1 炎症性死亡研究</a>提供的是其他细胞/实验背景的先验，不是本项目 CD4 的实测机制。<a href="https://www.ncbi.nlm.nih.gov/gene/12363/">NCBI Casp4/Casp11 名称说明</a>用于避免将同一鼠基因重复计算。</p>
<h2>统计与文件</h2><p>独立小鼠 n=6；精确样本级置换，池化720种，疾病组内36种。全部非空样本保留，无至少20细胞门槛。相关性家族按状态/种子/模型内的全部可测目标基因定义；固定组成为单独完整家族；人类每状态为完整目标基因家族。随机种子不是样本。所有名义、BH、方向、不可评估结果和来源表完整保存。</p>'''
for n in ['STAT1_gene_associations.tsv','gene_disease_DE_original.tsv','fixed_composition_STAT1_gene.tsv','mechanistic_ligand_STAT1.tsv','all_nominal_upstream_alternatives.tsv','alternative_receptor_coverage.tsv','source_complex_identity_audit.tsv','CD4_executor_identity_audit.tsv','identity_marker_counts.tsv','human_gene_effects.tsv','signed_prior_edges.tsv','posthoc_two_segment_candidates.tsv']:
 b+=f'<p><a href="tables/{n}">{n}</a></p>'
(O/'REPORT_INTERNAL.html').write_text('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>STAT1焦亡补充分析v5</title><style>body{font:16px/1.75 "Microsoft YaHei",sans-serif;max-width:1180px;margin:36px auto;padding:0 22px;color:#203247}h1,h2{color:#193f62}h2{margin-top:32px}table{border-collapse:collapse;display:block;overflow:auto;font-size:12px}td,th{padding:6px 9px;border:1px solid #ccd6df}img{width:100%}a{color:#086b99}</style>'+b+'</html>',encoding='utf-8')
(O/'README.txt').write_text('v5内部分析。打开REPORT_INTERNAL.html。新发现：STAT1-GSDMD转录关联；关键未决：来源背景RNA、配体特异性、患者方向和蛋白执行。不是已证实焦亡机制。原结果未修改。\n',encoding='utf-8')
(O/'provenance/final_integrity.json').write_text(json.dumps({'status':'PASS','source_hash_checks':len(sources),'new_BH_families_verified':bhchecks,'strict_identity_sensitivity':'not evaluable with six mice because Saline_2 has zero retained cells','causal_chain':'not established'},indent=2),encoding='utf-8');print('REPORT_COMPLETE',len(sources),bhchecks,'AMBIENT_INPUT_MATCHES',len(paths));print('\n'.join(paths[:12]))