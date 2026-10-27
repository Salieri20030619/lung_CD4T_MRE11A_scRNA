from pathlib import Path
import pandas as pd,numpy as np,json,base64,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup
O=Path('C:/codex_test/output/rora_followup_audit_v1/rora_followup_audit_20261016_v1');I=Path('C:/codex_test/output/silicosis_independent_v1/independent_screen_20261016_v1');F=O.parent
read=lambda p:pd.read_csv(O/p,sep='\t')
assert (O/'exit_status').read_text().strip()=='0' and json.loads((O/'COMPLETE.json').read_text())['status']=='COMPLETE' and 'EXTRA_COMPLETE' in (O/'extra.log').read_text()
c=read('correlation_summary.tsv');w=read('within_sample_correlations.tsv');d=read('disease_comparisons.tsv');s=read('sample_expression.tsv');a=read('mechanism_alias_checked.tsv');r=read('time_RORA_raw_counts.tsv');ov=read('score_overlap_audit.tsv');h=pd.read_csv(I/'intersection/nominal_intersection.tsv',sep='\t').query("feature=='RORA'")
primary=c.query("score_mode=='target_excluded' and model=='full' and module in ['Autophagy_GO','Mitophagy']").copy();assert primary.n_samples.eq(6).all() and primary.P.ge(.05).all()
assert not ov.query("mode=='target_excluded'").removed_members.notna().any() and not ov.query("mode=='target_excluded'").removed_controls.notna().any()
print('ALIAS CHECK\n',a.query("scope=='Core'").to_string(index=False))
print('TREG ADJUST\n',d.query("scope=='Treg'").to_string(index=False))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':160})
fig,ax=plt.subplots(1,2,figsize=(11,4.1),layout='constrained');scopes=['Naive','ActCD4','Treg']
for j,scope in enumerate(scopes):
 b=s[s.scope==scope]
 for group,color,offset in [('Saline','#248d9c',-.13),('Silica','#c56531',.13)]:
  vals=b[b.group==group].RORA.to_numpy();ax[0].scatter(np.full(3,j+offset),vals,c=color,s=42,label=group if j==0 else None);ax[0].plot([j+offset-.07,j+offset+.07],[vals.mean()]*2,c=color,lw=2)
ax[0].set(xticks=range(3),xticklabels=scopes,ylabel='Sample mean RORA: log1p(CP10k)',title='A. Three independent mice per group');ax[0].legend(frameon=False)
for k,(model,color,label) in enumerate([('unadjusted','#8ca3ad','Unadjusted'),('depth_activation_tissue','#e0ae65','Depth + state proxies'),('plus_pseudotime','#247584','Also pseudotime')]):
 vals=[d[(d.scope==scope)&(d.model==model)].difference.iloc[0] for scope in scopes];ax[1].barh(np.arange(3)+(k-1)*.23,vals,height=.21,color=color,label=label)
ax[1].axvline(0,color='black',lw=.7);ax[1].set(yticks=range(3),yticklabels=scopes,xlabel='Silica minus saline: log1p(CP10k)',title='B. Adjusted difference (descriptive sensitivity)');ax[1].legend(frameon=False,fontsize=8)
fig.savefig(F/'rora_state_audit.png');plt.close(fig)
fig,ax=plt.subplots(1,2,figsize=(11,4.1),layout='constrained')
for j,mod in enumerate(['Autophagy_GO','Mitophagy']):
 b=w.query("score_mode=='target_excluded' and model=='full'");b=b[b.module==mod]
 for k,scope in enumerate(scopes):
  for group,color,offset,marker in [('Saline','#248d9c',-.11,'o'),('Silica','#c56531',.11,'^')]:
   vals=b[(b.scope==scope)&(b.group==group)].partial_r.to_numpy();ax[j].scatter(k+offset+np.linspace(-.04,.04,len(vals)),vals,c=color,s=40,marker=marker,label=group if k==0 else None)
 ax[j].axhline(0,c='gray',lw=.7);ax[j].set(xticks=range(3),xticklabels=scopes,ylabel='Within-mouse adjusted correlation',title=mod,ylim=(-.5,.5));ax[j].legend(frameon=False)
fig.savefig(F/'rora_autophagy_audit.png');plt.close(fig)
def table(df):return df.to_html(index=False,border=0,float_format=lambda v:f'{v:.5g}')
def img(n):return '<img style="width:100%" src="data:image/png;base64,'+base64.b64encode((F/n).read_bytes()).decode()+'">'
issues=[
 {'编号':'RORA-01','分类／程度':'待澄清／高','定位':'disease_comparisons.tsv：ActCD4、Treg；correlation_summary.tsv：State2','证据':'ActCD4 原始差异 +0.09185，完全调整后 −0.11231；Treg +0.84583→+0.43459。ActCD4 与组织适应代理评分相关 r=0.18538，P=0.03125，FDR=0.1875。','影响':'ActCD4 信号不能区分于一般状态变化；Treg 保留暴露关联效应，但无其他疾病对照，不能建立矽肺特异性。','建议':'保留 Treg 暴露相关候选定位；避免称为已证实的矽肺特异靶点。','复核验收':'独立样本验证 Treg 内表达／蛋白与功能；若声称特异性，应加入其他肺炎症或纤维化疾病对照。'},
 {'编号':'RORA-02','分类／程度':'证据边界／高','定位':'correlation_summary.tsv：target_excluded/full，Autophagy_GO 与 Mitophagy；score_overlap_audit.tsv','证据':'六样本平均偏相关：ActCD4 自噬 −0.09492、线粒体自噬 +0.06140；Treg +0.03368、+0.04909，均 P≥0.05。仅调整深度及去除重叠标记敏感性也未通过。RORA 未被纳入任何程序成员或背景。','影响':'当前数据不支持把 RORA 定义为自噬／线粒体自噬驱动候选；不是由 RORA 进入评分产生的循环相关。','建议':'机制表述止于待检验假设；不套用心肌 CAV3 或 COPD p53–DRAM 通路。','复核验收':'在目标 CD4 亚群干预 RORA，测量自噬／线粒体自噬通量，并用下游救援验证方向。'},
 {'编号':'RORA-03','分类／程度':'新增解释风险／高','定位':'time_RORA_raw_counts.tsv：D00_R1–R3；原交集表 RORA 的 D28_D7、D56_D28 对比','证据':'D0 三样本均 0 UMI（13、8、29 个细胞）；暴露后均检出。D28−D7 log2FC=−0.46844，P=0.24113；D56−D28 log2FC=−0.34429，P=0.40764。','影响':'总体时间 FDR=0.00964 主要支持基线与暴露后差异，不能当作 D7→D28 的显著进展变化；低细胞量基线零检出影响效应量解释。','建议':'将时间支持描述为暴露前后变化；若客户明确要求 D7／D28 对比，RORA 不通过该 P<0.05 门槛。','复核验收':'使用足量 CD4/Treg 的独立时间系列，并设置同期对照；预先明确总体时间检验或 D28−D7 对比。'}]
pd.DataFrame(issues).to_csv(O/'audit_findings.tsv',sep='\t',index=False,encoding='utf-8')
parts=['''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>RORA 两项进一步核查</title><style>body{font-family:Microsoft YaHei,Arial,sans-serif;max-width:1160px;margin:35px auto;line-height:1.7;color:#243b4d}table{border-collapse:collapse;font-size:12px;width:100%;margin:18px 0}th,td{border:1px solid #cad5df;padding:6px;text-align:left;vertical-align:top}th{background:#edf3f7}h1,h2{color:#145c70}a{color:#156283}.box{padding:15px;background:#edf3f7}</style><h1>RORA：矽肺特异性与自噬机制的进一步核查</h1><p>2026-10-16 · 独立核查 v1 · 原筛选、图件及报告保持原样。</p><p class="box"><b>结论：</b>RORA 可保留为 Treg 暴露相关／状态相关候选；当前数据未建立矽肺特异性，也不支持自噬驱动机制。总体时间显著主要来自 D0 与暴露后差异，而 D7→D28 未达到名义显著。该判断不撤销它按原“总体时间＋拟时序”口径进入候选的事实。</p>''']
parts+=['<h2>1. Treg 保留上升效应，ActCD4 缺少独立疾病差异</h2>',img('rora_state_audit.png'),'<p>每个点为一个独立生物学样本，组均值按样本等权。差异单位为平均 log1p(CP10k)，不是 log2 倍数变化。全 CD4 的原始差异为 +0.74948，控制亚群构成、测序深度、状态代理评分及拟时序后为 +0.05830；此敏感性结果不能定量解释为某个因素的因果贡献。Treg 调整后仍为正方向，三只硅尘样本的调整均值全部高于三只对照，精确 P=0.10。</p><p>3 对 3 双侧精确置换的最小 P 为 0.10；因此不能用“未过 0.05”否定 Treg 的方向和效应。另一方面，ActCD4 原始差异小且经调整变号；与组织适应代理评分为名义正相关。这些结果使一般状态变化成为需要保留的解释。</p>',table(d[['scope','model','difference','P','BH','loo_sign_agreement']])]
parts+=['<h2>2. 未检出与自噬程序一致的独立联系</h2>',img('rora_autophagy_audit.png'),'<p>逐样本控制文库深度、活化／组织适应代理评分和拟时序线性／二次项后，计算 RORA 与程序评分的偏相关。随后对六个样本 Fisher z 等权汇总并作 64 种符号翻转检验；不以细胞数构造 P 值。下表保留原始 P 和按预设家族计算的 FDR。</p>',table(primary[['scope','module','mean_Fisher_r','P','BH','positive_samples','negative_samples']]),'<p>仅控制深度、去除评分中与协变量重叠的标记基因两种敏感性结果同样未出现自噬关联 P&lt;0.05。RORA 本身既不在程序基因中，也不在背景基因中；去除 RORA 的评分与原分数一致。基因与程序相关阴性不能证明不存在任何非转录机制，但不提供现有主线所需的正向机制证据。</p>']
parts+=['<h2>3. 时间支持是基线—暴露后差异，不是已证实的晚期进展</h2>',table(r[['sample','day','cells','RORA_UMI','RORA_raw_CPM']]),table(h[h.scope=='Core'][['feature','time_P','time_FDR','D7_D0_logFC','D28_D0_logFC','D56_D0_logFC','D28_D7_logFC','D28_D7_PValue','D28_D7_FDR','D56_D28_logFC','D56_D28_PValue','D56_D28_FDR']]),'<p>D0 三份低细胞量 TLike 样本均未检出 RORA；不能把零检出写成生物学上绝对不表达。暴露后九份样本均检出，支持暴露前后差异。D28−D7 与 D56−D28 均为下降方向但未达原始 P&lt;0.05。时间侧为 TLike，无法归因于 Treg，也不能区分细胞组成、年龄及暴露时长。</p>']
parts+=['''<h2>4. 文献核查进一步限制了机制外推</h2><ul><li>2021 年肺部等组织的研究已报告 Rora 随 CD4 活化与非淋巴组织环境变化，构成本项目的一般状态解释。<a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC8130942/">原始研究</a>。</li><li>2021 年结肠炎研究通过 CD4 T 细胞 Rora 缺失和相关实验，支持 RORα 调控 mTORC1、T 细胞存活与炎症。该论文不能替代矽肺中的自噬通量证据，也说明“RORA 参与 CD4 炎症”本身不是新发现。<a href="https://doi.org/10.1016/j.celrep.2021.109682">Cell Reports 原始研究</a>。</li><li>2025 年肝纤维化研究是在 ATG5 缺失的 CD4 T 细胞中观察到 Rora 等 Th17 相关基因上调。它支持自噬状态与 T 细胞程序有关，不能逆推为 RORA 位于自噬上游。<a href="https://www.nature.com/articles/s41467-025-59218-y.pdf">Nature Communications 原始研究，Fig. 2d</a>。</li><li>心肌 RORα–CAV3 线粒体自噬、COPD RORα–p53–DRAM 属于其他细胞／疾病的机制证据。本批 CD4 的 CAV3 在全部 2,311 个细胞中计数为零，不能直接移植该轴；零检出不排除低于测序灵敏度的表达。<a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC8626585/">心肌研究</a>；<a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC9066967/">COPD 研究</a>。</li></ul><p>本次未检出完成矽肺 CD4–RORA–自噬因果链条的直接原始研究；这是一项检索范围内的观察，不是“从未报道”的保证。</p>''']
parts+=['<h2>5. 编号核查意见及验收要求</h2>',table(pd.DataFrame(issues))]
parts+=['''<h2>方法、文件定位与未验证事项</h2><p>输入为独立筛选 v1 的固定计数、细胞标签、程序成员／背景和 K15_R1 拟时序，所有源文件 SHA256 保存在 input_manifest.tsv；时间源文件记录在 time_input_manifest.tsv。RORA 来自原始计数，按原始文库量 CP10k 后 log1p。活化代理为 CD44、ICOS、CD38、TNFRSF4、TNFRSF9 的未加权均值；组织适应代理为 ITGAE、ITGA1、CXCR6、RGS1，均不含 RORA。它们只是本次混杂敏感性分析的代理变量，并非经过验证的分类器。</p><p>疾病差异：未调整均值与两套调整均值分别按六样本等权比较。调整斜率由各样本内中心化协变量估计，并使用每个样本相等的总权重；每个样本标准化到共同协变量均值。Core 额外控制亚群指示变量。每套模型在四个范围的比较内 BH。调整分析依赖线性／二次形式及共同斜率假设，不能证明因果独立或证明状态变量是混杂而非中介。</p><p>关联：仅使用六样本共同拟时序区间。每个范围／样本计算残差 Pearson r，最低样本数为全部六个；六个 Fisher z 的符号翻转依赖零假设下对称性。三个亚群×两个自噬程序组成六项主要检验家族；Core 两项另行校正；其他六程序以及两个状态代理另设家族，完整结果见 correlation_summary.tsv。各种分数／调整配置分别记录为敏感性分析，不以其中最小 P 选取结果。未新增疾病数据或蛋白实验，未验证自噬通量、调控靶点结合及救援效应。</p><p>机制检测表存在物种命名差异的严格符号检索结果，必须以 mechanism_alias_checked.tsv 的 TP53/TRP53、PRKN/PARK2 别名复核为准，不把未匹配人类符号写成小鼠基因缺失。原计数和旧报告均未修改。</p><p>服务器证据目录：/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/02_analysis/rora_followup_audit_20261016_v1。此文件为中文内部核查说明，不是客户版本或交付授权。</p></html>''']
out=F/'RORA_followup_audit_v1.html';out.write_text('\n'.join(parts),encoding='utf-8');soup=BeautifulSoup(out.read_text(encoding='utf-8'),'html.parser');assert len(soup.find_all('img'))==2 and len(soup.find_all('table'))==5 and '\ufffd' not in soup.get_text()
# Check core claims against source values and preserve exact output fingerprint.
assert a.query("target=='CAV3' and scope=='Core'").positive_cells.sum()==0
assert r.query('day==0').RORA_UMI.eq(0).all() and r.query('day>0').RORA_UMI.gt(0).all()
assert len(c.query("score_mode=='target_excluded' and model=='full' and module in ['Autophagy_GO','Mitophagy'] and P<.05"))==0
(F/'verification.json').write_text(json.dumps({'html_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'source_exit':0,'main_complete':True,'extra_complete':True,'numeric_claims_checked':True,'html_tables':5,'embedded_figures':2},indent=2),encoding='utf-8');print('REPORT_VERIFIED',out)
