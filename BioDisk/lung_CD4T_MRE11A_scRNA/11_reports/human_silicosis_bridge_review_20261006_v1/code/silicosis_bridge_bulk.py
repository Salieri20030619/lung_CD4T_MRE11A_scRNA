from pathlib import Path
import datetime, hashlib, json, platform, sys
import numpy as np
import pandas as pd
import scipy
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

P = Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
O = P / '11_reports/human_silicosis_bridge_review_20261006_v1'
assert (O / 'provenance/FROZEN_SPEC.txt').exists()
rawmeta = O / 'public_inputs/GSE264182_family.soft'
expr = O / 'public_inputs/GSE264182_log_normalized_counts_filtered_genes.csv.gz'
records, rec = [], None
for line in rawmeta.read_text().splitlines():
    if line.startswith('^SAMPLE = '):
        rec = {'GSM': line.split(' = ', 1)[1]}
        records.append(rec)
    elif rec is not None and line.startswith(('!Sample_title = ', '!Sample_description = ')):
        key, value = line.split(' = ', 1)
        rec[key.removeprefix('!Sample_')] = value
    elif rec is not None and line.startswith('!Sample_characteristics_ch1 = '):
        key, value = line.split(' = ', 1)[1].split(': ', 1)
        rec[key] = value
m = pd.DataFrame(records)
assert len(m) == 19 and m.GSM.nunique() == 19
assert set(m.Sex) == {'Male'}
assert m['group'].value_counts().to_dict() == {'Complicated Silicosis': 8, 'Simple Silicosis': 6, 'No Radiological Disease': 5}
m['matrix_id'] = m.description.replace({'SD2': 'SD02'})
assert m.matrix_id.nunique() == 19
m['mapping_rule'] = np.where(m.description == 'SD2', 'Explicit leading-zero correction SD2 -> SD02', 'Exact GEO description')
m['subgroup'] = m['group'].replace({'Simple Silicosis': 'Simple', 'Complicated Silicosis': 'Complicated'})
m.loc[m.title.str.contains('Chronic Cough', regex=False), 'subgroup'] = 'CoughCtrl'
m.loc[m.title.str.contains('Silica Exposed No Disease', regex=False), 'subgroup'] = 'ExposedCtrl'
assert m.subgroup.value_counts().to_dict() == {'Complicated': 8, 'Simple': 6, 'ExposedCtrl': 3, 'CoughCtrl': 2}
m['disease'] = m.subgroup.isin(['Simple', 'Complicated'])
m['source'] = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE264182'
m.to_csv(O / 'results/GSE264182_sample_mapping.tsv', sep='\t', index=False)
d = pd.read_csv(expr, index_col=0)
assert d.shape == (14407, 19), d.shape
assert d.index.is_unique and d.columns.is_unique
assert set(d.columns) == set(m.matrix_id), (set(d.columns) - set(m.matrix_id), set(m.matrix_id) - set(d.columns))
assert set(['IKZF2', 'STAT1']).issubset(set(d.index)), d.index[:5].tolist()
v = m.merge(d.loc[['IKZF2', 'STAT1']].T.rename_axis('matrix_id').reset_index(), on='matrix_id', validate='one_to_one')
assert np.isfinite(v[['IKZF2', 'STAT1']].to_numpy()).all()
v.to_csv(O / 'results/GSE264182_two_candidate_donor_values.tsv', sep='\t', index=False)
rows = []
for contrast in ['Primary_AllControls', 'Sensitivity_ExposedControls']:
    case = v[v.disease]
    ctrl = v[~v.disease] if contrast == 'Primary_AllControls' else v[v.subgroup == 'ExposedCtrl']
    assert (len(case), len(ctrl)) == (14, 5 if contrast == 'Primary_AllControls' else 3)
    for gene in ['IKZF2', 'STAT1']:
        x, y = ctrl[gene].to_numpy(), case[gene].to_numpy()
        vx, vy = x.var(ddof=1), y.var(ddof=1)
        a, b = vx/len(x), vy/len(y)
        se = np.sqrt(a+b)
        df = (a+b)**2 / (a*a/(len(x)-1) + b*b/(len(y)-1))
        effect = y.mean() - x.mean()
        statistic = effect / se
        p = 2 * stats.t.sf(abs(statistic), df)
        check = stats.ttest_ind(y, x, equal_var=False)
        assert np.isclose(p, check.pvalue, atol=1e-14)
        ci = stats.t.ppf(.975, df) * se
        rows.append(dict(contrast=contrast, gene=gene, n_silicosis=len(y), n_control=len(x),
                         silicosis_mean_log2CPM=float(y.mean()), control_mean_log2CPM=float(x.mean()),
                         silicosis_variance=float(vy), control_variance=float(vx),
                         effect_log2CPM=float(effect), CI95_low=float(effect-ci), CI95_high=float(effect+ci),
                         standard_error=float(se), Welch_df=float(df), Welch_t=float(statistic), P=float(p),
                         method='Welch two-sided t test; Welch-Satterthwaite t 95% CI',
                         effect_definition='Silicosis donor mean minus comparator donor mean on submitted log2 CPM scale'))
out = pd.DataFrame(rows)
for contrast, g in out.groupby('contrast'):
    pv = g.P.to_numpy()
    order = np.argsort(pv)
    adj = np.minimum.accumulate((pv[order]*len(pv)/np.arange(1, len(pv)+1))[::-1])[::-1]
    restored = np.empty(len(pv)); restored[order] = np.minimum(adj, 1)
    out.loc[g.index, 'FDR_BH'] = restored
out['BH_family_size'] = 2
out['BH_family'] = out.contrast + '_IKZF2_STAT1'
out.to_csv(O / 'results/GSE264182_two_candidate_tests.tsv', sep='\t', index=False)
des = v.melt(id_vars=['GSM', 'matrix_id', 'subgroup'], value_vars=['IKZF2', 'STAT1'], var_name='gene', value_name='log2CPM').groupby(['gene', 'subgroup']).agg(n=('log2CPM', 'size'), mean=('log2CPM', 'mean'), sd=('log2CPM', 'std'), minimum=('log2CPM', 'min'), maximum=('log2CPM', 'max')).reset_index()
des['interpretation'] = 'Descriptive only; no subgroup hypothesis tests'
des.to_csv(O / 'results/GSE264182_subgroups_descriptive.tsv', sep='\t', index=False)
plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':10, 'axes.spines.top':False, 'axes.spines.right':False, 'pdf.fonttype':42})
fig, axes = plt.subplots(1, 2, figsize=(11, 5.8))
groups = ['CoughCtrl', 'ExposedCtrl', 'Simple', 'Complicated']
colors = ['#64748b', '#0891b2', '#e79a23', '#c65356']
for ax, gene in zip(axes, ['IKZF2', 'STAT1']):
    for k, (group, color) in enumerate(zip(groups, colors)):
        z = v[v.subgroup == group].sort_values(gene)
        jitter = np.linspace(-.16, .16, len(z))
        ax.scatter(k+jitter, z[gene], c=color, s=42, edgecolors='white', linewidth=.6, zorder=3)
        ax.plot([k-.24, k+.24], [z[gene].mean()]*2, color=color, linewidth=2.4)
    ax.set_xticks(range(4), ['CoughCtrl\nn=2', 'ExposedCtrl\nn=3', 'Simple\nn=6', 'Complicated\nn=8'])
    ax.set_xlim(-.5, 3.6)
    ax.set_ylabel('Submitted TMM log2 CPM')
    ax.set_title(gene, loc='left', fontweight='bold', fontsize=14)
    ax.grid(axis='y', color='#e2e8f0', linewidth=.7, zorder=0)
    pp = out[(out.gene == gene) & (out.contrast == 'Primary_AllControls')].iloc[0]
    ss = out[(out.gene == gene) & (out.contrast == 'Sensitivity_ExposedControls')].iloc[0]
    ax.text(0, -.25, f'Primary 14 vs 5: Δ={pp.effect_log2CPM:+.3f}; P={pp.P:.3g}; BH={pp.FDR_BH:.3g}\nExposed 14 vs 3: Δ={ss.effect_log2CPM:+.3f}; P={ss.P:.3g}; BH={ss.FDR_BH:.3g}', transform=ax.transAxes, fontsize=9, va='top')
fig.suptitle('GSE264182 | Independent silicosis cohort · mixed BAL cells', fontsize=14, y=.98)
fig.text(.5, .9, 'Each point is one donor; horizontal lines show subgroup means.', ha='center', fontsize=10, color='#475569')
fig.subplots_adjust(top=.82, bottom=.28, left=.075, right=.98, wspace=.3)
for suffix in ['png', 'pdf']:
    fig.savefig(O / ('figures/GSE264182_two_candidates.'+suffix), dpi=220, bbox_inches='tight')
plt.close(fig)
(O / 'figures/GSE264182_two_candidates_caption_zh.txt').write_text('每点为独立供者，短横线为分层均值。CoughCtrl：无影像疾病的慢性咳嗽对照（2）；ExposedCtrl：硅暴露但无影像疾病（3）；Simple：单纯矽肺（6）；Complicated：复杂矽肺（8）。纵轴为作者提供的 TMM log2 CPM。主比较为全部矽肺14人对全部无影像疾病5人，敏感性为同14例对3位暴露无病对照。双侧Welch检验，两个候选各比较内单独构成BH家族（m=2）；效应为矽肺减对照的供者均值差。分层仅展示，未选择有利亚组检验。该图为BAL混合细胞证据，不能定位到CD4内在变化或验证巨噬细胞通讯。')
audit = dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), dimensions=list(d.shape), donors=19,
             metadata_groups=m.subgroup.value_counts().to_dict(), candidates=['IKZF2', 'STAT1'], tests=4,
             missing_candidate_values=0, mapped_one_to_one=True, correction='GSM8213669: SD2 -> SD02',
             software=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__, matplotlib=matplotlib.__version__),
             metadata_sha256=hashlib.sha256(rawmeta.read_bytes()).hexdigest(), expression_sha256=hashlib.sha256(expr.read_bytes()).hexdigest(),
             caveats=['Mixed BAL cells: cell composition unavailable', 'All male', 'No donor-level age, smoking, batch or cell proportion adjustment', 'Controls are heterogeneous; exposure-only sensitivity has n=3', 'Welch inference relies on small-sample distributional assumptions; broad CIs retained', 'RNA abundance is not TF activity, causality, or cell-specific mechanism'])
(O / 'provenance/GSE264182_analysis_audit.json').write_text(json.dumps(audit, indent=2))
print(out.to_string(index=False))
print(des.to_string(index=False))
print('COMPLETED', O, flush=True)
