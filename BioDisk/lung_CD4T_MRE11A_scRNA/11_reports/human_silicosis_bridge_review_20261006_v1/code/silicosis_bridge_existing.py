from pathlib import Path
import hashlib, json, shutil
import numpy as np
import pandas as pd

P = Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
O = P / '11_reports/human_silicosis_bridge_review_20261006_v1'
assert (O / 'provenance/FROZEN_SPEC.txt').exists()
R = P / '11_reports/GSE174725_annotation_validation_20261005_v1/results'
src = R / 'fixed_validation_donor_scores_SERVER_ONLY.tsv.gz'
testsrc = R / 'all_fixed_validation_tests.tsv'
d = pd.read_csv(src, sep='\t')
t = pd.read_csv(testsrc, sep='\t')
d = d[(d.scope == 'Core') & ((d.family == 'receiver_program') | ((d.family == 'receiver_RNA') & d.endpoint.isin(['IKZF2', 'STAT1'])))].copy()
assert d.endpoint.nunique() == 13
assert not d.duplicated(['sample', 'cell_filter', 'family', 'endpoint']).any()
rows, summaries = [], []
for (filt, fam, ep), g in d.groupby(['cell_filter', 'family', 'endpoint']):
    for cut in [1, 20]:
        v = g[(g.cells >= cut) & np.isfinite(g.value)]
        x = v.loc[v.group == 'Exposed', 'value'].values
        y = v.loc[v.group == 'Silicosis', 'value'].values
        full = float(y.mean() - x.mean())
        old = t[(t.scope == 'Core') & (t.cell_filter == filt) & (t.family == fam) & (t.endpoint == ep) & (t.min_cells == cut)]
        assert len(old) == 1
        assert np.isclose(full, old.iloc[0].effect, rtol=0, atol=1e-12), (filt, ep, cut)
        tag = dict(scope='Core', cell_filter=filt, min_cells=cut, family=fam, endpoint=ep)
        effects = []
        for sample in v['sample']:
            w = v[v['sample'] != sample]
            xx = w.loc[w.group == 'Exposed', 'value'].values
            yy = w.loc[w.group == 'Silicosis', 'value'].values
            effect = float(yy.mean() - xx.mean()) if len(xx) and len(yy) else np.nan
            effects.append(effect)
            rows.append(dict(tag, omitted=sample, n_exposed=len(xx), n_silicosis=len(yy), full_effect=full, leave_one_effect=effect, status='DESCRIPTIVE_ONLY_NO_NEW_INFERENCE'))
        e = np.array(effects)
        summaries.append(dict(tag, n_exposed=len(x), n_silicosis=len(y), effect=full,
                              source_P_exact=old.iloc[0].P_exact, source_FDR_BH=old.iloc[0].FDR_BH,
                              source_BH_family_size=old.iloc[0].BH_family_size,
                              leave_one_min=float(np.nanmin(e)), leave_one_max=float(np.nanmax(e)),
                              leave_one_positive=int((e > 0).sum()), leave_one_negative=int((e < 0).sum()),
                              leave_one_total=int(np.isfinite(e).sum()), all_leave_one_same_sign=bool(np.all(np.sign(e) == np.sign(full))),
                              case_control_pairs_lower=int((y[:, None] < x[None, :]).sum()),
                              case_control_pairs_higher=int((y[:, None] > x[None, :]).sum()),
                              case_control_pairs_total=len(x)*len(y),
                              interpretation='Leave-one and pairwise signs are descriptive; pairwise contrasts are not independent observations'))
out = pd.DataFrame(summaries)
out.to_csv(O / 'results/GSE174725_Core_leave_one_summary.tsv', sep='\t', index=False)
pd.DataFrame(rows).to_csv(O / 'results/GSE174725_Core_leave_one_details.tsv', sep='\t', index=False)
assert len(out) == 52 and len(rows) == 234
srcs = [src, testsrc, P / '11_reports/human_alignment_20261005_v1/results/stage1_donor_design.tsv']
(O / 'provenance/reused_input_manifest.json').write_text(json.dumps([dict(path=str(p), bytes=p.stat().st_size, sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in srcs], indent=2))
shutil.copy2(srcs[2], O / 'results/original_RNA19_donor_design.tsv')
print('G174_LATEST', len(out), len(rows))
print(out[(out.cell_filter == 'All') & (out.min_cells == 1)].drop(columns=['scope', 'cell_filter', 'min_cells', 'interpretation']).to_string(index=False))
print('TF_SENSITIVITIES')
print(out[out.family == 'receiver_RNA'].drop(columns=['interpretation']).to_string(index=False))
