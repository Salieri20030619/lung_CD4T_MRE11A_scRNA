from pathlib import Path
import ast, datetime, hashlib, json, shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

P = Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
O = P / '11_reports/human_silicosis_bridge_review_20261006_v1'
v = pd.read_csv(O / 'results/GSE264182_two_candidate_donor_values.tsv', sep='\t')
out = pd.read_csv(O / 'results/GSE264182_two_candidate_tests.tsv', sep='\t')
# Render the unchanged saved numbers with the final plot code, without refitting tests.
tree = ast.parse((O / 'code/silicosis_bridge_bulk.py').read_text())
nodes = []
active = False
for n in tree.body:
    segment = ast.get_source_segment((O / 'code/silicosis_bridge_bulk.py').read_text(), n) or ''
    if segment.startswith('plt.rcParams.update('): active = True
    if active: nodes.append(n)
    if segment == 'plt.close(fig)': break
exec(compile(ast.Module(body=nodes, type_ignores=[]), '<figure only from saved results>', 'exec'))
shutil.copy2(P / '11_reports/silicosis_bridge_prepare.py', O / 'code/silicosis_bridge_prepare.py')
g174 = pd.read_csv(O / 'results/GSE174725_Core_leave_one_summary.tsv', sep='\t')
assert len(out) == 4 and out.P.notna().all() and out.FDR_BH.notna().all()
assert (out.BH_family_size == 2).all()
assert (out.CI95_low < 0).all() and (out.CI95_high > 0).all()
assert len(v) == 19 and v.GSM.nunique() == 19
assert len(g174) == 52
srcs = [P / '11_reports/GSE296682_patient_localization_20261005_v1/provenance/final_interpretation.json',
        P / '11_reports/GSE296682_patient_localization_20261005_v1/results/fixed_candidate_two_patient_coverage_summary.tsv',
        P / '11_reports/GSE174725_annotation_validation_20261005_v1/results/frozen_scope_donor_counts.tsv']
evidence = [dict(path=str(f), bytes=f.stat().st_size, sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in srcs]
(O / 'provenance/additional_reused_evidence.json').write_text(json.dumps(evidence, indent=2))
assessment = {
    'utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'status': 'COMPLETE_BOUNDED_EVALUATION_STOP_EXPANSION',
    'GSE264182_tests': out.to_dict('records'),
    'GSE264182_nominal_significant': int((out.P < .05).sum()),
    'GSE264182_BH_significant': int((out.FDR_BH < .05).sum()),
    'GSE264182_interpretation_zh': '真实矽肺的独立BAL混合细胞疾病桥接已计算，但两基因均未获得明确差异支持。STAT1两个对照口径均为负向，暴露对照敏感性效应缩小；IKZF2接近零且依赖对照口径换向。所有区间跨零，不等于排除真实效应。',
    'GSE174725_interpretation_zh': 'Core STAT1 RNA下降方向在All/NoFlags及至少1/20细胞四规格的留一诊断均保持。IKZF2留一换向。现有程序全保留，不能把方向性描述当统计复现。',
    'GSE174725_descriptive_rows': {'summary': 52, 'leave_one': 234},
    'GSE174725_pairing_zh': '成熟BALF_AM与Core两端均至少20细胞仅1位暴露对照，无矽肺病例；单核细胞不能替代成熟巨噬细胞。',
    'GSE296682_interpretation_zh': '两位矽肺患者无疾病阴性对照；仅补定位可检出性。作者CD4标签存在低深度/谱系混合不确定，不能验证疾病差异或完整受体复合物。',
    'disease_labels_zh': '原RNA19仍为9位CWP与10位对照，非纯矽肺。CTL117已包含且最终CD4=461；缺少CTL105的直接依据是最终CD4=0，早期QC有LOW_RETENTION标记。',
    'recommended_patch_zh': '建议以独立矽肺队列的有限验证附表/短段说明补充数据来源，并如实保留无显著支持及跨对照口径敏感性。不要把主队列改为矽肺，不主张扩大阈值或程序搜索。',
    'cannot_patch_zh': ['CD4亚群内在疾病机制', 'IKZF2/STAT1调控子活性', '成熟巨噬细胞到CD4的同供者疾病通讯链', 'RNA变化的因果性'],
    'PDF_ZIP_modified': False,
    'Library_uploaded': False,
    'new_SCENIC_communication_program_screen': False,
    'blockers': [],
    'limitations': ['GSE264182 mixed BAL and limited control donors; no age/smoking/batch/cell proportion covariates', 'GSE174725 3+2 donors and uneven cell coverage', 'GSE296682 2 patients without controls'],
    'validation': ['Strict 19-donor matrix/metadata bijection including explicit SD2/SD02 correction', '14407 by 19 matrix; candidates unique and finite', 'Welch formula agrees with scipy.ttest_ind(equal_var=False)', 'BH m=2 per prespecified contrast', 'All 52 existing Core full effects reproduce source to 1e-12', 'All original RNA19 diagnoses verified from official GEO', 'Final point figure checked, no model changes during rendering']
}
(O / 'provenance/ASSESSMENT.json').write_text(json.dumps(assessment, indent=2, ensure_ascii=False))
for f in (O / 'code').glob('*.py'): compile(f.read_text(), str(f), 'exec')
manifest = []
for f in sorted(O.rglob('*')):
    if f.is_file() and f.name != 'OUTPUT_MANIFEST.tsv':
        manifest.append(dict(path=str(f.relative_to(O)), bytes=f.stat().st_size, sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
pd.DataFrame(manifest).to_csv(O / 'OUTPUT_MANIFEST.tsv', sep='\t', index=False)
print('COMPLETE',len(manifest),'files',sum(r['bytes'] for r in manifest),'bytes',O)
print('ASSESSMENT_SHA256',hashlib.sha256((O / 'provenance/ASSESSMENT.json').read_bytes()).hexdigest())
