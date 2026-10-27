from pathlib import Path
import re, hashlib, json, shutil
import pandas as pd

P = Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
O = P / '11_reports/human_silicosis_bridge_review_20261006_v1'
design = pd.read_csv(O / 'results/original_RNA19_donor_design.tsv', sep='\t')
rows = []
for acc, source in [('GSE135893', 'G135'), ('GSE227136', 'G227'), ('GSE288887', 'G288')]:
    txt = (O / 'public_inputs' / (acc+'_family.soft')).read_text()
    for block in txt.split('^SAMPLE = ')[1:]:
        lines = block.splitlines()
        title = next(x.split(' = ',1)[1] for x in lines if x.startswith('!Sample_title = '))
        char = [x.split(' = ',1)[1] for x in lines if x.startswith('!Sample_characteristics_ch1 = ')]
        diag = '; '.join(x for x in char if x.startswith(('disease state', 'diagnosis')))
        mid = re.search(r'VU_?(HD|ILD)_?0*(\d+)', title, flags=re.I)
        if not mid:
            continue
        sample = ('CTL' + mid.group(2).zfill(3)) if mid.group(1).upper() == 'HD' else ('CWP' + str(int(mid.group(2))))
        match = design[(design['sample'] == sample) & (design.series == source)]
        if len(match) != 1:
            continue
        if source == 'G288':
            assert 'lung explants from subjects with coal-workers pneumoconiosis' in txt
            observed = 'CWP; series overall design (individual sample characteristics list tissue only)'
        elif match.iloc[0]['group'] == 'CWP':
            assert 'CWP' in diag
            observed = diag
        else:
            assert 'control' in diag.lower()
            observed = diag
        rows.append(dict(sample=sample, source=source, GSM=lines[0], title=title, submitted_diagnosis=observed,
                         current_group=match.iloc[0]['group'], match='VERIFIED',
                         note='Multiple regions/libraries from one donor remain one donor',
                         source_url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+acc))
out = pd.DataFrame(rows)
assert set(out['sample']) == set(design['sample']), set(design['sample'])-set(out['sample'])
out.to_csv(O / 'results/original_RNA19_GEO_diagnosis_audit.tsv', sep='\t', index=False)
countfile = P / '10_tables/CD4_final/CD4_cells_by_subject.tsv'
c = pd.read_csv(countfile, sep='\t')
assert c.loc[c.subject_id == 'CTL105', 'final_cd4_cells'].iloc[0] == 0
assert c.loc[c.subject_id == 'CTL117', 'final_cd4_cells'].iloc[0] == 461
shutil.copy2(countfile, O / 'provenance/original_final_CD4_cells_by_subject.tsv')
meta = P / '03_data_audit/corrected_v2/GSE288887_sample_metadata.tsv'
q = P / '03_data_audit/formal_preprocessing_v1/GSE288887_formal_sample_qc.tsv'
for f in [meta, q]:
    d = pd.read_csv(f, sep='\t')
    mask = d.astype(str).apply(lambda col: col.str.contains('CTL105|CTL117|VUHD105|VUHD117', regex=True)).any(axis=1)
    d[mask].to_csv(O / 'provenance' / (f.stem+'_105_117.tsv'), sep='\t', index=False)
conclusion = dict(current_RNA19_source_counts=design.groupby(['series', 'group']).size().to_dict().__str__(),
    total_current_donors=19, current_controls=10, current_CWP=9,
    G227_controls=8, G227_CWP=2, G135_controls=2, G288_CWP=7,
    CTL117=dict(in_current_RNA19=True, final_CD4_cells=461, original_sample_type='MNC'),
    CTL105=dict(in_current_RNA19=False, final_CD4_cells=0, original_group='Control', QC_flag='LOW_RETENTION',
               interpretation='Verified reason for absence from final CD4 roster: zero final CD4 cells. Earlier low-retention QC flag is documented; do not claim it was the sole cause or an explicit donor-level exclusion decision.'),
    retained_diagnosis='Coal workers pneumoconiosis (CWP), not relabeled as pure silicosis',
    mutations='No original donor/sample was added, removed, or relabeled',
    sources=[dict(path=str(f), sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in [countfile, meta, q]])
(O / 'provenance/original_RNA19_roster_audit.json').write_text(json.dumps(conclusion, indent=2))
print(json.dumps(conclusion, indent=2))
print('ALL_19_DIAGNOSES_VERIFIED', len(out), 'sample/library metadata records')
