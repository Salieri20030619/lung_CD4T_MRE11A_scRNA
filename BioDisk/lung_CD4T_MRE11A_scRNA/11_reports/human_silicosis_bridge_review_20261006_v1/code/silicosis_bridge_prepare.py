from pathlib import Path
import datetime, gzip, hashlib, json, urllib.request

P = Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
O = P / '11_reports/human_silicosis_bridge_review_20261006_v1'
for d in ['code', 'provenance', 'public_inputs', 'results', 'figures']:
    (O / d).mkdir(parents=True, exist_ok=True)
spec = '''Frozen before reading GSE264182 expression; 2026-10-06.
Aim: bounded independent disease bridge, human BAL mixed-cell bulk RNA, not CD4-intrinsic or macrophage communication validation.
Only candidates IKZF2 and STAT1. No genome-wide screening, programs, SCENIC, or communication analysis in GSE264182.
Use submitted TMM log2 CPM values directly, never as integer counts; independent donor is unit.
Verify all 19 samples against GEO sample description/characteristics; SD2 -> SD02 only through explicit metadata mapping.
Primary contrast: all 14 silicosis (6 simple + 8 complicated) minus all 5 radiologically disease-free controls (3 silica-exposed + 2 chronic cough).
Sensitivity contrast: same 14 cases minus only the 3 silica-exposed disease-free controls.
Fixed method: unweighted difference of donor means; Welch two-sided t test and Welch-Satterthwaite t-based 95% CI, allowing unequal variances. Report sample sizes, means, variances, df, effect, CI, P and BH.
BH family: two candidate P values in primary; separate two-candidate family in exposure-control sensitivity. No alternative tests selected after results; no new permutation tests or cutoffs.
Simple and complicated groups: donor points and descriptive means only, no subtype hypothesis tests or favorable subgroup selection.
No age, smoking, batch, or cell-composition adjustment without actual donor-level metadata. All male.
All finite results, including null/contradictory results, are retained. Bulk mixture, small controls, and heterogeneous comparator types limit interpretation.
Separate GSE174725 diagnostic reuses frozen donor scores: Core only, IKZF2/STAT1 RNA plus all 11 existing frozen programs, All/NoFlags and min_cells 1/20 exactly as existing analysis. Compute leave-one-donor unweighted disease-minus-exposure mean differences only. Do not calculate new P/BH or new scores. Retain all results; groups with one remaining donor are descriptive only. Join original P/BH, check full effects reproduce source tests.
GSE296682 remains localization-only (two silicosis patients, no disease-free comparator).
No edits to existing PDF/ZIP; no addition/relabeling of original RNA19 samples.
'''
sp = O / 'provenance/FROZEN_SPEC.txt'
if sp.exists():
    assert sp.read_text() == spec, 'Existing specification differs: stop'
else:
    sp.write_text(spec)
    (O / 'provenance/spec_timestamp.json').write_text(json.dumps({'utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'sha256': hashlib.sha256(sp.read_bytes()).hexdigest(), 'expression_read': False}, indent=2))
print('FROZEN_SPEC', str(sp), flush=True)
receipts = []
for acc in ['GSE264182', 'GSE288887', 'GSE227136', 'GSE135893']:
    url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/' + acc[:-3] + 'nnn/' + acc + '/soft/' + acc + '_family.soft.gz'
    dest = O / 'public_inputs' / (acc + '_family.soft.gz')
    if not dest.exists():
        with urllib.request.urlopen(url, timeout=35) as h: raw = h.read(2000001)
        assert len(raw) <= 2000000
        gzip.decompress(raw)
        dest.write_bytes(raw)
    raw = dest.read_bytes()
    receipts.append(dict(file=dest.name, url=url, bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest()))
    txt = gzip.decompress(raw).decode()
    (O / 'public_inputs' / (acc + '_family.soft')).write_text(txt)
    print('GEO_DOWNLOADED', acc, len(raw), flush=True)
url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264182/suppl/GSE264182_log_normalized_counts_filtered_genes.csv.gz'
dest = O / 'public_inputs/GSE264182_log_normalized_counts_filtered_genes.csv.gz'
if not dest.exists():
    with urllib.request.urlopen(url, timeout=35) as h: raw = h.read(10000001)
    assert len(raw) <= 10000000
    gzip.decompress(raw)
    dest.write_bytes(raw)
raw = dest.read_bytes()
receipts.append(dict(file=dest.name, url=url, bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest()))
(O / 'provenance/download_manifest.json').write_text(json.dumps(receipts, indent=2))
print('EXPRESSION_DOWNLOADED', len(raw), flush=True)
for line in (O / 'public_inputs/GSE264182_family.soft').read_text().splitlines():
    if line.startswith(('^SAMPLE', '!Series_title', '!Series_summary', '!Series_overall_design', '!Sample_title', '!Sample_characteristics_ch1', '!Sample_description', '!Sample_data_processing')):
        print(line)
print('EXPRESSION_HEADER', gzip.decompress(raw).decode().splitlines()[0], flush=True)
