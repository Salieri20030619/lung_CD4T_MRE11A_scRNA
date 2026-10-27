from pathlib import Path
import subprocess,shutil
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/GSE174725_reaudit_20261005_v1')
B=O/'provenance/first_attempt_before_blank_row_fix';B.mkdir(exist_ok=True)
for n in ['CSV_structure_summary.tsv','count_vector_comparison_all_suffixes.tsv','paired_lineage_feasibility.tsv','existing_Core_count_integrity.tsv']:
 p=O/'results'/n
 if p.exists() and not (B/n).exists():shutil.copy2(p,B/n)
(B/'failure.txt').write_text('First attempt stopped at GSM5324868 finite-value assertion because a documented all-empty CSV row was parsed as NaN. No expression imputed. Corrected parser skips only a row with missing gene and all missing values; asserts exactly one such row for this sample. First two sample aggregate results preserved. Validation rerun once to save complete per-cell trace/checkpoints and validate parser correction.\n')
with (O/'provenance/validation.log').open('w') as f:
 p=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python',str(O/'code/gse174725_validate.py')],stdout=f,stderr=subprocess.STDOUT)
print('EXIT',p.returncode);print((O/'provenance/validation.log').read_text()[-7000:]);raise SystemExit(p.returncode)
