from pathlib import Path
import sys
P=Path(sys.argv[1]).parents[1];J=P/'11_reports/GSE296682_patient_localization_20261005_v1'
for f in ['provenance/official_fixed_expression_provenance.json','provenance/metadata_verified.json','code/g296_extract_range.py','code/g296_remote_h5.py']:
 print(f);print((J/f).read_text()[:13500])