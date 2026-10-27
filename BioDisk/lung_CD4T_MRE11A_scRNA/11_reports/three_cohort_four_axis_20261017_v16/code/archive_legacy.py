from pathlib import Path
import json,zipfile,hashlib
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/three_cohort_four_axis_20261017_v16';records=[]
with zipfile.ZipFile(O/'resources/legacy_v13_evidence.zip','w',zipfile.ZIP_DEFLATED) as z:
 for i,f in enumerate(sorted((P/'11_reports/methods_complete_20261016_v13/supplement').glob('S17*.tsv')),1):
  name=f'S17_table_{i:02d}.tsv';z.writestr(name,f.read_bytes());records.append(dict(archive_name=name,original_path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
 z.writestr('original_paths_and_hashes.json',json.dumps(records,ensure_ascii=False,indent=2).encode('utf-8'))
print('LEGACY',len(records))