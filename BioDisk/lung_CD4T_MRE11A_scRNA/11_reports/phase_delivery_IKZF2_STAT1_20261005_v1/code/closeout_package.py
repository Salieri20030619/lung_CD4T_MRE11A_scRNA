from pathlib import Path
import json,hashlib,zipfile
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/phase_delivery_IKZF2_STAT1_20261005_v1'
files=[f for root in ['code','results','provenance','deliverables','figures'] for f in (O/root).rglob('*') if f.is_file() and f.suffix!='.zip']
manifest=[dict(path=str(f.relative_to(O)),size=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in files]
(O/'provenance/final_file_manifest.json').write_text(json.dumps(manifest,indent=2))
z=O/'deliverables/CD4_IKZF2_STAT1_closeout_SERVER_ONLY.zip'
with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED) as a:
 for f in files:a.write(f,f.relative_to(O))
with zipfile.ZipFile(z) as a:assert a.testzip() is None;count=len(a.namelist())
r=dict(status='COMPLETE',new_models=0,automatic_followup='STOPPED',server_directory=str(O),zip=str(z),zip_files=count,zip_size=z.stat().st_size,zip_sha256=hashlib.sha256(z.read_bytes()).hexdigest(),raw_data_uploaded=False)
(O/'provenance/FINAL_STATUS.json').write_text(json.dumps(r,indent=2));print(json.dumps(r))
