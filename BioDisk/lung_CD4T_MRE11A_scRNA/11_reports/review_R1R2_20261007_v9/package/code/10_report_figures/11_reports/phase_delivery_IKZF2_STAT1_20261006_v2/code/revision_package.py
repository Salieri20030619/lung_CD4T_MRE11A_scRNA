from pathlib import Path
import json,hashlib,zipfile
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/phase_delivery_IKZF2_STAT1_20261006_v2'
expected={'CD4_IKZF2_STAT1_stage_delivery_zh.pdf':'38616059799547fc94019169d804cf0773ef521cfe415c5fd299e55cb4974c8c','CD4_IKZF2_STAT1_stage_delivery_zh.docx':'94f8e4f48031af42a2ad4cf4f19edb6c941f48aeac3e2644b0a6edd7c2a77e7b'}
for name,want in expected.items():
 p=O/'deliverables'/name;got=hashlib.sha256(p.read_bytes()).hexdigest();assert got==want,(name,got)
qa=json.loads((O/'provenance/QA_render.json').read_text());assert qa['pages']==26 and qa['automated_bounds_pass'] and not qa['empty_fdr_columns'] and qa['visual_review'].startswith('PASS')
files=[f for root in ['code','results','provenance','deliverables','figures'] for f in (O/root).rglob('*') if f.is_file() and f.suffix!='.zip' and f.name not in ['final_file_manifest.json','FINAL_STATUS.json']]
manifest=[dict(path=str(f.relative_to(O)),size=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in files]
m=O/'provenance/final_file_manifest.json';m.write_text(json.dumps(manifest,indent=2));files.append(m)
z=O/'deliverables/CD4_IKZF2_STAT1_revision_SERVER_ONLY.zip'
with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED) as a:
 for f in files:a.write(f,f.relative_to(O))
with zipfile.ZipFile(z) as a:assert a.testzip() is None;count=len(a.namelist())
library=json.loads((O/'provenance/library_final_receipt_identity.json').read_text())
assert len(library)==2 and all(x['status']=='succeeded' and x['current_version_number']==1 for x in library)
r=dict(status='COMPLETE',pages=26,tables=40,figures=7,references=12,new_models=0,new_network_runs=0,new_bootstrap_runs=0,automatic_followup='STOPPED',server_directory=str(O),reports=expected,library=[{k:x[k] for k in ['library_file_id','file_name','current_version_number','status']} for x in library],zip=str(z),zip_files=count,zip_size=z.stat().st_size,zip_sha256=hashlib.sha256(z.read_bytes()).hexdigest(),raw_data_uploaded=False)
(O/'provenance/FINAL_STATUS.json').write_text(json.dumps(r,indent=2));print(json.dumps(r))
