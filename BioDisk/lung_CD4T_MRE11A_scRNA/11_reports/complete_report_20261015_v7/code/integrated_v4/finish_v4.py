from pathlib import Path
import sys,json,hashlib,csv,zipfile
import pandas as pd
O=Path(sys.argv[1]);checks=[]
for name in ['provenance/mouse_input_manifest.tsv','05_human/input_manifest.tsv']:
 d=pd.read_csv(O/name,sep='\t');print(name,list(d.columns))
 pathcol=next(c for c in ['path','source','file'] if c in d.columns);hashcol=next(c for c in ['sha256','SHA256'] if c in d.columns)
 for _,r in d.iterrows():
  p=Path(r[pathcol]);assert p.is_file(),str(p);h=hashlib.sha256(p.read_bytes()).hexdigest();assert h==r[hashcol],str(p);checks.append(str(p))
status={p.name:json.loads(p.read_text())['exit_code'] for p in (O/'provenance').glob('*.exit.json')};assert all(v==0 for v in status.values()),status
(O/'provenance/final_integrity.json').write_text(json.dumps({'status':'PASS','unchanged_source_files':len(checks),'recorded_script_exit_codes':status,'scope':'new v4 outputs and referenced source files; causal chain remains unconfirmed'},indent=2),encoding='utf-8')
(O/'PROJECT_COMPLETE.json').write_text(json.dumps({'execution':'complete','version':'20261015_v4','internal_report':'REPORT_INTERNAL.html','scientific_outcome':'partial exploratory evidence; original ligand-TF-pyroptosis chain not established','trajectory_required':False,'old_results_overwritten':False},indent=2),encoding='utf-8')
print('SOURCE_INTEGRITY_PASS',len(checks))