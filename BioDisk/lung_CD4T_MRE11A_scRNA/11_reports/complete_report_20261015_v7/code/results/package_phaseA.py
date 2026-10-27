from pathlib import Path
import json,hashlib,csv,zipfile,datetime
O=Path(__file__).resolve().parent.parent
assert json.loads((O/'provenance/exit_status.json').read_text(encoding='utf-8'))['exit_code']==0
assert json.loads((O/'provenance/finalize.py.exit.json').read_text(encoding='utf-8'))['exit_code']==0
assert json.loads((O/'provenance/final_validation.json').read_text(encoding='utf-8'))['status']=='PASS'
# Check read-back dimensions and source integrity, no new calculations.
for r in csv.DictReader((O/'tables/input_manifest_SHA256.tsv').open(encoding='utf-8'),delimiter='\t'):
 assert hashlib.sha256(Path(r['path']).read_bytes()).hexdigest()==r['sha256']
status={'phase':'A','status':'COMPLETE','components':['input_and_reuse_inventory','historical_screening_flow','Core_ActCD4_volcanoes','historical_six_evidence_and_sample_points','within_subtype_program_comparisons','independent_statistical_and_visual_validation'],'new_target_selection':'NOT_STARTED','customer_delivery':'NOT_GENERATED','GSE288887_primary_discovery':'EXCLUDED','PRJNA1289989_download':'NOT_STARTED','time_support_blockers':['GSE274372 D0/D28 CD4/TLike coverage under existing cell gate','exact animal/pool-to-library mapping unresolved'],'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(O/'PHASE_A_COMPLETE.json').write_text(json.dumps(status,indent=2),encoding='utf-8')
readme='''Plan 1 phase A: internal working outputs

Open phaseA_working_report.html. Figures are relative links under figures/.

Completed: input/reuse inventory; unchanged 20260914 historical DE reconstruction;
Core/ActCD4 volcanoes; historical-six effect/P/BH and library-point figures;
Naive/ActCD4/Treg within-state disease changes for seven frozen transcript panels.

Actual statistics are in tables/phenotype_group_comparisons.tsv.
Missing P/BH means descriptive-only or insufficient replication, never P=1.
BG1 Treg has only one control library at the 20-cell gate.
Library-to-animal/pool mapping is unresolved; new inference is exploratory.

Key tables:
- cohort_roles_and_reuse.tsv, required_report_reuse.tsv
- input_manifest_SHA256.tsv, later_phase_input_inventory.tsv
- historical_screening_counts.tsv, historical_six_DE_all_scopes.tsv
- phenotype_donor_inputs.tsv, phenotype_group_comparisons.tsv
- subtype_differential_response.tsv, BH_family_inventory.tsv

Validation:
8 frozen DE count checks; 77 exact permutation results independently checked
against SciPy; 6 valid BH families checked; six PNG figures visually inspected.
Original source hashes rechecked unchanged. User-supplied ZIP PDF/HTML hashes
match its 20261007_R1R2 manifest. All scripts and run exits are preserved.
First attempt stopped at an overly strict assumption that all six genes had
been tested in every subtype. Missing Naive results were explicitly restored as
not-tested records; no P or BH was imputed.

Reproduction (writes only this revision directory):
1. batch_bench Python code/phaseA.py
2. batch_bench Python code/finalize.py
3. batch_bench Python code/package_phaseA.py
Existing upstream inputs remain external dependencies; no raw matrices bundled.

Next-stage gene intersection, pseudotime and real-time selection have NOT run.
This package is an internal checkpoint, not a Chinese customer report or final
target shortlist. Do not represent historical integrated judgement as an
algorithmically derived intersection.
'''
(O/'README.txt').write_text(readme,encoding='utf-8')
rows=[]
for f in sorted(O.rglob('*')):
 if f.is_file() and f.name not in ['OUTPUT_MANIFEST_SHA256.tsv','package_phaseA.py.log','package_phaseA.py.exit.json']:
  rows.append({'path':str(f.relative_to(O)),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
with (O/'OUTPUT_MANIFEST_SHA256.tsv').open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['path','bytes','sha256'],delimiter='\t');w.writeheader();w.writerows(rows)
zpath=O.parent/(O.name+'.zip')
with zipfile.ZipFile(zpath,'w',compression=zipfile.ZIP_DEFLATED) as z:
 for row in rows:z.write(O/row['path'],row['path'])
 z.write(O/'OUTPUT_MANIFEST_SHA256.tsv','OUTPUT_MANIFEST_SHA256.tsv')
with zipfile.ZipFile(zpath) as z:assert z.testzip() is None
print(json.dumps({'zip':str(zpath),'bytes':zpath.stat().st_size,'sha256':hashlib.sha256(zpath.read_bytes()).hexdigest(),'manifest_files':len(rows)}))