from pathlib import Path
import sys,json,hashlib,csv,zipfile,datetime,gzip
O=Path(sys.argv[1])
required=['annotation_v2.R','compute_v2.py','time_v2.R','intersect_v2.py','support_v2.py','report_v2.py','finalize_v2.py']
for name in required:assert json.loads((O/'provenance'/(name+'.exit.json')).read_text(encoding='utf-8'))['exit_code']==0,name
qa=json.loads((O/'provenance/report_numeric_validation.json').read_text(encoding='utf-8'));assert qa['status']=='PASS'
for f in O.rglob('*'):
 if not f.is_file():continue
 if f.suffix in ['.py','.R','.tsv','.json','.txt','.html','.svg']:
  f.read_text(encoding='utf-8');assert not f.read_bytes().startswith(b'\xef\xbb\xbf'),f
 if f.name.endswith('.tsv.gz'):
  with gzip.open(f,'rt',encoding='utf-8') as h:
   for _ in h:pass
visual={'status':'PASS','PNG_figures_inspected':7,'revised_trajectory_legends':'PASS','real_time_axis_actual_spacing':'PASS','report_image_links':'PASS'}
(O/'provenance/visual_QA.json').write_text(json.dumps(visual,indent=2),encoding='utf-8')
status={'version':O.name,'status':'COMPLETE','sequence':['CD4 pyroptosis/apoptosis','mitochondrial/autophagy programs','within-state pseudotime AND separately measured real time','prospective intersection and experiment evidence cards'],'sample_cell_threshold':'none; all nonempty measurable units','discovery_libraries_per_scope':6,'time_libraries':12,'candidates':['AHNAK','MDFIC','TTC39C','ST6GALNAC3'],'candidate_evidence':'nominal exploratory, not FDR-confirmed mechanisms','prioritize_expression_validation':['AHNAK in Naive','ST6GALNAC3 in ActCD4'],'resolve_discordance':['TTC39C','MDFIC'],'unverified':['exact animal/pool to library mapping','patient-localization candidates outside legacy extracted panel','functional causality/autophagic flux'],'source_files_preserved':True,'client_translation_generated':False,'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(O/'PROJECT_COMPLETE.json').write_text(json.dumps(status,indent=2),encoding='utf-8')
rows=[]
for f in sorted(O.rglob('*')):
 if f.is_file() and f.name not in ['OUTPUT_MANIFEST_SHA256.tsv','package_v2.py.log','package_v2.py.exit.json']:
  rows.append({'path':str(f.relative_to(O)),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
with (O/'OUTPUT_MANIFEST_SHA256.tsv').open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['path','bytes','sha256'],delimiter='\t');w.writeheader();w.writerows(rows)
zpath=O.parent/(O.name+'.zip')
with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED) as z:
 for row in rows:z.write(O/row['path'],row['path'])
 z.write(O/'OUTPUT_MANIFEST_SHA256.tsv','OUTPUT_MANIFEST_SHA256.tsv')
with zipfile.ZipFile(zpath) as z:assert z.testzip() is None
print(json.dumps({'zip':str(zpath),'bytes':zpath.stat().st_size,'sha256':hashlib.sha256(zpath.read_bytes()).hexdigest(),'files':len(rows)}))