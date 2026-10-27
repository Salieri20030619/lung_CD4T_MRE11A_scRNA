from pathlib import Path
import json,hashlib,zipfile,datetime
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/manuscript_integrated_20261006_v3'
receipt=json.loads((O/'provenance/library_upload_result.json').read_text(encoding='utf8'))
results=receipt['results'];assert len(results)==2 and all(r['status']=='succeeded' and r['current_version_number']==2 for r in results)
qa=json.loads((O/'provenance/QA_result.json').read_text(encoding='utf8'))
assert qa['pages']==26 and not qa['out_of_bounds'] and not qa['internal_term_hits'] and not qa['empty_BH_columns']
for name,digest in qa['sha256'].items():assert hashlib.sha256((O/'deliverables'/name).read_bytes()).hexdigest()==digest
status={'date_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'server':'osv01','output':str(O),'pages':26,'main_figures':7,'supplementary_figures':1,'main_tables':1,'supplementary_tables':14,'references':15,'new_scientific_models':0,'numeric_checks':'34/34 passed','internal_terms':0,'text_overflow':0,'library_results':[{k:r.get(k) for k in ['library_file_id','file_id','status','current_version_number','file_name','local_metadata_applied']} for r in results],'remaining_provenance_gaps':['Human MSigDB release number not independently recovered; archived members retained','Human legacy UCell maxRank not explicitly frozen in original call'],'skill_errors':['documents/PDF package not available in delegated environment; Word and PyMuPDF QA completed']}
(O/'provenance/FINAL_STATUS.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf8')
z=O/'deliverables/CD4_integrated_reproduction_server_only.zip'
with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as zz:
 for f in sorted(O.rglob('*')):
  if f.is_file() and f!=z and f.suffix!='.zip':zz.write(f,f.relative_to(O))
status['zip_bytes']=z.stat().st_size;status['zip_sha256']=hashlib.sha256(z.read_bytes()).hexdigest()
(O/'provenance/FINAL_STATUS.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf8');print(json.dumps(status,ensure_ascii=True))
