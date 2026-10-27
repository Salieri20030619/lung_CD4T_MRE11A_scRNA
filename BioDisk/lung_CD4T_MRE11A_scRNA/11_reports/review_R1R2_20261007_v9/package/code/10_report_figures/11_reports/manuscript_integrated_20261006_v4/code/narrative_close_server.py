from pathlib import Path
import json, hashlib, zipfile, datetime

P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
O=P/'11_reports/manuscript_integrated_20261006_v4'
assert O.resolve().is_relative_to(P.resolve())
receipt=json.loads((O/'provenance/library_upload_result.json').read_text(encoding='utf8'))
results=receipt['results']
assert len(results)==2 and all(r['status']=='succeeded' and r['current_version_number']==3 for r in results)
assert {r['library_file_id'] for r in results}=={'libfile_51bec083de2c81919ea20e968043ae82','libfile_fc4b3de706cc81918c05a5e5612aa70d'}
qa=json.loads((O/'provenance/QA_result.json').read_text(encoding='utf8'))
audit=json.loads((O/'provenance/editorial_audit.json').read_text(encoding='utf8'))
assert qa['pages']==25 and not qa['out_of_bounds'] and not qa['internal_term_hits'] and not qa['empty_BH_columns']
assert audit['tables_identical']==15 and audit['source_data_unchanged'] and not audit['old_prose_numerals_for_review']
for name,digest in qa['sha256'].items():
    assert hashlib.sha256((O/'deliverables'/name).read_bytes()).hexdigest()==digest
status={
    'date_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'server':'osv01','output':str(O),'pages':25,
    'main_figures':8,'supplementary_figures':1,'main_tables':1,'supplementary_tables':14,'references':15,
    'new_scientific_models':0,'numeric_checks':'34/34 passed','table_data_preserved':15,
    'prior_prose_numerals_preserved':True,'source_data_unchanged':True,
    'visual_review':qa.get('visual_review'),'cross_reference_review':qa.get('cross_reference_review'),
    'internal_terms':0,'text_overflow':0,'library_status':'completed',
    'library_results':[{k:r.get(k) for k in ['library_file_id','file_id','status','current_version_number','file_name','local_metadata_applied']} for r in results],
    'deliverable_sha256':qa['sha256'],
    'remaining_provenance_gaps':['Human MSigDB release number not independently recovered; archived members retained','Human legacy UCell maxRank not explicitly frozen in original call'],
    'delivery_blockers':[]
}
status_path=O/'provenance/FINAL_STATUS.json'
status_path.write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf8')
z=O/'deliverables/CD4_narrative_reproduction_server_only.zip'
with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as zz:
    for f in sorted(O.rglob('*')):
        if f.is_file() and f!=z and f.suffix!='.zip':
            zz.write(f,f.relative_to(O))
with zipfile.ZipFile(z) as zz:
    assert zz.testzip() is None
    assert json.loads(zz.read('provenance/FINAL_STATUS.json'))['library_status']=='completed'
status['zip_bytes']=z.stat().st_size
status['zip_sha256']=hashlib.sha256(z.read_bytes()).hexdigest()
status_path.write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps(status,ensure_ascii=True))
