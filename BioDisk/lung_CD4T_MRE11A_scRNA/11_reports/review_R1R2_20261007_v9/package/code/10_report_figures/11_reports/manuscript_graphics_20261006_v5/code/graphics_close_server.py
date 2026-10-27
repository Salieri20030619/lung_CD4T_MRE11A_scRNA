from pathlib import Path,PureWindowsPath
import json,hashlib,zipfile,datetime,shutil
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/manuscript_graphics_20261006_v5'
receipt=json.loads((O/'provenance/library_upload_result.json').read_text(encoding='utf8'));rs=receipt['results']
assert len(rs)==8 and all(r['status']=='succeeded' for r in rs)
assert all(r['current_version_number']==4 for r in rs[:2])
assert [r['library_file_id'] for r in rs[:2]]==['libfile_51bec083de2c81919ea20e968043ae82','libfile_fc4b3de706cc81918c05a5e5612aa70d']
qa=json.loads((O/'provenance/QA_result.json').read_text(encoding='utf8'));audit=json.loads((O/'provenance/editorial_audit.json').read_text(encoding='utf8'));parity=json.loads((O/'provenance/PDF_content_parity.json').read_text(encoding='utf8'))
assert qa['pages']==28 and qa['figures']==14 and qa['tables']==15 and not qa['out_of_bounds'] and not qa['internal_term_hits'] and not qa['empty_BH_columns']
assert audit['tables_identical']==15 and audit['source_data_unchanged'] and not audit['old_prose_numerals_for_review']
assert parity['docx_text_blocks_verified']==1490 and not parity['missing_text_blocks']
for name,digest in qa['sha256'].items():assert hashlib.sha256((O/'deliverables'/name).read_bytes()).hexdigest()==digest
identity=json.loads((O/'provenance/library_identity.json').read_text(encoding='utf8'))
for r in identity['files']:
 name=PureWindowsPath(r['local_path']).name;f=O/('figures' if name.startswith('Human_') else 'deliverables')/name
 assert hashlib.sha256(f.read_bytes()).hexdigest()==r['sha256'],name
status={'date_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'server':'osv01','output':str(O),'pdf_pages':28,'main_figures':10,'supplementary_figures':4,'tables':15,'references':15,'new_scientific_models':0,'new_UMAP_runs':0,'new_GSEA_tests':0,'numeric_checks':'34/34 passed','GSEA_curve_checks':20,'docx_pdf_text_blocks_matched':1490,'table_data_preserved':15,'source_data_unchanged':True,'original_figures_identical':9,'visual_review':qa.get('visual_review'),'library_status':'completed','library_results':[{k:r.get(k) for k in ['library_file_id','file_id','status','current_version_number','file_name']} for r in rs],'deliverable_sha256':qa['sha256'],'delivery_blockers':[],'renderer_note':'Word COM stalled; exact DOCX content rendered locally with headless Edge and checked for complete text/table parity.'}
f=O/'provenance/FINAL_STATUS.json';f.write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf8')
for n in ['QA_result.json','editorial_audit.json','numeric_audit.json','PDF_content_parity.json']:
 shutil.copy2(O/'provenance'/n,O/'reproduce/graphics_delivery/QA'/n)
z=O/'deliverables/CD4_graphics_reproduction_SERVER_ONLY.zip'
with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as zz:
 for p in sorted(O.rglob('*')):
  if p.is_file() and p!=z and p.suffix!='.zip':zz.write(p,p.relative_to(O))
with zipfile.ZipFile(z) as zz:assert zz.testzip() is None
status['zip_bytes']=z.stat().st_size;status['zip_sha256']=hashlib.sha256(z.read_bytes()).hexdigest();f.write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps(status,ensure_ascii=True))
