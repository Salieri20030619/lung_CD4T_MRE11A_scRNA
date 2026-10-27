from pathlib import Path
import shutil,json,hashlib
B=Path(__file__).parent;D=B/'narrative_delivery';D.mkdir(exist_ok=True);(D/'figures').mkdir(exist_ok=True)
for name in ['report_evidence_aggregate.json','document_manifest.json']:
 pass
for f in (B/'manuscript_delivery/figures').glob('*.png'):shutil.copy2(f,D/'figures'/f.name)
for name in ['added_evidence_aggregate.json','figure_label_audit.json']:
 shutil.copy2(B/'manuscript_delivery'/name,D/name)
baseline={'source_document_sha256':hashlib.sha256((B/'manuscript_delivery/CD4_IKZF2_STAT1_integrated_manuscript_zh.docx').read_bytes()).hexdigest(),'source_version':2,'source_tables_sha256':{name:hashlib.sha256(path.read_bytes()).hexdigest() for name,path in [('E',B/'revision_delivery/report_evidence_aggregate.json'),('F',B/'manuscript_delivery/added_evidence_aggregate.json')]},'new_scientific_models':0}
(D/'baseline.json').write_text(json.dumps(baseline,indent=2),encoding='utf8')
print(json.dumps({'directory':str(D),'baseline_saved':True}))
