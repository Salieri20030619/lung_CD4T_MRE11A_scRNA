from pathlib import Path
import shutil,json,hashlib
B=Path(__file__).parent;D=B/'graphics_delivery';D.mkdir(exist_ok=True)
for x in ['provenance','figures','QA']:(D/x).mkdir(exist_ok=True)
for n in ['added_evidence_aggregate.json','old_blocks.json']:
 shutil.copy2(B/'narrative_delivery'/n,D/n)
shutil.copytree(B/'narrative_delivery/figures',D/'figures',dirs_exist_ok=True)
files=['narrative_delivery/CD4_IKZF2_STAT1_integrated_manuscript_zh.docx','narrative_delivery/CD4_IKZF2_STAT1_integrated_manuscript_zh.pdf','narrative_delivery/added_evidence_aggregate.json','revision_delivery/report_evidence_aggregate.json']
(D/'provenance/baseline_v3.json').write_text(json.dumps({f:hashlib.sha256((B/f).read_bytes()).hexdigest() for f in files},indent=2),encoding='utf8')
p=(B/'narrative_results.py').read_text(encoding='utf8')
print(p[p.index("h('2.10"):p.index("h('2.11")].encode('ascii','backslashreplace').decode())
