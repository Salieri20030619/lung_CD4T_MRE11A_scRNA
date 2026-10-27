from pathlib import Path
import json,hashlib,os
O=Path('/mnt/BioData/bio/wangmingyuan/lung_CD4T_MRE11A_scRNA/18_macrophage_supplement_20261019')
c=json.loads((O/'completion.json').read_text(encoding='utf-8'))
assert hashlib.sha256((O/'delivery/report_v10.pdf').read_bytes()).hexdigest()==c['sha256']
assert hashlib.sha256((O/'delivery/v10_supplement_materials.zip').read_bytes()).hexdigest()==c['materials_sha256']
assert json.loads((O/'logs/mouse_complete.json').read_text())['status']=='success'
assert (O/'logs/targeted_nichenet_complete.txt').read_text().strip()=='success'
assert json.loads((O/'logs/silicosis_BAL_complete.json').read_text())['status']=='completed_coverage_and_descriptive_only'
p=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/00_project/LATEST_DEATH_PRIORITY.json')
old=json.loads(p.read_text(encoding='utf-8'));(O/'previous_LATEST_DEATH_PRIORITY.json').write_text(json.dumps(old,ensure_ascii=False,indent=2),encoding='utf-8')
new=dict(old);new.update(version='v10_mainline',report=str(O/'delivery/report_v10.pdf'),report_pdf=str(O/'delivery/report_v10.pdf'),previous_version='v9_mainline',macrophage_supplement=str(O),human_validation_scope='silicosis_only; CWP excluded',supplement_status='complete',report_pages=37,report_sha256=c['sha256'])
# Preserve legacy pending analyses instead of implying that an unrelated gap was completed.
new['pending_trajectory_gene_GO_KEGG']=True
tmp=p.with_suffix('.v10.tmp');tmp.write_text(json.dumps(new,ensure_ascii=False,indent=2),encoding='utf-8');os.replace(tmp,p)
(O/'server_verified.json').write_text(json.dumps(dict(status='passed',report_sha256=c['sha256'],archive_sha256=c['materials_sha256'],latest_pointer=str(p),previous_version_preserved=True),indent=2),encoding='utf-8')
print('SERVER_ARCHIVE_AND_POINTER_VERIFIED')
