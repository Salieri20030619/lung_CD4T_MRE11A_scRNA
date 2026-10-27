from pathlib import Path
import json,os
O=Path('/mnt/BioData/bio/wangmingyuan/lung_CD4T_MRE11A_scRNA/12_scTenifoldKnk_20261018');P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/00_project')
status=json.loads((O/'delivery_complete.json').read_text());assert status['status']=='complete'
f=P/'LATEST_DEATH_PRIORITY.json';old=f.read_bytes();backup=O/'LATEST_DEATH_PRIORITY_before_v4.json'
if not backup.exists():backup.write_bytes(old)
a=dict(version='v4',status='complete',completion=str(O/'delivery_complete.json'),report=str(O/'delivery/report_v4.pdf'),baseline='v1/v2/v3 preserved',data_root='/mnt/BioData/bio/wangmingyuan/lung_CD4T_MRE11A_scRNA',scTenifoldKnk=dict(networks=24,single_knockouts=36,QA_checks=133,death_module_tests=72,death_module_nominal=0),previous=json.loads(old))
text=json.dumps(a,ensure_ascii=False,indent=2);(O/'LATEST_DEATH_PRIORITY.json').write_text(text,encoding='utf-8')
tmp=P/'LATEST_DEATH_PRIORITY_v4.tmp';tmp.write_text(text,encoding='utf-8');os.replace(tmp,f)
assert json.loads(f.read_text())['version']=='v4'
print('LATEST POINTER UPDATED TO V4')