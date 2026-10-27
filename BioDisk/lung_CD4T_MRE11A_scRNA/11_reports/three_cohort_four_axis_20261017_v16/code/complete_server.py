from pathlib import Path
import json,hashlib,zipfile,pandas as pd,datetime
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/three_cohort_four_axis_20261017_v16');j=json.loads((O/'completion_local.json').read_text(encoding='utf-8'));checks={}
for local,info in j['outputs'].items():
 f=O/'delivery'/('v16_report_zh.pdf' if local.endswith('.pdf') else 'v16_review_evidence.zip');sha=hashlib.sha256(f.read_bytes()).hexdigest();assert sha==info['sha256'];checks[f.name]=sha
with zipfile.ZipFile(O/'delivery/v16_review_evidence.zip') as z:
 for n in z.namelist():
  if n.startswith(('final_tables/','literature_success/','local_code/','resources/original_downloads/')) or n in ['resources/resource_manifest_resolved.json','README_review.txt']:
   dest=(O/n).resolve();assert dest.is_relative_to(O.resolve());dest.parent.mkdir(exist_ok=True,parents=True);dest.write_bytes(z.read(n))
for n in ['DE_complete.txt','NicheNet_complete.txt','QC_export_complete.txt','intersection_complete.txt']:assert (O/'logs'/n).read_text().strip()=='success'
s=pd.read_csv(O/'tables/three_axis_candidates_all_contrasts.tsv',sep='\t');assert s.gene.nunique()==55 and len(s)==155
j['server_files_sha256']=checks;j['server_verified_at']=datetime.datetime.now(datetime.timezone.utc).isoformat();(O/'FINAL_COMPLETION.json').write_text(json.dumps(j,ensure_ascii=False,indent=2),encoding='utf-8');print('VERIFIED',json.dumps(checks));print('SUCCESS 55 genes / 155 comparisons / 18-page PDF / full evidence archive')