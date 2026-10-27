from pathlib import Path
import json,hashlib,zipfile,subprocess,shutil
B=Path(__file__).resolve().parent;D=B/'r1r2_revision_20261007';K=D/'package'
K.mkdir(parents=True,exist_ok=True)
old=B/'logic_revision_20261007';baseline=next(old.glob('*最终逻辑核查版.zip'))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(baseline)=='de47b5d5e22f93f3c1f700007064186c09cdb336d6952ceff01ca16d8fe2f1eb'
with zipfile.ZipFile(baseline) as z:
    for n in z.namelist():
        dest=(K/n).resolve();assert dest.is_relative_to(K.resolve()) and not Path(n).is_absolute()
        dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(z.read(n))
v5=K/'report_build/history/library_v5';v7=K/'report_build/history/library_v7';cur=K/'report_build/current'
for p in [v5/'QA',v7/'QA',v7/'inputs',v7/'outputs',v7/'code',cur/'inputs',cur/'sources',cur/'QA',cur/'code']:p.mkdir(parents=True,exist_ok=True)
for n in ['PDF_content_parity.json','QA_result.json','numeric_revision_checks.json','revision_manifest.json','package_summary.json']:
    shutil.copy2(K/'provenance'/n,v5/'QA'/n)
for n in ['README_zh.md','documentation/运行顺序与依赖.md']:
    target=v5/Path(n).name;shutil.copy2(K/n,target)
shutil.copy2(K/'provenance/package_validation.json',v7/'QA/package_validation_from_v7_zip.json')
shutil.copy2(B/'review_delivery/provenance/review_render_source.html',v7/'inputs/review_render_source.html')
shutil.copy2(old/'provenance/original_blocks.json',v7/'inputs/original_blocks.json')
shutil.copy2(old/'provenance/editorial_changes.json',v7/'outputs/editorial_changes.json')
shutil.copy2(old/'provenance/logic_render_source.html',v7/'outputs/logic_render_source.html')
shutil.copy2(next(old.glob('*最终逻辑核查版.pdf')),v7/'outputs/library_v7.pdf')
for n in ['logic_edit_report.py','logic_render_pdf.py','logic_finalize_package.py']:shutil.copy2(B/n,v7/'code'/n)
for f in (old/'QA').rglob('*'):
    if f.is_file() and f.name!='renderer_log.txt':
        target=v7/'QA'/f.relative_to(old/'QA');target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(f,target)
# Read-only source extraction, keyed by full endpoint/model identity. No model fit or BH calculation.
code=r'''
from pathlib import Path
import json,hashlib,socket
import pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
f=P/'11_reports/human_first_20261004_v1/results/stage2_evidence_all_candidates.tsv.gz'
d=pd.read_csv(f,sep='\t',low_memory=False)
selector={'filter':'All','model':'Source','test':'Nominal05','direction':'Up','ligand':'IL15'}
q=d[d.receiver.isin(['Core','CD4-IL7R']) & d.sender.isin(['MacAM','MonoMac','MacWide'])].copy()
for key,val in selector.items():q=q[q[key].eq(val)]
cols=['receiver','sender',*selector,'aupr_corrected','rank','status','association_n','association_n_case','association_n_control','association_rho','association_P','association_FDR','association_status']
assert len(q)==6 and not q.duplicated(['receiver','sender',*selector]).any()
S=P/'11_reports/review_logic_20261007_v8'
paths=['provenance/editorial_changes.json','provenance/logic_render_source.html','provenance/original_blocks.json','QA/PDF_content_parity.json','QA/final_logic_audit.json','code/logic_edit_report.py','code/logic_render_pdf.py','code/logic_finalize_package.py']
out={'host':socket.gethostname(),'source_file':str(f),'source_sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'source_bytes':f.stat().st_size,'selectors':selector,'rows':json.loads(q[cols].to_json(orient='records',double_precision=15)), 'server_v7_inventory':[{'relative_path':x,'bytes':(S/x).stat().st_size,'sha256':hashlib.sha256((S/x).read_bytes()).hexdigest()} for x in paths]}
print(json.dumps(out))
'''
flags=['-F',r'C:\Users\31937\.ssh\config','-o','BatchMode=yes','-o','StrictHostKeyChecking=yes','-o','ClearAllForwardings=yes','-o','ConnectTimeout=15']
r=subprocess.run(['ssh',*flags,'sgs10','/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python -'],input=code,text=True,encoding='utf-8',capture_output=True,timeout=60)
assert r.returncode==0,r.stderr
e=json.loads(r.stdout);assert e['host']=='osv01' and e['source_sha256']=='6d4cfde6e328c298f20c9b6cda97c3101758f4662400502e4827769d5d30a9b3'
(cur/'inputs/IL15_table2_source_rows.json').write_text(json.dumps(e,ensure_ascii=False,indent=2),encoding='utf-8')
for record in e['server_v7_inventory']:
    rel=record['relative_path'];local=(B/Path(rel).name) if rel.startswith('code/') else old/rel
    assert local.stat().st_size==record['bytes'] and sha(local)==record['sha256'],rel
print('BASELINE_STAGED',len(list(K.rglob('*'))),'SOURCE_ROWS',len(e['rows']),'SERVER_V7_MATCHED',len(e['server_v7_inventory']),flush=True)
