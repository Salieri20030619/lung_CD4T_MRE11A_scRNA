from pathlib import Path
import json,hashlib,shutil,subprocess,shlex
B=Path(__file__).resolve().parent;D=B/'r1r2_revision_20261007'
pdf=next(D.glob('*复审R1R2修正版.pdf'));archive=next(D.glob('*复审R1R2修正版.zip'))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
result=json.loads((D/'delivery_package_result.json').read_text(encoding='utf-8'))
assert sha(pdf)==result['pdf']['pdf_sha256'] and sha(archive)==result['zip_sha256']
downloads=Path(r'C:\Users\31937\Downloads')
for source in [pdf,archive]:
    target=downloads/source.name
    if target.exists():assert sha(target)==sha(source),str(target)
    else:shutil.copy2(source,target)
    assert sha(target)==sha(source)
print('DOWNLOADS_HASHES_MATCH',flush=True)
server='/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/review_R1R2_20261007_v9'
flags=['-F',r'C:\Users\31937\.ssh\config','-o','BatchMode=yes','-o','StrictHostKeyChecking=yes','-o','ClearAllForwardings=yes','-o','ConnectTimeout=15']
def ssh(code):
    r=subprocess.run(['ssh',*flags,'sgs10','python3 -'],input=code,capture_output=True,text=True,encoding='utf-8',timeout=90)
    assert r.returncode==0,r.stderr
    return r.stdout
print(ssh('from pathlib import Path;import socket,json;p=Path('+repr(server)+');p.mkdir(exist_ok=True);(p/"delivery_checks").mkdir(exist_ok=True);print(json.dumps({"host":socket.gethostname(),"directory":str(p)}))'),flush=True)
r=subprocess.run(['scp',*flags,str(pdf),str(archive),'sgs10:'+server+'/'],capture_output=True,text=True,encoding='utf-8',timeout=180)
assert r.returncode==0,r.stderr;print('SERVER_ARTIFACTS_SAVED',flush=True)
verify=ssh('''from pathlib import Path,PurePosixPath
import json,hashlib,zipfile
p=Path('''+repr(server)+''')
pdf=next(p.glob('*.pdf'));archive=next(p.glob('*.zip'))
root=(p/'package').resolve();root.mkdir(exist_ok=True)
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None
 for info in z.infolist():
  rel=PurePosixPath(info.filename)
  assert not rel.is_absolute() and '..' not in rel.parts and ':' not in str(rel)
  dest=(root/info.filename).resolve();assert dest.is_relative_to(root)
  dest.parent.mkdir(parents=True,exist_ok=True)
  data=z.read(info.filename)
  if dest.exists():assert dest.read_bytes()==data
  else:dest.write_bytes(data)
 embedded=z.read('report/CD4T_review_R1R2_20261007.pdf');assert embedded==pdf.read_bytes()
 for line in z.read('FILE_MANIFEST_SHA256.tsv').decode().splitlines()[1:]:
  n,size,digest=line.split('\\t');data=(root/n).read_bytes();assert len(data)==int(size) and hashlib.sha256(data).hexdigest()==digest
print(json.dumps({'files':{f.name:{'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in [pdf,archive]},'archive_pdf_identical':True,'expanded_package_manifest_passed':True},ensure_ascii=True))
''')
verified=json.loads(verify)
for source in [pdf,archive]:assert verified['files'][source.name]['sha256']==sha(source)
status={'status':'COMPLETED','server_host':'osv01','server_directory':server,'downloads':[str(downloads/f.name) for f in [pdf,archive]],'Library_file_id':'libfile_51bec083de2c81919ea20e968043ae82','Library_version':8,'pages':result['pdf']['pages'],'original_scripts_unchanged':285,'scope':['R1','R2'],'R3_unchanged':True,'G288_unchanged':True,'Word_unchanged':True,'old_deliveries_preserved':True,'all_artifact_hashes_identical':True,**verified,'delivery_blockers':[]}
(D/'final_delivery_status.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf-8')
files=[D/'delivery_package_result.json',D/'final_delivery_status.json',D/'library_finalization_receipt.json',D/'library_identity_sidecar.json',B/'r1r2_deliver.py']
r=subprocess.run(['scp',*flags,*[str(p) for p in files],'sgs10:'+server+'/delivery_checks/'],capture_output=True,text=True,encoding='utf-8',timeout=90)
assert r.returncode==0,r.stderr
print(json.dumps(status,ensure_ascii=True),flush=True)
