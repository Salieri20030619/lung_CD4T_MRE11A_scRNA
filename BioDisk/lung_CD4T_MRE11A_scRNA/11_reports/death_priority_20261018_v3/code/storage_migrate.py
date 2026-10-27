from pathlib import Path
import os,json,subprocess,hashlib,time,shutil,stat
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/death_priority_20261018_v3';D=Path('/mnt/BioData/bio/wangmingyuan/lung_CD4T_MRE11A_scRNA')
assert D.parent.resolve()==Path('/mnt/BioData/bio/wangmingyuan') and D.parent.stat().st_uid==os.getuid()
D.mkdir(exist_ok=True)
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def listing(p):
 out={}
 for f in p.rglob('*'):
  if f.is_symlink():out[str(f.relative_to(p))]={'symlink':os.readlink(f)}
  elif f.is_file():out[str(f.relative_to(p))]={'bytes':f.stat().st_size,'sha256':digest(f)}
 return out
records=[]
for rel in ['00_raw','11_reports/death_priority_20261018_v3/external']:
 src=P/rel;dst=D/rel
 if src.is_symlink():
  assert src.resolve()==dst.resolve();print('ALREADY_MOVED',src,flush=True);continue
 assert src.resolve().is_relative_to(P.resolve()) and dst.resolve().is_relative_to(D.resolve())
 dst.parent.mkdir(parents=True,exist_ok=True)
 before=listing(src);print('COPY',rel,len(before),sum(v.get('bytes',0) for v in before.values()),flush=True)
 subprocess.run(['rsync','-a','--',str(src)+'/',str(dst)+'/'],check=True)
 after=listing(dst);assert before==after,'destination verification failed'
 # Ensure no source mutation during copy.
 assert before==listing(src),'source changed during migration'
 backup=src.with_name(src.name+'.migration_verified_20261018');assert not backup.exists()
 src.rename(backup)
 try:src.symlink_to(dst,target_is_directory=True)
 except BaseException:backup.rename(src);raise
 assert src.resolve()==dst.resolve() and listing(src)==before
 records.append({'source':str(src),'destination':str(dst),'files':len(before),'bytes':sum(v.get('bytes',0) for v in before.values()),'sha256_verified':True,'manifest':before})
 (O/'storage_migration.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
 assert backup.resolve().is_relative_to(P.resolve()) and backup.name.endswith('.migration_verified_20261018')
 shutil.rmtree(backup)
 print('MOVED_AND_VERIFIED',rel,flush=True)
tmp=D/'tmp/death_priority_v3';tmp.mkdir(parents=True,exist_ok=True)
# Existing empty temp directory can be redirected; never relocate files in active use.
old=O/'tmp'
if not old.is_symlink() and not any(old.iterdir()):old.rmdir();old.symlink_to(tmp,target_is_directory=True)
print('TASK_TEMP',old.resolve(),flush=True)
# Read-only accounting of this user's system-volume files and deleted-open files.
rootdev=Path('/').stat().st_dev;own=[]
for base in ['/tmp','/var/tmp']:
 for f in Path(base).iterdir():
  try:
   if f.stat().st_uid!=os.getuid():continue
   if f.is_file():size=f.stat().st_size
   elif f.is_dir():size=int(subprocess.check_output(['du','-sb','--',str(f)],stderr=subprocess.DEVNULL).split()[0])
   else:continue
   own.append({'path':str(f),'bytes':size})
  except (OSError,subprocess.CalledProcessError):continue
deleted=[]
for pr in Path('/proc').iterdir():
 if not pr.name.isdigit():continue
 try:
  if pr.stat().st_uid!=os.getuid():continue
  for fd in (pr/'fd').iterdir():
   try:
    st=fd.stat();target=os.readlink(fd)
    if st.st_dev==rootdev and st.st_size>100*1024**2:deleted.append({'pid':int(pr.name),'fd':fd.name,'bytes':st.st_size,'target':target})
   except OSError:continue
 except OSError:continue
(O/'system_volume_readonly_audit.json').write_text(json.dumps({'own_tmp_bytes':sum(v['bytes'] for v in own),'largest_own_tmp':sorted(own,key=lambda x:-x['bytes'])[:20],'large_open_system_files':deleted},indent=2),encoding='utf-8')
print('SYSTEM_AUDIT',sum(v['bytes'] for v in own),deleted[:10],flush=True)
