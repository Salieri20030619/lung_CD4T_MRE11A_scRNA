from pathlib import Path
import os,json,time,subprocess,hashlib,shutil,stat
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/death_priority_20261018_v3')
D=Path('/mnt/BioData/bio/wangmingyuan/system_volume_archive_20261018');D.mkdir(exist_ok=True)
uid=os.getuid();active=set();maps=[]
for p in Path('/proc').iterdir():
 if not p.name.isdigit():continue
 try:
  if p.stat().st_uid!=uid:continue
  for f in [p/'cwd',p/'exe',*(p/'fd').iterdir()]:
   try:active.add(os.readlink(f).removesuffix(' (deleted)'))
   except OSError:pass
  for line in (p/'maps').read_text().splitlines():
   if '/' in line:active.add('/'+line.split('/',1)[1])
  for arg in (p/'cmdline').read_bytes().split(b'\0'):
   if arg.startswith(b'/'):active.add(arg.decode('utf-8'))
 except (OSError,UnicodeDecodeError):pass
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as hnd:
  for b in iter(lambda:hnd.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def manifest(p):
 fs=[p] if p.is_file() else list(p.rglob('*'))
 out={}
 for f in fs:
  st=f.lstat()
  if st.st_uid!=uid or not (stat.S_ISREG(st.st_mode) or stat.S_ISDIR(st.st_mode)):raise ValueError('mixed owner or special file')
  if stat.S_ISREG(st.st_mode):out[str(f.relative_to(p)) if f!=p else '.']=(st.st_size,st.st_mtime_ns,sha(f))
 return out
records=[];skipped=[]
for base in ['/tmp','/var/tmp']:
 for src in sorted(Path(base).iterdir()):
  try:
   if src.is_symlink() or src.lstat().st_uid!=uid:continue
   if not(src.is_file() or src.is_dir()):continue
   if time.time()-src.stat().st_mtime<86400:skipped.append((str(src),'recent'));continue
   if any(a==str(src) or a.startswith(str(src)+'/') for a in active):skipped.append((str(src),'active_reference'));continue
   before=manifest(src)
   if not before:continue
   size=sum(v[0] for v in before.values())
   if size<1024*1024:continue
   dst=D/src.relative_to('/');dst.parent.mkdir(parents=True,exist_ok=True)
   if dst.exists():skipped.append((str(src),'destination_exists'));continue
   if src.is_dir():subprocess.run(['rsync','-a','--',str(src)+'/',str(dst)+'/'],check=True)
   else:shutil.copy2(src,dst)
   assert manifest(dst)==before and manifest(src)==before
   backup=src.with_name(src.name+'.move_verified_20261018');assert not backup.exists()
   src.rename(backup)
   try:src.symlink_to(dst,target_is_directory=dst.is_dir())
   except BaseException:backup.rename(src);raise
   assert src.resolve()==dst.resolve() and manifest(dst)==before
   if backup.is_dir():shutil.rmtree(backup)
   else:backup.unlink()
   records.append({'source':str(src),'destination':str(dst),'bytes':size,'files':len(before),'verified':True})
   (O/'system_tmp_migration.json').write_text(json.dumps({'moved':records,'skipped':skipped},indent=2),encoding='utf-8')
   print('MOVED',str(src),size,flush=True)
  except (OSError,ValueError) as e:skipped.append((str(src),str(e)))
out={'status':'complete','moved':records,'skipped':skipped,'bytes_relocated':sum(v['bytes'] for v in records)}
(O/'system_tmp_migration.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print('COMPLETE',len(records),out['bytes_relocated'],flush=True)
