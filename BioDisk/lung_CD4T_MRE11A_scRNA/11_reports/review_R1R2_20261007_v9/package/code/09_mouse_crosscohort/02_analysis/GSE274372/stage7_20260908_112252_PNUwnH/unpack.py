import os,gzip,shutil
from pathlib import Path
root=Path(os.environ['ROOT']);run=Path(os.environ['RUN'])
src=root/'00_raw/GSE274372/GSE274372_cds_object.rds.gz'
if not src.is_file(): raise FileNotFoundError(src)
current=src
for level in range(5):
 with current.open('rb') as f: magic=f.read(12)
 print('AUTHOR_LAYER=',level,' MAGIC_HEX=',magic.hex(),flush=True)
 if magic[:2]!=b'\x1f\x8b': break
 if level==4: raise RuntimeError('More than four gzip layers')
 dest=run/f'author_unpacked_{level}.bin'
 with gzip.open(current,'rb') as fi, dest.open('wb') as fo:
  shutil.copyfileobj(fi,fo,length=1024*1024)
 current=dest
# Only a run-local copy is passed to R. R handles any remaining supported compression.
shutil.copyfile(current,run/'author_payload.bin')
print('AUTHOR_UNPACK_PASS',flush=True)
