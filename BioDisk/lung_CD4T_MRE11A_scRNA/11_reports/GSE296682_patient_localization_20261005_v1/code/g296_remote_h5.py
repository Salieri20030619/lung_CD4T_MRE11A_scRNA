from pathlib import Path
import io,subprocess,json,hashlib,time
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/GSE296682_patient_localization_20261005_v1')
class OfficialRangeFile(io.RawIOBase):
 def __init__(self):
  super().__init__();self.pos=0;self.size=642115063;self.block=262144;self.cache=O/'public_cache/h5ad_range_blocks';self.cache.mkdir(exist_ok=True);self.url='https://zenodo.org/records/18626573/files/BAL_FINAL.h5ad?download=1';self.requests=0;self.log=O/'provenance/range_blocks.jsonl'
 def readable(self):return True
 def seekable(self):return True
 def tell(self):return self.pos
 def seek(self,offset,whence=0):
  self.pos=int(offset if whence==0 else self.pos+offset if whence==1 else self.size+offset);return self.pos
 def getblock(self,i):
  start=i*self.block;end=min(self.size-1,start+self.block-1);p=self.cache/f'{start}_{end}.bin';expected=end-start+1
  if not p.exists():
   current=sum(x.stat().st_size for x in self.cache.glob('*.bin'));assert current+expected<=96*1024**2,'96MB bounded transfer budget reached; inspect before continuing'
   part=p.with_suffix('.part');r=subprocess.run(['curl','--fail','-L','--connect-timeout','15','--max-time','60','--max-filesize',str(expected),'--silent','--show-error','--range',f'{start}-{end}','--output',str(part),'--write-out','%{http_code}',self.url],capture_output=True,text=True,timeout=65)
   assert r.returncode==0 and r.stdout=='206',(r.returncode,r.stdout,r.stderr,start,end);assert part.stat().st_size==expected;part.rename(p);self.requests+=1
   with self.log.open('a') as f:f.write(json.dumps(dict(start=start,end=end,bytes=expected,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),url=self.url,http=206))+'\n')
   if self.requests%5==0:print('RANGE_PROGRESS',self.requests,'new_blocks','cached_MB',round((current+expected)/1024**2,2),flush=True)
  return p.read_bytes()
 def read(self,n=-1):
  if n<0:n=self.size-self.pos
  n=min(int(n),self.size-self.pos)
  if n<=0:return b''
  # Reuse already verified initial 64KB probe without an additional request.
  initial=O/'public_cache/BAL_FINAL.h5ad.range0_65535'
  if self.pos+n<=65536 and initial.exists():
   with initial.open('rb') as f:f.seek(self.pos);data=f.read(n)
   self.pos+=len(data);return data
  pieces=[]
  while n:
   i=self.pos//self.block;off=self.pos%self.block;b=self.getblock(i);take=min(n,len(b)-off);pieces.append(b[off:off+take]);self.pos+=take;n-=take
  return b''.join(pieces)
 def readinto(self,b):
  d=self.read(len(b));b[:len(d)]=d;return len(d)

if __name__=='__main__':
 import h5py
 f=OfficialRangeFile()
 with h5py.File(f,'r') as h:
  result=[]
  def visit(name,o):
   if name.count('/')<2:result.append(dict(path=name,kind=type(o).__name__,shape=str(getattr(o,'shape','')),dtype=str(getattr(o,'dtype','')),attrs={k:str(v) for k,v in o.attrs.items()}))
  h.visititems(visit)
 (O/'provenance/official_h5_range_schema.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
