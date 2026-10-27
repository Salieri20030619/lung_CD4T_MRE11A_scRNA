cat > /mnt/BioDisk/wangmingyuan/screen6_structure.py <<'PY'
from pathlib import Path
import gzip,tarfile,json,csv
import h5py
root=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
out=root/'02_analysis/screen6_structure';out.mkdir(parents=True,exist_ok=True)
records=[]
def preview(f):
 lines=[]
 for _ in range(4):
  raw=f.readline()
  if not raw:break
  s=raw.decode('utf-8',errors='replace').rstrip('\r\n');parts=s.split('\t') if '\t' in s else s.split(',')
  lines.append({'fields':len(parts),'head':s[:1200],'tail':s[-300:] if len(s)>1200 else ''})
 return lines
def dec(v):
 return v.decode(errors='replace') if isinstance(v,bytes) else str(v)
for g in ['GSE184854','GSE183682','GSE296682','GSE295463','GSE174725','GSE159354']:
 rows=[]
 for p in sorted((root/'00_raw'/g).iterdir()):
  if not p.is_file() or p.name.endswith('.part') or '.invalid_' in p.name:continue
  r={'file':p.name,'bytes':p.stat().st_size}
  try:
   if p.name.endswith('.tar'):
    with tarfile.open(p) as t:
     members=[m for m in t if m.isfile()];r['members']=[{'name':m.name,'bytes':m.size} for m in members]
     r['previews']={}
     for m in members[:6]:
      if m.name.endswith(('.gz','.txt','.tsv','.csv','.mtx')):
       f=t.extractfile(m)
       if m.name.endswith('.gz'):f=gzip.GzipFile(fileobj=f)
       r['previews'][m.name]=preview(f);f.close()
   elif p.name.endswith(('.h5ad','.h5')):
    with h5py.File(p,'r') as h:
     r['keys']=list(h.keys());r['groups']={}
     for key in ['obs','var','raw/var','layers','X','raw/X']:
      if key not in h:continue
      z=h[key];r['groups'][key]={'keys':list(z.keys()) if isinstance(z,h5py.Group) else [],'shape':str(z.attrs.get('shape',getattr(z,'shape','')))}
      if key in ['obs','var','raw/var']:
       r['groups'][key]['columns']={}
       for name,d in z.items():
        if isinstance(d,h5py.Dataset):
         item={'shape':str(d.shape),'head':[dec(v) for v in d[:3]] if d.ndim==1 else []}
         if key!='obs' and d.ndim==1 and d.dtype.kind in 'OSU':
          vals=[dec(v) for v in d[:]];item['target_matches']=[v for v in vals if v.upper() in ['MRE11','MRE11A','ENSMUSG00000031928','ENSG00000020922']]
         r['groups'][key]['columns'][name]=item
        elif isinstance(d,h5py.Group) and 'categories' in d:
         vals=d['categories'];r['groups'][key]['columns'][name]={'categories':[dec(v) for v in vals[:100]],'n_categories':len(vals)}
   elif p.name.endswith(('.gz','.txt','.csv','.tsv','.mtx')) and 'family.soft' not in p.name:
    with (gzip.open(p,'rb') if p.name.endswith('.gz') else p.open('rb')) as f:r['preview']=preview(f)
  except Exception as e:r['inspection_error']=str(e)
  rows.append(r)
 (out/(g+'_structure.json')).write_text(json.dumps(rows,ensure_ascii=False,indent=2))
 records.append({'dataset':g,'files':len(rows),'inspection_errors':sum('inspection_error' in r for r in rows)})
 print(g,records[-1],flush=True)
(out/'summary.json').write_text(json.dumps(records,indent=2))
with tarfile.open(root/'tmp/screen6_structure_return.tar.gz','w:gz') as t:
 for p in sorted(out.glob('*.json')):t.add(p,arcname=p.name)
print('SCREEN6_STRUCTURE_PASS' if not any(r['inspection_errors'] for r in records) else 'SCREEN6_STRUCTURE_REVIEW_ERRORS',flush=True)
PY
conda run -n batch_bench python -m py_compile /mnt/BioDisk/wangmingyuan/screen6_structure.py &&
conda run --no-capture-output -n batch_bench python -u /mnt/BioDisk/wangmingyuan/screen6_structure.py
