from pathlib import Path
import csv, gzip, io, json, tarfile, math
root=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
out=root/'02_analysis/screen6_G174_empty_audit'
out.mkdir(parents=True,exist_ok=True)
results=[]
with tarfile.open(root/'00_raw/GSE174725/GSE174725_RAW.tar') as t:
 for member in t:
  if not member.isfile() or not member.name.endswith('.csv.gz'):continue
  r=dict(file=member.name,rows=0,blank_rows=0,width_errors=0,empty_fields=0,other_bad_fields=0,examples=[],target_rows=[])
  with t.extractfile(member) as b,gzip.GzipFile(fileobj=b) as z,io.TextIOWrapper(z) as h:
   reader=csv.reader(h);header=next(reader)
   r['header_fields']=len(header);r['empty_header_columns']=[j for j,v in enumerate(header) if not v.strip()]
   for line,row in enumerate(reader,2):
    r['rows']+=1
    if not row or not any(v.strip() for v in row):r['blank_rows']+=1;continue
    if len(row)!=len(header):
     r['width_errors']+=1
     if len(r['examples'])<12:r['examples'].append(dict(line=line,kind='width',fields=len(row),head=row[:3]))
    empty=0;bad=0;positive=0;total=0
    for j,value in enumerate(row[1:],1):
     kind=None
     if not value.strip():empty+=1;kind='empty'
     else:
      try:
       v=float(value)
       if not math.isfinite(v):raise ValueError('nonfinite')
       positive+=int(v>0);total+=v
      except ValueError:bad+=1;kind='nonnumeric'
     if kind and len(r['examples'])<12:r['examples'].append(dict(line=line,gene=row[0],column=j,cell=header[j] if j<len(header) else '',kind=kind,value=value[:60]))
    r['empty_fields']+=empty;r['other_bad_fields']+=bad
    if row[0].strip().upper() in ['MRE11','MRE11A']:
     r['target_rows'].append(dict(gene=row[0],empty=empty,other_bad=bad,positive=positive,sum_observed=total))
  results.append(r);print(json.dumps(r,ensure_ascii=False),flush=True)
assert len(results)==5, f'Expected 5 patient CSVs, found {len(results)}'
(out/'empty_audit.json').write_text(json.dumps(results,ensure_ascii=False,indent=2))
with tarfile.open(root/'tmp/screen6_G174_empty_audit.tar.gz','w:gz') as t:t.add(out/'empty_audit.json',arcname='empty_audit.json')
print('G174_EMPTY_AUDIT_COMPLETE',flush=True)
