import os,json,tarfile,hashlib
from pathlib import Path
import pandas as pd
import numpy as np
out=Path(os.environ['OUT']);root=Path(os.environ['ROOT']);p6=Path(os.environ['PREV6'])
a=pd.read_csv(out/'tables/human_TF_activity.tsv',sep='\t');e=pd.read_csv(out/'tables/human_TF_effects.tsv',sep='\t');r=pd.read_csv(out/'tables/human_TF_associations.tsv',sep='\t')
assert len(a)>0 and {'Core','Primary'}<=set(a.scope)
assert not a.duplicated(['scope','filter','level','tf','sample']).any()
assert np.isfinite(a.activity).all()
assert (r.rho.dropna().abs()<=1+1e-10).all()
assert not e.duplicated(['scope','filter','level','tf','model']).any()
assert set(r[r.kind=='Module'].exclusion)=={'Drop7PlusModule'}
# Mouse and human activity estimates are not numerically pooled.
m=pd.read_csv(p6/'tables/TF_group_effects.tsv',sep='\t');m=m[(m.scope=='ActCD4')&(m['mode']=='Drop7')][['filter','level','tf','delta','LOO_min','LOO_max']].rename(columns={'delta':'mouse_delta','LOO_min':'mouse_LOO_min','LOO_max':'mouse_LOO_max'})
z=e.merge(m,on=['filter','level','tf'],how='left',validate='many_to_one');ok=np.isfinite(z.coefficient)&np.isfinite(z.mouse_delta)
z['same_direction']=pd.Series(pd.NA,index=z.index,dtype='boolean');z.loc[ok,'same_direction']=np.sign(z.loc[ok,'coefficient'])==np.sign(z.loc[ok,'mouse_delta'])
z.to_csv(out/'tables/mouse_human_TF_direction.tsv',sep='\t',index=False,na_rep='NA')
status=dict(status='PASS',stage='mechanism_step7_human_TF',output=str(out),human_subjects=int(a['sample'].nunique()),human_scopes=sorted(a.scope.unique()),TFs=sorted(a.tf.unique()),pending=['human_macrophage_object_verification_and_analysis','final_mechanism_literature_review'])
(out/'status.json').write_text(json.dumps(status,indent=2))
archive=root/'tmp/lung_mechanism_step7_return.tar.gz';part=archive.with_suffix('.tar.gz.part')
with tarfile.open(part,'w:gz') as t:
 for p in sorted(out.rglob('*')):
  if p.is_file() and not {'objects','__pycache__'}&set(p.relative_to(out).parts):t.add(p,arcname=str(p.relative_to(out)))
with tarfile.open(part) as t:
 for m in t:
  if m.isfile():
   f=t.extractfile(m)
   while f.read(1024*1024):pass
part.replace(archive)
print('RESULT_STRUCTURE_CHECK=PASS',flush=True)
print('RETURN='+str(archive),flush=True)
print('RETURN_BYTES='+str(archive.stat().st_size),flush=True)
print('RETURN_SHA256='+hashlib.sha256(archive.read_bytes()).hexdigest(),flush=True)
print('LUNG_MECHANISM_STEP7_PASS',flush=True)

