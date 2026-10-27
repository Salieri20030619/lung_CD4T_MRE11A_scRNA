import os,json,hashlib,tarfile
from pathlib import Path
import pandas as pd
import numpy as np
out=Path(os.environ['OUT']);root=Path(os.environ['ROOT'])
def read(n):return pd.read_csv(out/'tables'/f'{n}.tsv',sep='\t')
a=read('TF_sample_activity');e=read('TF_group_effects');r=read('TF_associations');c=read('TF_coverage')
assert not a.duplicated(['key','level','mode','tf','sample']).any()
assert np.isfinite(a.activity).all() and (a.n_targets>=10).all()
assert (a.groupby(['key','level','mode','tf']).size()==6).all()
assert {'CD4_Core_All','CD4_Core_BG1','CD4_ActCD4_All','CD4_ActCD4_BG1'} <= set(a.key)
assert (r.rho.dropna().abs()<=1+1e-10).all()
assert set(r[r.kind=='Module'].exclusion)=={'Drop7PlusModule'}
assert set(r[r.kind=='Candidate'].exclusion)=={'Drop7'}
# Wide, complete evidence matrix; no favorable TF or module is removed.
z=e[e['mode']=='Drop7'].copy()
z['effect_LOO_same_direction']=((z.delta>0)&(z.LOO_min>0))|((z.delta<0)&(z.LOO_max<0))
for kind,features in [('Candidate',['IL1R2','PENK','TTC39C','FOXP1','TNFRSF9','LAG3','MRE11']),('AMlike',['IL15','TNF']),('Module',['MitoISR','OXPHOS','MitoROS','Mitophagy','PyroExec','TAct','TReg','Inhib'])]:
 for feature in features:
  t=r[(r.kind==kind)&(r.feature==feature)][['key','level','tf','rho','LOO_min','LOO_max','LOO_n']].copy()
  prefix=f'{kind}_{feature}_';t=t.rename(columns={v:prefix+v for v in ['rho','LOO_min','LOO_max','LOO_n']})
  z=z.merge(t,on=['key','level','tf'],how='left',validate='one_to_one')
z['abs_delta']=z.delta.abs();z=z.sort_values(['key','level','abs_delta'],ascending=[True,True,False]);z.to_csv(out/'tables/TF_evidence_matrix.tsv',sep='\t',index=False,na_rep='NA')
e1=e[e['mode']=='Full'];e2=e[e['mode']=='Drop7'];b=e1.merge(e2,on=['key','scope','filter','level','tf'],suffixes=('_Full','_Drop7'),validate='one_to_one')
b['direction_agreement']=np.sign(b.delta_Full)==np.sign(b.delta_Drop7)
b.to_csv(out/'tables/candidate_exclusion_sensitivity.tsv',sep='\t',index=False)
files=['TF_sample_activity','TF_group_effects','TF_associations','TF_target_edges','candidate_prior_coverage','TF_evidence_matrix']
for n in files:assert (out/'tables'/f'{n}.tsv').is_file()
status=dict(status='PASS',stage='mechanism_step6_TF_network',output=str(out),previous=os.environ['PREV3'],scope_n=int(a.key.nunique()),TF_n_ABC=int(a[a.level=='ABC'].tf.nunique()),activity_rows=len(a),pending=['human_macrophage_support','human_TF_program_replication','final_mechanism_literature_review'])
(out/'status.json').write_text(json.dumps(status,indent=2))
archive=root/'tmp/lung_mechanism_step6_return.tar.gz';part=archive.with_suffix('.tar.gz.part')
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
print('LUNG_MECHANISM_STEP6_PASS',flush=True)

