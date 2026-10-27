import os,json,hashlib
from pathlib import Path
for var,stage in [('PREV1','mechanism_step1'),('PREV3','mechanism_step3_nichenet'),('PREV4','mechanism_step4_receptor_audit')]:
 p=Path(os.environ[var]);s=json.loads((p/'status.json').read_text());assert s['status']=='PASS' and s['stage']==stage
p3=Path(os.environ['PREV3']);cache=Path(os.environ['CACHE'])
manifest=json.loads((p3/'resource_manifest.json').read_text())
record=next(v for v in manifest if v['file']=='ligand_target_matrix_nsga2r_final_mouse.rds')
p=cache/record['file'];assert p.is_file() and p.stat().st_size==record['bytes']
algorithm,expected=record['checksum'].split(':',1);h=hashlib.new(algorithm)
with p.open('rb') as f:
 for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
assert h.hexdigest()==expected,'NicheNet cached prior checksum mismatch'
print('INPUT_AND_PRIOR_CHECK=PASS',flush=True)

