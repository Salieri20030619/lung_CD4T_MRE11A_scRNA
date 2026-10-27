import json,hashlib,tarfile
from pathlib import Path
import numpy as np,pandas as pd
from PIL import Image
root=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
out=root/'02_analysis/scenic4_20261002_091545_iazXf1'
audit=root/'tmp/scenic4_repair_20261002_1205'
status=json.loads((out/'status.json').read_text())
assert status['status']=='PASS' and status['stage']=='lung_scenic4_external_transfer'
assert status['inferential_p'] is False and status['tf']==['IKZF2','STAT1']
names=['frozen_regulon_sample_activity.tsv','transfer_coverage.tsv','external_time_matched_TF_deltas.tsv','external_phenotype_sample.tsv','external_time_matched_phenotype_deltas.tsv','external_gene_detection.tsv','discovery_shared_universe_deltas.tsv','direction_replication.tsv']
tables={n:pd.read_csv(out/'tables'/n,sep='\t') for n in names}
assert all(len(d)>0 for d in tables.values())
e=tables['external_time_matched_TF_deltas.tsv']
assert set(e.seed)=={17,29,43} and set(e['filter'])=={'All','BG1'} and set(e.day)=={7,56}
assert set(e.feature)=={'IKZF2','STAT1'}
assert np.isfinite(e[['control','silica','delta']]).all().all()
assert np.allclose(e.delta,e.silica-e.control)
assert (e.control_cells>=20).all() and (e.silica_cells>=20).all()
assert set(e.evidence)=={'one_library_per_condition_time_no_inferential_P'}
a=tables['frozen_regulon_sample_activity.tsv']
external=a[a.cohort=='External']
assert external.groupby(['filter','scope','seed','day','group']).size().max()==1
coverage=tables['transfer_coverage.tsv']
assert len(coverage)==12 and (coverage.shared_genes>0).all()
obj=list((out/'objects').glob('*_shared_auc.csv.gz'))
assert len(obj)==6
object_shapes={}
for f in obj:
    z=pd.read_csv(f,index_col=0)
    assert z.index.is_unique and set(z.columns)=={'IKZF2','STAT1'} and np.isfinite(z).all().all()
    object_shapes[f.name]=list(z.shape)
for f in (out/'figures').glob('*.png'):
    with Image.open(f) as im: im.verify()
pdf=out/'figures/ExtTF.pdf'
assert pdf.read_bytes().startswith(b'%PDF-')
ret=root/'tmp/lung_scenic4_return.tar.gz'
members=[]
with tarfile.open(ret,'r:gz') as tf:
    for member in tf:
        if member.isfile():
            f=tf.extractfile(member)
            content=f.read()
            assert content==(out/member.name).read_bytes(),member.name
            members.append(member.name)
assert all('tables/'+n in members for n in names)
assert 'status.json' in members and 'figures/ExtTF.pdf' in members and 'figures/ExtTF.png' in members
v=tables['direction_replication.tsv']
v=v[v.scope=='ActCD4']
summary=v.groupby(['filter','day','feature']).agg(seeds=('seed','count'),delta_min=('delta','min'),delta_max=('delta','max'),discovery_min=('discovery_delta','min'),discovery_max=('discovery_delta','max'),matching_seed_directions=('direction_matches_discovery','sum')).reset_index()
result=dict(validation='PASS',row_counts={n:len(d) for n,d in tables.items()},objects=object_shapes,archive_bytes=ret.stat().st_size,archive_sha256=hashlib.sha256(ret.read_bytes()).hexdigest(),archive_member_count=len(members),actcd4=summary.to_dict('records'))
(audit/'verified_results.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
print('TRANSFER_COVERAGE')
print(coverage.to_csv(sep='\t',index=False))
print('EXTERNAL_SAMPLE_CELLS')
print(external[external.seed==17][['filter','scope','sample','group','day','cells']].to_csv(sep='\t',index=False))
print('ACTCD4_PHENOTYPE_DELTAS')
p=tables['external_time_matched_phenotype_deltas.tsv']
print(p[p.scope=='ActCD4'].to_csv(sep='\t',index=False))
print('POST_RUN_VERIFICATION=PASS')
