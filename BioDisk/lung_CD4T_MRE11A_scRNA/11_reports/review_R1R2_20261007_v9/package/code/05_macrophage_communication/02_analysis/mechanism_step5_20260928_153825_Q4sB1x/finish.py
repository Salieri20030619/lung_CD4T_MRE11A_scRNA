import os,json,tarfile,hashlib
from pathlib import Path
import pandas as pd
out=Path(os.environ['OUT']);root=Path(os.environ['ROOT'])
a=pd.read_csv(out/'tables/ligand_activities.tsv',sep='\t');s=pd.read_csv(out/'tables/geneset_status.tsv',sep='\t')
assert len(s)==64
cols=['receiver','filter','gate','test','direction','omit']
assert a[cols].drop_duplicates().shape[0]==(s.status=='PASS').sum()
assert not a.duplicated(cols+['test_ligand']).any()
assert (a['rank']>=1).all()
a[a.test_ligand.isin(['Il15','Tnf','Cd86','Tgfb1','Il1b','Tnfsf9','Il18','Cxcl16'])].sort_values(cols+['rank']).to_csv(out/'tables/focused_ligand_ranks.tsv',sep='\t',index=False)
key=['receiver','filter','gate','test','direction','test_ligand']
x=a[a.omit=='Full'][key+['rank','aupr_corrected']].merge(a[a.omit=='Drop7'][key+['rank','aupr_corrected']],on=key,how='outer',suffixes=('_Full','_Drop7'),validate='one_to_one')
x['top20_both']=(x.rank_Full<=20)&(x.rank_Drop7<=20)
x.to_csv(out/'tables/candidate_removal_sensitivity.tsv',sep='\t',index=False,na_rep='NA')
status=dict(status='PASS',stage='mechanism_step5_target_gate_sensitivity',output=str(out),previous=os.environ['PREV3'],sets_run=int((s.status=='PASS').sum()),sets_skipped=int((s.status!='PASS').sum()))
(out/'status.json').write_text(json.dumps(status,indent=2))
archive=root/'tmp/lung_mechanism_step5_return.tar.gz';part=archive.with_suffix('.tar.gz.part')
with tarfile.open(part,'w:gz') as t:
    for p in sorted(out.rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts:t.add(p,arcname=str(p.relative_to(out)))
with tarfile.open(part) as t:
    for m in t:
        if m.isfile():
            f=t.extractfile(m)
            while f.read(1024*1024):pass
part.replace(archive)
print('RETURN='+str(archive),flush=True)
print('RETURN_BYTES='+str(archive.stat().st_size),flush=True)
print('RETURN_SHA256='+hashlib.sha256(archive.read_bytes()).hexdigest(),flush=True)
print('LUNG_MECHANISM_STEP5_PASS',flush=True)

