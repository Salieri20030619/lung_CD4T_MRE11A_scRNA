import os,json,hashlib,tarfile
from pathlib import Path
import pandas as pd
out=Path(os.environ['OUT']);root=Path(os.environ['ROOT'])
a=pd.read_csv(out/'tables/annotation_sample_expression.tsv',sep='\t');m=pd.read_csv(out/'tables/sample_mapping_status.tsv',sep='\t');coverage=pd.read_csv(out/'tables/panel_coverage.tsv',sep='\t')
assert not a.duplicated(['column','label','sample','gene']).any()
assert a.pct.dropna().between(0,100).all()
assert a.IL15_IL15RA_pct.dropna().between(0,100).all()
assert (m.mapped_cells+m.unknown_cells==m.total_cells).all()
status=dict(status='PASS',stage='mechanism_step8_human_mac_input_audit',output=str(out),input=os.environ['HUMAN_OBJECT'],annotation_columns=int(a['column'].nunique()),mapped_subjects=int(a['sample'].nunique()),total_cells=int(m.total_cells.iloc[0]),mapped_cells=int(m.mapped_cells.iloc[0]),panel_genes_present=int(coverage.present.sum()),next='review_macrophage_identity_and_matched_subjects_before_sender_statistics')
(out/'status.json').write_text(json.dumps(status,indent=2))
archive=root/'tmp/lung_mechanism_step8_return.tar.gz';part=archive.with_suffix('.tar.gz.part')
with tarfile.open(part,'w:gz') as t:
 for p in sorted(out.rglob('*')):
  if p.is_file() and not {'objects','__pycache__'}&set(p.relative_to(out).parts):t.add(p,arcname=str(p.relative_to(out)))
with tarfile.open(part) as t:
 for m in t:
  if m.isfile():
   f=t.extractfile(m)
   while f.read(1024*1024):pass
part.replace(archive)
print('RETURN='+str(archive),flush=True)
print('RETURN_BYTES='+str(archive.stat().st_size),flush=True)
print('RETURN_SHA256='+hashlib.sha256(archive.read_bytes()).hexdigest(),flush=True)
print('LUNG_MECHANISM_STEP8_PASS',flush=True)

