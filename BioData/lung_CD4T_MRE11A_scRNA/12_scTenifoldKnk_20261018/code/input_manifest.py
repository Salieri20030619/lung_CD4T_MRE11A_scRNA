from pathlib import Path
import hashlib,json
O=Path('/mnt/BioData/bio/wangmingyuan/lung_CD4T_MRE11A_scRNA/12_scTenifoldKnk_20261018');P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/three_cohort_four_axis_20261017_v16/inputs')
rows=[]
for ds in ['GSE314943','GSE183682']:
 for suf in ['.cellcounts.npz','.genes.txt','.cells.tsv.gz']:
  f=P/(ds+'__Core'+suf);h=hashlib.sha256()
  with f.open('rb') as r:
   for b in iter(lambda:r.read(8*1024*1024),b''):h.update(b)
  rows.append(dict(path=str(f),bytes=f.stat().st_size,sha256=h.hexdigest()))
for f in sorted((O/'inputs').glob('*/*')):
 h=hashlib.sha256(f.read_bytes()).hexdigest();rows.append(dict(path=str(f),bytes=f.stat().st_size,sha256=h))
(O/'input_manifest.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
print('INPUT_MANIFEST',len(rows))