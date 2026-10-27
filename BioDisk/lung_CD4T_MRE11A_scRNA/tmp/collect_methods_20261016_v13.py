from pathlib import Path
import tarfile,json,hashlib
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
paths=[
'02_analysis/GSE314943/stage2_20260905_190530_BeH7tA',
'02_analysis/GSE314943/stage3_20260905_214228_VPtd0r',
'02_analysis/GSE314943/stage4_20260905_215942_3XPis4',
'02_analysis/GSE314943/stage5_20260905_222457_a2dApT',
'02_analysis/GSE314943/stage6_20260905_231807_L44XdT',
'02_analysis/GSE314943/stage9_20260906_080538_wN3sjO',
'02_analysis/GSE314943/stage11_20260906_100314_pALwAB',
'02_analysis/scenic1_20260929_103513_xRoEnr',
'02_analysis/scenic2_20261001_092655_ZKNWD6',
'11_reports/upstream_human_support_20261015_v6/code',
'11_reports/STAT1_pyro_bridge_20261015_v5/code',
'11_reports/customer_sequence_20261015_v2/code',
'11_reports/GSE174725_annotation_validation_20261005_v1',
'11_reports/GSE296682_patient_localization_20261005_v1']
files=[]
for rel in paths:
 d=P/rel
 if not d.exists():continue
 for sub in [d,d/'code',d/'scripts',d/'provenance']:
  if sub.exists():
   for f in sub.iterdir():
    if f.is_file() and f.suffix in ['.py','.R','.sh','.txt','.json'] and f.stat().st_size<2000000:files.append(f)
files=sorted(set(files));a=P/'tmp/methods_sources_20261016_v13.tar.gz'
if a.exists():raise RuntimeError('Archive already exists')
with tarfile.open(a,'w:gz') as tar:
 for f in files:tar.add(f,arcname=str(f.relative_to(P)))
print(json.dumps({'files':len(files),'archive':str(a),'sha256':hashlib.sha256(a.read_bytes()).hexdigest()}))
