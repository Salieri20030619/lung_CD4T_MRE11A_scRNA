from pathlib import Path
import pandas as pd,json,zipfile,hashlib,shutil
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/review_revision_20261017_v14';Q=O/'supplement_v14';BASE=P/'11_reports/methods_complete_20261016_v13'
# Avoid re-labelling existing results using the smaller new sensitivity family.
f=Q/'S27_control_overlap_associations.tsv';a=pd.read_csv(f,sep='\t');a.loc[a.variant=='original_score','BH']=float('nan');a.loc[a.variant=='original_score','family']='original_reference_only_see_original_full_family';a['BH_note']='BH only for new clean-control sensitivity; original P/rho reference, original full-family BH in S08/S21';a.to_csv(f,sep='\t',index=False,na_rep='NA')
f=P/'11_reports/integrated_macrophage_CD4_20261015_v4/06_support/Core_fixed_composition_TF_program.tsv';shutil.copy2(f,Q/'S21_complete_TF_fixed.tsv')
manifest={str(f.relative_to(BASE)):hashlib.sha256(f.read_bytes()).hexdigest() for f in BASE.rglob('*') if f.is_file()};(O/'provenance/baseline_remote_sha256.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
archive=O/'review_analysis_bundle.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
 for d in ['supplement_v14','provenance']:
  for f in (O/d).rglob('*'):
   if f.is_file():z.write(f,f.relative_to(O))
print('BUNDLE',archive.stat().st_size,hashlib.sha256(archive.read_bytes()).hexdigest())