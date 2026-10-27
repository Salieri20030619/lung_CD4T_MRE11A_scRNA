from pathlib import Path
import json,gzip,csv,subprocess
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/GSE174725_reaudit_20261005_v1'
for n in ['code','results','provenance','deliverables']:(O/n).mkdir(parents=True,exist_ok=True)
spec='''Finite format/mapping audit, no disease model or clustering. Preserve all raw files. Reuse prior 5-CSV recovery excluding the single audited missing-value cell and skipping one blank row. Validate archive identity, matrix orientation/unit, raw barcode suffix mapping and cross-library barcode reuse. Compare every shared unambiguous literal gene (and separately deterministic R make.names aliases if needed) for every mapped CSV cell against pooled counts; record exact, <= and > differences, never rescale/impute to force agreement. Gene labels changed to dates or duplicate names remain ambiguous unless count-vector identity uniquely establishes mapping; no guessed symbol repair. Audit existing CD4 annotation lineage; quantify per-sample CD4 and provisional macrophage-marker-supported cells using fixed marker gates (T>=2 of CD3D/E/G/TRAC/TRBC1/2; myeloid>=2 of LYZ/LST1/TYROBP/FCER1G/AIF1/CTSS and macrophage>=1 of C1QA/B/C/MARCO/MRC1/CD68/CSF1R; exclude T and MS4A1/CD79A/EPCAM positive competitors for conservative count). Marker gate is feasibility, not final annotation; optional broad myeloid count separate. No full communication run. Large data remain server, no new download.\n'''
(O/'provenance/FROZEN_AUDIT_SPEC.txt').write_text(spec)
print(subprocess.run(['free','-h'],capture_output=True,text=True).stdout)
for name in ['02_analysis/screen2_TNK_20260911_103239_nfW0yu','02_analysis/screen2_Trefine_nuzVQu','02_analysis/screen2_CD4_target_wz1lfx']:
 d=P/name;print('CODE_DIR',name,[(x.name,x.stat().st_size) for x in d.iterdir() if x.is_file()])
 for p in d.glob('*.py'):
  t=p.read_text();print('CODE_RELEVANT',p)
  lines=t.splitlines()
  for i,l in enumerate(lines):
   if any(w in l for w in ['def load','GSE174725','csv','read_csv','make_unique','groupby','names_make','recover']):print('\n'.join(lines[max(0,i-2):i+7]))
for name in ['logs/02_corrected_filtered_qc_audit.log','logs/screen6_G174_audit.log','logs/screen6_G174_recover.log','logs/GSE174725_suppC2_resume.log']:
 p=P/name;(O/'provenance'/p.name).write_bytes(p.read_bytes())
for p in (P/'00_raw/GSE174725').iterdir():print('RAW',p.name,p.stat().st_size)
