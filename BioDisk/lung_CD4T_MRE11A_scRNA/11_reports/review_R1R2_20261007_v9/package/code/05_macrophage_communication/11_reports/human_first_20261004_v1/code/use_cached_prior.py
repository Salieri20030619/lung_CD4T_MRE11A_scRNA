from pathlib import Path
import csv,hashlib,json,os,signal,datetime,shutil
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/human_first_20261004_v1')
C=Path('/mnt/BioDisk/wangmingyuan/IL17_Candida_integrated_v1/07_synthesis/supplementary_bioinformatics_v1/a6_consensus_regulatory_chain_v1/01_nichenet_resources_v1')
rows=[]
for r in csv.DictReader((C/'source_manifest_v1.tsv').open(),delimiter='\t'):
 f=C/r['file'];h=hashlib.sha256(f.read_bytes()).hexdigest();assert h==r['sha256'] and f.stat().st_size==int(r['bytes']);rows.append(dict(path=str(f),species='human',version='NicheNet-v1',source=r['url'],sha256=h,bytes=f.stat().st_size))
for p in Path('/proc').iterdir():
 if not p.name.isdigit():continue
 try:
  args=(p/'cmdline').read_bytes().split(b'\0')
  if args and args[0]==b'curl' and str(O/'resources/ligand_target_matrix_nsga2r_final.rds.part').encode() in args:os.kill(int(p.name),signal.SIGTERM);print('STOPPED_OWN_V2_DOWNLOAD',p.name)
 except (FileNotFoundError,PermissionError,ProcessLookupError):pass
(O/'provenance/human_prior_manifest.json').write_text(json.dumps(rows,indent=2))
(O/'provenance/prior_choice_before_analysis.txt').write_text('Before human DE/NicheNet computation: reuse verified locally installed official HUMAN NicheNet v1 (Zenodo3260758) LR+LT pair. Attempted v2 download was too slow and stopped; no v2 results computed. Human v1 and existing mouse v2 are different priors; cross-species rank agreement is contextual, not direct comparable numeric validation. No endpoint or expression gate changed. '+datetime.datetime.now().isoformat())
print('CACHED_HUMAN_V1_VERIFIED',rows)
