from pathlib import Path
import os,json,hashlib
os.environ['OPENBLAS_NUM_THREADS']='1'
import pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');B=P/'11_reports';O=B/'independent_candidates_revision_20261005_v1'
for s in ['code','results','provenance','figures','deliverables']:(O/s).mkdir(parents=True,exist_ok=True)
pool='IL15 TNF IKZF2 STAT1 CLSTN1 ADAM15 PLAU DSC2 CXCL2 ICAM2 VTN VEGFB GSTP1 RTN4 IL18'.split()
(O/'provenance/REVISION_SCOPE.txt').write_text('User correction: each fixed candidate is an independent proposal. No requirement for IL15/TNF -> IKZF2/STAT1 or IKZF2 -> STAT1 cascade. Withdraw feasibility downgrades based solely on absent cross-candidate links. Preserve valid within-gene disease, phenotype/program, TF activity, communication prediction, donor association and independent localization results at their own evidence levels. Fixed15: '+','.join(pool)+'. Phase1 reuses existing tables without rerunning tests. Distinguish nominal/FDR, prediction ranks/no P, unmeasured/insufficient/conditional effects. Roles need not share inapplicable tests. Macrophage annotation, clinical-source confounding, genuine specificity and sample limitations remain. New limited gap analysis only after tracing existing evidence; no significance-driven grids. Only report/necessary figures local+Library; matrices, standalone donor/cell tables and ZIP stay server.\n')
inventory=[]
for folder in B.iterdir():
 if not folder.is_dir() or folder==O:continue
 for p in folder.rglob('*'):
  if p.is_file() and p.name.endswith(('.tsv','.tsv.gz','.csv','.csv.gz','.json','.txt')) and not any(x in p.parts for x in ['public_cache','dependencies']):
   inventory.append(dict(report=folder.name,path=str(p),bytes=p.stat().st_size))
pd.DataFrame(inventory).to_csv(O/'provenance/existing_evidence_inventory.tsv',sep='\t',index=False)
for report in ['technical_delivery_20261002_ikzf2_stat1_v1','stepwise_diagnostic_20261002_v1','missing_links_20261003_v1','mouse_sender_20261002_v1','dual_line_TNF_IL18_20261004_v1','breakpoint_TNF_IL18_20261004_v1','bootstrap_TNF_IL18_20261004_v1']:
 files=[x for x in inventory if x['report']==report and x['path'].endswith(('.tsv','.tsv.gz')) and not any(t in Path(x['path']).name.lower() for t in ['cell','member','edge','target_weight','sample_scores','donor_scores']) and (report!='technical_delivery_20261002_ikzf2_stat1_v1' or '/revision_v2/final_release/data/' in x['path'] or '/revision_v2/final_release/source_tables/' in x['path'])]
 print('REPORT',report,'TABLES',len(files),flush=True)
 for x in files:
  if x['bytes']<40_000_000:print(x['path'].replace(str(B)+'/',''),x['bytes'],flush=True)
print('INVENTORY_DONE',len(inventory),flush=True)
