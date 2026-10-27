from pathlib import Path
import json,re,subprocess
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/manuscript_graphics_20261006_v5'
for d in ['code','provenance','figures','tables','deliverables']:(O/d).mkdir(parents=True,exist_ok=True)
proc=subprocess.run(['ps','-eo','pid,etime,args'],capture_output=True,text=True).stdout.splitlines()
print('ACTIVE_PROJECT_PROCESSES', '\n'.join(x for x in proc if str(P) in x and 'graphics_inventory' not in x and 'ps -eo' not in x))
files=[];codes=[]
for f in P.rglob('*'):
 if not f.is_file() or '/manuscript_' in str(f) or '/reproduce/' in str(f):continue
 s=str(f.relative_to(P))
 if f.suffix.lower() in ['.rds','.rda','.rdata','.h5ad','.h5seurat'] or re.search('umap|embedding|rank|geneset|term2gene|gmt',f.name,re.I):
  files.append({'path':s,'bytes':f.stat().st_size})
 if f.suffix.lower() in ['.r','.py'] and ('pathway' in s.lower() or 'stage7_' in s or re.search('umap|annotat|seurat',f.name,re.I)):codes.append(s)
(O/'provenance/object_inventory.json').write_text(json.dumps(files,indent=2));(O/'provenance/code_inventory.json').write_text(json.dumps(codes,indent=2))
print('OBJECTS',json.dumps([x for x in files if not any(y in x['path'] for y in ['/SCENIC/','int/','cache/','regulon','bootstrap','prior'])]))
print('CODES',json.dumps(codes))
print('ENVIRONMENTS',[p.name for p in (P.parent.parent/'miniconda3/envs').glob('*')])
