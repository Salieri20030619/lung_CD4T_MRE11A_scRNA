from pathlib import Path
import json,subprocess,hashlib,datetime
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/human_first_20261004_v1';cache=O/'resources';cache.mkdir(exist_ok=True)
j=json.loads((P/'00_raw/nichenet_v2_mouse/zenodo_7074291.json').read_text());records={x['key']:x for x in j['files']};rows=[]
for name in ['lr_network_human_21122021.rds','ligand_target_matrix_nsga2r_final.rds']:
 r=records[name];f=cache/name
 if not f.exists():
  print('FETCH',name,flush=True)
  subprocess.run(['curl','-L','--fail','--silent','--show-error','--connect-timeout','20','--max-time','900','--retry','1','-o',str(f)+'.part','https://zenodo.org/records/7074291/files/'+name+'?download=1'],check=True)
  f.with_suffix('.rds.part').rename(f)
 h=hashlib.md5(f.read_bytes()).hexdigest();assert 'md5:'+h==r['checksum'] and f.stat().st_size==r['size']
 rows.append(dict(file=str(f),species='human',source='https://zenodo.org/records/7074291',version='NicheNet-v2',md5=h,bytes=f.stat().st_size));print('VERIFIED',name,flush=True)
(O/'provenance/human_prior_manifest.json').write_text(json.dumps(rows,indent=2));print('HUMAN_PRIOR_READY',flush=True)
