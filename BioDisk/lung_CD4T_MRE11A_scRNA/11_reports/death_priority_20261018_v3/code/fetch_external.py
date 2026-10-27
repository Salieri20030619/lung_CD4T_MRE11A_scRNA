from pathlib import Path
import requests,gzip,re,json,hashlib,concurrent.futures
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/death_priority_20261018_v3/external')
P.mkdir(parents=True,exist_ok=True)
def fetch(acc):
 d=P/acc;d.mkdir(exist_ok=True)
 base='https://ftp.ncbi.nlm.nih.gov/geo/series/'+acc[:-3]+'nnn/'+acc
 for suffix,name in [('/soft/'+acc+'_family.soft.gz','family.soft.gz'),('/suppl/','supplement_listing.html')]:
  p=d/name
  if not p.exists():
   r=requests.get(base+suffix,timeout=120);r.raise_for_status();p.write_bytes(r.content)
 txt=gzip.decompress((d/'family.soft.gz').read_bytes()).decode('utf-8')
 (d/'family.soft').write_text(txt,encoding='utf-8')
 print(acc, '\n'.join(l for l in txt.splitlines() if l.startswith(('!Series_title','!Series_overall_design','!Series_supplementary_file','!Sample_title','!Sample_characteristics_ch1','!Sample_description'))),flush=True)
for acc in ['GSE121883','GSE253090','GSE324038','GSE159354','GSE184854','GSE295463','GSE296682']:
 try:fetch(acc)
 except Exception as e:print(acc,'ERROR',repr(e),flush=True)
