from pathlib import Path
import subprocess,json,hashlib
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/GSE296682_patient_localization_20261005_v1');p=O/'public_cache/BAL_FINAL.h5ad.range0_65535';url='https://zenodo.org/records/18626573/files/BAL_FINAL.h5ad?download=1'
r=subprocess.run(['curl','--fail','-L','--connect-timeout','15','--max-time','45','--max-filesize','65536','--silent','--show-error','--range','0-65535','--output',str(p),'--write-out','%{http_code}',url],capture_output=True,text=True,timeout=50)
x=dict(returncode=r.returncode,http=r.stdout,stderr=r.stderr,bytes=p.stat().st_size if p.exists() else 0,sha256=hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None);(O/'provenance/range_probe.json').write_text(json.dumps(x,indent=2));print(json.dumps(x))
