from pathlib import Path
import json,subprocess
O=Path(__file__).parent;P='/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA'
r=subprocess.run(['ssh','-T','sgs10','find '+P+'/11_reports '+P+'/10_tables '+P+'/06_CD4T_analysis '+P+'/09_figures -type f'],capture_output=True,check=True,creationflags=subprocess.CREATE_NO_WINDOW)
paths=r.stdout.decode('utf-8').splitlines();(O/'server_file_inventory.json').write_text(json.dumps(paths,indent=2),encoding='utf-8')
for p in paths:
 if any(s in p.lower() for s in ['front_half','trajectory','pseudotime','slingshot','monocle','gsea','enrich','umap','composition']) and any(p.endswith(s) for s in ['.tsv','.csv','.json','.py','.R','.h5ad']):print(p)