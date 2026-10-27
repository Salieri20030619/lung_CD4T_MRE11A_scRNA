from pathlib import Path
import requests,json,hashlib,subprocess,concurrent.futures
O=Path(r'C:\Users\31937\lung_three_cohort_20261017_v16');(O/'resources').mkdir(exist_ok=True)
urls={'CellChatDB.mouse.rda':'https://raw.githubusercontent.com/jinworks/CellChat/main/data/CellChatDB.mouse.rda'}
for p in ['R-MMU-109606','R-MMU-5620971']:urls[p+'.json']='https://reactome.org/ContentService/data/participants/'+p+'/referenceEntities'
def fetch(kv):
 n,u=kv
 try:
  r=requests.get(u,timeout=(15,60));r.raise_for_status();(O/'resources'/n).write_bytes(r.content);return dict(file=n,url=u,sha256=hashlib.sha256(r.content).hexdigest(),bytes=len(r.content))
 except Exception as e:return dict(file=n,url=u,error=str(e))
with concurrent.futures.ThreadPoolExecutor(3) as ex:r=list(ex.map(fetch,urls.items()))
(O/'resources/download_manifest.json').write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps(r,indent=2))
for f in (O/'resources').iterdir():
 subprocess.run(['scp','-q',str(f),'sgs10:/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/three_cohort_four_axis_20261017_v16/resources/'],check=True,creationflags=subprocess.CREATE_NO_WINDOW)