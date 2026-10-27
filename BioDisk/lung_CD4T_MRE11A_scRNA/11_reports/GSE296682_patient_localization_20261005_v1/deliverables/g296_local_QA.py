from pathlib import Path
import json,hashlib,base64,io
from bs4 import BeautifulSoup
from PIL import Image
D=Path(__file__).parent/'g296_delivery';p=D/'GSE296682_patient_localization_20261005_report_zh.html';s=BeautifulSoup(p.read_text(encoding='utf-8'),'html.parser');t=s.get_text();assert '\ufffd' not in t;assert '348004' in t and 'IL15RA' in t and '没有计算P/FDR' in t;assert not s.find_all('script');assert len(s.find_all('img'))==2
for im,name in zip(s.find_all('img'),['01_author_label_marker_review.png','02_fixed_candidate_localization.png']):
 data=base64.b64decode(im['src'].split(',')[1]);Image.open(io.BytesIO(data)).verify();(D/name).write_bytes(data)
assert len(s.find_all('table'))>=12;assert all(len(x.find_all('tr'))>1 for x in s.find_all('table'))
q=BeautifulSoup(str(s),'html.parser');h=q.find_all('h2')[1]
for x in list(h.find_next_siblings()):x.decompose()
h.decompose();(D/'qa_header.html').write_text(str(q),encoding='utf-8');out=dict(status='PASS',report_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),report_bytes=p.stat().st_size,tables=len(s.find_all('table')),embedded_figures=2,external_assets=False,visual_QA='pending');(D/'QA_local.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))
