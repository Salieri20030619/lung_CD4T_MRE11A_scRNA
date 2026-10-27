from pathlib import Path
import json,hashlib,base64,io
from bs4 import BeautifulSoup
from PIL import Image
D=Path(__file__).parent/'g174_annotation_delivery';p=D/'GSE174725_annotation_validation_20261005_report_zh.html';s=BeautifulSoup(p.read_text(encoding='utf-8'),'html.parser');text=s.get_text()
assert '\ufffd' not in text and '444' in text and '180' in text and '肺泡巨噬' in text
assert not s.find_all('script');assert len(s.find_all('img'))==2
for im in s.find_all('img'):Image.open(io.BytesIO(base64.b64decode(im['src'].split(',')[1]))).verify()
assert all(len(t.find_all('tr'))>1 for t in s.find_all('table'))
assert p.stat().st_size>500000
q=BeautifulSoup(str(s),'html.parser');h=q.find('h2')
for x in list(h.find_next_siblings()):x.decompose()
h.decompose();(D/'qa_header.html').write_text(str(q),encoding='utf-8');out=dict(status='PASS',sha256=hashlib.sha256(p.read_bytes()).hexdigest(),tables=len(s.find_all('table')),embedded_figures=2,external_assets=False,visual_QA='pending');(D/'QA_local.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))
