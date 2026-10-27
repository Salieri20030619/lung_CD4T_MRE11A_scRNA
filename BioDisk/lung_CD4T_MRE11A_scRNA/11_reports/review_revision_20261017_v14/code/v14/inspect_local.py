from pathlib import Path
from bs4 import BeautifulSoup
import hashlib,json
p=Path('C:/Users/31937/lung_review_20261017_v14/baseline_v13');s=BeautifulSoup((p/'report_v13.html').read_text(encoding='utf-8'),'html.parser')
for i,sec in enumerate(s.find_all('section')):
 print(i,sec.find(['h1','h2']).get_text())
 if i==0 or any(sec.find(['h1','h2']).get_text().startswith(x) for x in ['4.4','5.2','7 ','8 ','9 ','参考']):print(sec.get_text(' ',strip=True))
print('PACKAGES')
import reportlab,fitz,pandas,scipy
print('OK')