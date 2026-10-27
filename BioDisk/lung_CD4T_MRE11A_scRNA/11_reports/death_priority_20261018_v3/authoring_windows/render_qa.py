from pathlib import Path
import json,hashlib,subprocess
import pymupdf as fitz
from PIL import Image,ImageDraw
P=Path(r'C:\Users\31937\lung_death_priority_20261018_v3');info=json.loads((P/'report_delivery.json').read_text(encoding='utf-8'));pdf=Path(info['pdf']);q=P/'pdf_qa';q.mkdir(exist_ok=True)
doc=fitz.open(pdf);issues=[]
for i,page in enumerate(doc):
 for b in page.get_text('blocks'):
  if b[0]<39 or b[2]>page.rect.width-38 or b[1]<23 or b[3]>page.rect.height-14:issues.append([i+1,list(b[:4]),b[4][:60]])
text='\n'.join(p.get_text() for p in doc)
assert len(doc)==info['pages']
assert all(p.rect.height>p.rect.width for p in doc)
assert '旧网络' not in text and '旧投射' not in text
assert '\ufffd' not in text
assert all(g in text for g in ['HIF1A','IL2RA','ARL6IP1','0.003025','0.073858','11,559','41,604'])
(P/'report_text.txt').write_text(text,encoding='utf-8')
pop=Path(r'C:\Users\31937\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin\pdftoppm.exe')
subprocess.run([str(pop),'-r','110','-png',str(pdf),str(q/'page')],check=True,creationflags=subprocess.CREATE_NO_WINDOW)
files=sorted(q.glob('page-*.png'))
for start in range(0,len(files),6):
 c=Image.new('RGB',(1800,1500),'#dbe1e4');d=ImageDraw.Draw(c)
 for k,f in enumerate(files[start:start+6]):
  im=Image.open(f);im.thumbnail((590,705));xx=k%3*600;yy=k//3*750;c.paste(im,(xx+(590-im.width)//2,yy+25));d.text((xx+8,yy+5),f.stem,fill='black')
 c.save(q/f'contact_{start:02d}.jpg',quality=92)
out={'pages':len(doc),'boundary_issues':issues,'sha256':hashlib.sha256(pdf.read_bytes()).hexdigest(),'pdf_bytes':pdf.stat().st_size,'text_no_replacement_character':True}
(P/'pdf_QA.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(out,ensure_ascii=False))
