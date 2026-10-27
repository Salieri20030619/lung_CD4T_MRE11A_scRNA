from pathlib import Path
import json,subprocess,fitz
from PIL import Image,ImageOps,ImageDraw
O=Path(__file__).parent;d=json.loads((O/'delivery.json').read_text(encoding='utf-8'));R=O/'render';R.mkdir(exist_ok=True)
subprocess.run([r'C:\Users\31937\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin\pdftoppm.exe','-scale-to','1150','-png',d['pdf'],str(R/'page')],check=True,creationflags=subprocess.CREATE_NO_WINDOW)
doc=fitz.open(d['pdf']);checks=[]
for i,p in enumerate(doc):
 txt=p.get_text();assert '\ufffd' not in txt and len(txt)>100
 bad=[b[:4] for b in p.get_text('blocks') if b[0]<35 or b[2]>560 or b[1]<30 or (b[3]>799 and b[1]<800)]
 checks.append({'page':i+1,'text_chars':len(txt),'out_of_body_bounds':bad})
imgs=[]
for p in sorted(R.glob('page-*.png')):
 im=Image.open(p).convert('RGB');im.thumbnail((398,565));imgs.append(im)
sheet=Image.new('RGB',(398*2,595*((len(imgs)+1)//2)),'#d6dce1');dr=ImageDraw.Draw(sheet)
for i,im in enumerate(imgs):x=(i%2)*398;y=(i//2)*595;sheet.paste(im,(x,y+20));dr.text((x+8,y+4),str(i+1),fill='black')
sheet.save(R/'contact.png');(O/'pdf_QA_v19.json').write_text(json.dumps({'pages':len(doc),'checks':checks},indent=2),encoding='utf-8');print('pages',len(doc),'bounds',sum(len(x['out_of_body_bounds']) for x in checks))