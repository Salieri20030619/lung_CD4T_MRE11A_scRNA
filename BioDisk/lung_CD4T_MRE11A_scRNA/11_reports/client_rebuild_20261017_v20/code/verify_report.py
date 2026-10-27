from pathlib import Path
import json,subprocess,pymupdf as fitz
from PIL import Image,ImageDraw
O=Path(__file__).parent;d=json.loads((O/'report_delivery.json').read_text(encoding='utf-8'));p=Path(d['pdf']);doc=fitz.open(p);checks=[];images=0
for i,page in enumerate(doc):
 text=page.get_text();assert '\ufffd' not in text
 im=page.get_images(full=True);images+=len(im);bad=[]
 for b in page.get_text('blocks'):
  if b[0]<39 or b[2]>page.rect.width-37 or b[1]<24 or b[3]>page.rect.height-15:bad.append(b[:4])
 checks.append(dict(page=i+1,width=page.rect.width,height=page.rect.height,text_chars=len(text),images=len(im),out_of_bounds=bad))
assert images==22,images
R=O/'render/report';R.mkdir(parents=True,exist_ok=True)
subprocess.run([r'C:\Users\31937\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin\pdftoppm.exe','-scale-to','1500','-png',str(p),str(R/'page')],check=True,creationflags=subprocess.CREATE_NO_WINDOW)
files=sorted(R.glob('page-*.png'))
for j in range(0,len(files),8):
 sheet=Image.new('RGB',(1440,530*((min(8,len(files)-j)+1)//2)),'#dce4e8');draw=ImageDraw.Draw(sheet)
 for k,f in enumerate(files[j:j+8]):
  im=Image.open(f).convert('RGB');im.thumbnail((700,495));x=k%2*720;y=k//2*530;sheet.paste(im,(x,y+23));draw.text((x+8,y+5),f.stem,fill='black')
 sheet.save(R/('contact_'+str(j//8+1)+'.png'))
(O/'report_QA.json').write_text(json.dumps({'pages':len(doc),'figure_images':images,'checks':checks},ensure_ascii=False,indent=2),encoding='utf-8');print('PAGES',len(doc),'IMAGES',images,'BOUNDARY_FLAGS',sum(len(c['out_of_bounds']) for c in checks));print([(c['page'],c['text_chars'],c['images']) for c in checks])