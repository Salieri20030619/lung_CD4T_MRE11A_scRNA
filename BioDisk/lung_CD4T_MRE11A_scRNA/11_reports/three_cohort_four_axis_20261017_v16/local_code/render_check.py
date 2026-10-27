from pathlib import Path
import json,subprocess,fitz
from PIL import Image,ImageOps,ImageDraw
O=Path(r'C:\Users\31937\lung_three_cohort_20261017_v16');j=json.loads((O/'report_output.json').read_text(encoding='utf-8'));pdf=Path(j['pdf']);tmp=O/'render_final';tmp.mkdir(exist_ok=True)
doc=fitz.open(pdf);print('PAGES',len(doc));audit=[]
for i,p in enumerate(doc):
 txt=p.get_text();lines=txt.splitlines();print(i+1,lines[0][:90],len(txt));assert '\ufffd' not in txt;assert len(txt)>100
 for b in p.get_text('blocks'):
  if b[0]<35 or b[2]>p.rect.width-35 or b[1]<25 or b[3]>p.rect.height-15:audit.append((i+1,b[:4],b[4][:40]))
print('OUT_OF_BOUNDS',audit)
exe=Path(r'C:\Users\31937\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin\pdftoppm.exe')
r=subprocess.run([str(exe),'-r','90','-png',str(pdf),str(tmp/'page')],capture_output=True,creationflags=subprocess.CREATE_NO_WINDOW);r.check_returncode()
files=sorted(tmp.glob('page-*.png'))
for k in range(0,len(files),6):
 canvas=Image.new('RGB',(1000,2250),'#dce4e8');draw=ImageDraw.Draw(canvas)
 for n,f in enumerate(files[k:k+6]):
  im=Image.open(f);im.thumbnail((480,700));x=(n%2)*500+10;y=(n//2)*750+30;canvas.paste(im,(x,y));draw.text((x,y-20),f'Page {k+n+1}',fill='black')
 canvas.save(tmp/f'contact_{k//6+1}.png')
(O/'pdf_qa.json').write_text(json.dumps({'pages':len(doc),'out_of_bounds':audit,'rendered_pages':len(files)},ensure_ascii=False,indent=2),encoding='utf-8')