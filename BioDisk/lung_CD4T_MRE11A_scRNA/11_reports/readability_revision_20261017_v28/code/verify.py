from pathlib import Path
import fitz,json,re,hashlib,subprocess
from PIL import Image,ImageDraw
O=Path(r'C:\Users\31937\lung_readability_revision_20261017_v28');p=Path(json.loads((O/'report_delivery.json').read_text(encoding='utf-8'))['pdf']);d=fitz.open(p);r=O/'render';r.mkdir(exist_ok=True)
text='\n'.join(pg.get_text() for pg in d);(O/'extracted.txt').write_text(text,encoding='utf-8');bad=[]
for i,pg in enumerate(d):
 for b in pg.get_text('blocks'):
  if b[0]<39 or b[2]>pg.rect.width-37 or b[1]<24 or b[3]>pg.rect.height-15:bad.append([i+1,list(b[:4]),b[4][:100]])
figs=[]
for i,pg in enumerate(d):
 for f in re.findall(r'(?m)^图(\d+)。',pg.get_text()):figs.append([int(f),i+1])
assert [x[0] for x in figs]==list(range(1,25)),figs
assert '待排' not in text
assert all(pg.rect.height>pg.rect.width for pg in d)
print('PAGES',len(d),'BOUNDS',bad,'FIGURES',figs)
for n in [1,2,3,7,10,25,28,31,32,33,35,47,48,49,50]:
 if n<=len(d):print('PAGE',n,d[n-1].get_text().splitlines()[2])
exe=r'C:\Users\31937\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin\pdftoppm.exe'
subprocess.run([exe,'-scale-to','1500','-png',str(p),str(r/'page')],check=True,creationflags=subprocess.CREATE_NO_WINDOW)
files=sorted(r.glob('page-*.png'))
for start in range(0,len(files),6):
 sheet=Image.new('RGB',(1500,1450),'#ddd');dr=ImageDraw.Draw(sheet)
 for j,f in enumerate(files[start:start+6]):
  im=Image.open(f).convert('RGB');im.thumbnail((490,690));x=(j%3)*500;y=(j//3)*725;sheet.paste(im,(x,y+25));dr.text((x+8,y+5),f.stem,fill='black')
 sheet.save(r/f'contact_{start//6+1:02d}.png')
(O/'layout_QA.json').write_text(json.dumps({'pages':len(d),'bounds':bad,'figures':figs,'portrait_all':True,'visual_count':sum(len(pg.get_images()) for pg in d),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()},ensure_ascii=False,indent=2),encoding='utf-8')