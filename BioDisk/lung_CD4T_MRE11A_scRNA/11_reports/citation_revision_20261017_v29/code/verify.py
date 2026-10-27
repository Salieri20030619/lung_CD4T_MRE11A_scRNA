from pathlib import Path
import json,fitz,re,subprocess
from PIL import Image,ImageDraw
O=Path(r'C:\Users\31937\lung_citation_revision_20261017_v29');p=Path(json.loads((O/'report_delivery.json').read_text(encoding='utf-8'))['pdf']);d=fitz.open(p);bad=[];figs=[]
for i,pg in enumerate(d):
 for b in pg.get_text('blocks'):
  if b[0]<39 or b[2]>pg.rect.width-37 or b[1]<24 or b[3]>pg.rect.height-15:bad.append([i+1,list(b[:4]),b[4][:80]])
 for n in re.findall(r'(?m)^图(\d+)。',pg.get_text()):figs.append([int(n),i+1])
print('PAGES',len(d),'BOUNDS',bad,'VISUALS',sum(len(pg.get_images()) for pg in d))
for i,pg in enumerate(d):
 if i>=33:print(i+1,len(pg.get_text()),pg.get_text().splitlines()[2:3])
t='\n'.join(pg.get_text() for pg in d);assert '待排' not in t
print('REF10_PAGES',[i+1 for i,pg in enumerate(d) if '[10]' in pg.get_text()])
r=O/'render';r.mkdir(exist_ok=True);exe=r'C:\Users\31937\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin\pdftoppm.exe'
nums=[2,29,32,34]+list(range(37,len(d)+1))
for n in nums:subprocess.run([exe,'-f',str(n),'-l',str(n),'-scale-to','1500','-png',str(p),str(r/'page')],check=True,creationflags=subprocess.CREATE_NO_WINDOW)
files=sorted(r.glob('page-*.png'))
for start in range(0,len(files),6):
 sheet=Image.new('RGB',(1500,1450),'#ddd');dr=ImageDraw.Draw(sheet)
 for j,f in enumerate(files[start:start+6]):
  im=Image.open(f).convert('RGB');im.thumbnail((490,690));x=(j%3)*500;y=(j//3)*725;sheet.paste(im,(x,y+25));dr.text((x+8,y+5),f.stem,fill='black')
 sheet.save(r/f'contact_{start//6+1:02d}.png')
(O/'layout_QA.json').write_text(json.dumps({'pages':len(d),'bounds':bad,'figures':figs,'visuals':sum(len(pg.get_images()) for pg in d)},ensure_ascii=False,indent=2),encoding='utf-8');(O/'extracted.txt').write_text(t,encoding='utf-8')