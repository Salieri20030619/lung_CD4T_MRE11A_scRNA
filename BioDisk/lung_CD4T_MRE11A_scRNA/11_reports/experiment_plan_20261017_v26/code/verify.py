from pathlib import Path
import json, hashlib, fitz, subprocess
from PIL import Image, ImageOps, ImageDraw
w=Path(r'C:\Users\31937\lung_experiment_plan_20261017_v26')
m=json.loads((w/'report_delivery.json').read_text(encoding='utf-8')); p=Path(m['pdf'])
b=Path(r'C:\Users\31937\Downloads\矽肺_CD4候选研究_队列衔接修订_v25_20261017\矽肺免疫全景_CD4候选与验证_报告_v25.pdf')
d=fitz.open(p); old=fitz.open(b)
assert len(d)==42
assert all(d[i].get_text()==old[i].get_text() for i in range(34))
assert sum(len(x.get_images()) for x in d)==27
assert all(x.rect.height>x.rect.width for x in d)
assert hashlib.sha256(b.read_bytes()).hexdigest()=='6d1d715a9873f531c250a4c477bda496066834f333eecec62e8afb5ff8549c34'
bad=[]
for i,x in enumerate(d):
 for z in x.get_text('blocks'):
  if z[0]<39 or z[2]>x.rect.width-37 or z[1]<24 or z[3]>x.rect.height-15: bad.append([i+1,list(z[:4])])
assert not bad,bad
r=w/'render';r.mkdir(exist_ok=True)
exe=r'C:\Users\31937\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin\pdftoppm.exe'
subprocess.run([exe,'-f','35','-l','42','-scale-to','1600','-png',str(p),str(r/'page')],check=True,creationflags=subprocess.CREATE_NO_WINDOW)
files=sorted(r.glob('page-*.png'))
for start in [0,4]:
 sheet=Image.new('RGB',(1400,2020),'#dddddd')
 for j,f in enumerate(files[start:start+4]):
  im=Image.open(f).convert('RGB');im.thumbnail((690,975)); x=(j%2)*700;y=(j//2)*1010
  sheet.paste(im,(x,y+25));ImageDraw.Draw(sheet).text((x+10,y+5),f.stem,fill='black')
 sheet.save(r/f'contact_{start//4+1}.png')
q={'pages':len(d),'chapter_pages':[35,40],'unchanged_pages_1_34':True,'visuals':27,'portrait_all':True,'bounds_issues':bad,'v25_preserved':True,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
(w/'qa.json').write_text(json.dumps(q,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(q))