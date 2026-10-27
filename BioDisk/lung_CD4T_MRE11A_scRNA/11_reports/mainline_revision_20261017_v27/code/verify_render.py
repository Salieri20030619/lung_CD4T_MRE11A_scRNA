from pathlib import Path
import json,fitz,subprocess
from PIL import Image,ImageDraw
O=Path(r'C:\Users\31937\lung_mainline_revision_20261017_v27');p=Path(json.loads((O/'report_delivery.json').read_text(encoding='utf-8'))['pdf']);d=fitz.open(p);r=O/'render';r.mkdir(exist_ok=True)
text='\n'.join(x.get_text() for x in d);(O/'pdf_extracted.txt').write_text(text,encoding='utf-8')
for i,pg in enumerate(d):
 z=pg.get_text();print(i+1,len(z),z.splitlines()[2:4])
bad=[]
for i,pg in enumerate(d):
 for z in pg.get_text('blocks'):
  if z[0]<39 or z[2]>pg.rect.width-37 or z[1]<24 or z[3]>pg.rect.height-15:bad.append([i+1,list(z[:4])])
print('BOUNDS',bad);print('VISUALS',sum(len(pg.get_images()) for pg in d));print('BANNED',[x for x in ['历史基因','继续硬筛','二元创新性P','巨噬入口','核心必做','保留原8条','数据索引：','min10敏感性','旧版','v18','v19'] if x in text])
exe=r'C:\Users\31937\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin\pdftoppm.exe'
subprocess.run([exe,'-scale-to','1400','-png',str(p),str(r/'page')],check=True,creationflags=subprocess.CREATE_NO_WINDOW)
files=sorted(r.glob('page-*.png'))
for start in range(0,len(files),6):
 sheet=Image.new('RGB',(1500,1450),'#ddd');dr=ImageDraw.Draw(sheet)
 for j,f in enumerate(files[start:start+6]):
  im=Image.open(f).convert('RGB');im.thumbnail((490,690));x=(j%3)*500;y=(j//3)*725;sheet.paste(im,(x,y+25));dr.text((x+8,y+5),f.stem,fill='black')
 sheet.save(r/f'contact_{start//6+1:02d}.png')
(O/'layout_check.json').write_text(json.dumps({'pages':len(d),'bounds_issues':bad,'visuals':sum(len(pg.get_images()) for pg in d)},indent=2),encoding='utf-8')