from pathlib import Path
import fitz,json,hashlib
from PIL import Image,ImageOps,ImageDraw
B=Path(__file__).parent;D=B/'closeout_delivery';f=D/'CD4_IKZF2_STAT1_stage_delivery_zh.pdf';pdf=fitz.open(f);report=[];thumbs=[]
for i,pg in enumerate(pdf):
 text=pg.get_text();blocks=pg.get_text('dict')['blocks'];outside=[];fonts=set()
 for b in blocks:
  for l in b.get('lines',[]):
   for sp in l['spans']:
    fonts.add(sp['font']);x0,y0,x1,y1=sp['bbox']
    if x0<0 or y0<0 or x1>pg.rect.width+.5 or y1>pg.rect.height+.5:outside.append(sp['text'])
 pix=pg.get_pixmap(matrix=fitz.Matrix(1.15,1.15));pix.save(D/f'qa_page_{i+1:02d}.png')
 im=Image.open(D/f'qa_page_{i+1:02d}.png').convert('RGB');im.thumbnail((420,595));can=Image.new('RGB',(440,625),'#dce4e8');can.paste(im,((440-im.width)//2,10));ImageDraw.Draw(can).text((10,607),f'Page {i+1}',fill='black');thumbs.append(can)
 report.append(dict(page=i+1,characters=len(text),out_of_page=outside,fonts=sorted(fonts),beginning=text[:100]))
sheet=Image.new('RGB',(440*3,625*((len(thumbs)+2)//3)),'white')
for i,im in enumerate(thumbs):sheet.paste(im,((i%3)*440,(i//3)*625))
sheet.save(D/'QA_contact_sheet.jpg')
result=dict(pages=len(pdf),pages_detail=report,sha256=hashlib.sha256(f.read_bytes()).hexdigest(),automated_bounds_pass=not any(x['out_of_page'] for x in report))
(D/'QA_render.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(result,ensure_ascii=False))
