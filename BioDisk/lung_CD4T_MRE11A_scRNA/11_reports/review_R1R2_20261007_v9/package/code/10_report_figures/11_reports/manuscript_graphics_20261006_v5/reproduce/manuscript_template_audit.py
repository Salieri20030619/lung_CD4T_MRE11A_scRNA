from pathlib import Path
from collections import Counter
import pymupdf as fitz,json,re
from PIL import Image,ImageDraw
src=Path(r'C:\Users\31937\Downloads\26210\Phase327_摘要措辞修订_v3\26210_manuscript_immune_integrated_v3.pdf')
D=Path(__file__).parent/'manuscript_delivery';D.mkdir(exist_ok=True)
pdf=fitz.open(src);fonts=Counter();sizes=Counter();heads=[];page_info=[]
for i,page in enumerate(pdf):
 text=page.get_text();page_info.append({'page':i+1,'characters':len(text),'images':len(page.get_images()),'first_lines':text.splitlines()[:3]})
 for b in page.get_text('dict')['blocks']:
  for line in b.get('lines',[]):
   spans=line['spans'];t=''.join(s['text'] for s in spans).strip()
   for s in spans:fonts[s['font']]+=len(s['text']);sizes[round(s['size'],1)]+=len(s['text'])
   if t and (re.match(r'^(\d+(?:\.\d+)*\s+|摘要|关键词|参考文献|补充|图\s*\d|表\s*\d)',t) or any(s['size']>=15 for s in spans)):
    heads.append({'page':i+1,'text':t[:180]})
out={'source_read_only':str(src),'pages':len(pdf),'page_size':list(pdf[0].rect),'fonts':fonts.most_common(12),'sizes':sizes.most_common(12),'headings':heads,'page_info':page_info,'use':'presentation style only; exclude study conclusions and external-project evidence'}
(D/'template_style_audit.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:out[k] for k in ['pages','page_size','fonts','sizes','headings']},ensure_ascii=False))
print('FIRST_PAGE_PROSE_FOR_STYLE_ONLY\n'+pdf[0].get_text()[:4000])
selected=sorted(set([0,1,2,min(5,len(pdf)-1),min(10,len(pdf)-1),len(pdf)-1]));thumbs=[]
for i in selected:
 f=D/f'template_page_{i+1:02d}.png';pdf[i].get_pixmap(matrix=fitz.Matrix(.9,.9)).save(f)
 im=Image.open(f).convert('RGB');im.thumbnail((420,595));can=Image.new('RGB',(440,625),'#e4e8ec');can.paste(im,((440-im.width)//2,8));ImageDraw.Draw(can).text((8,606),f'Reference page {i+1}',fill='black');thumbs.append(can)
sheet=Image.new('RGB',(1320,1250),'white')
for j,im in enumerate(thumbs):sheet.paste(im,((j%3)*440,(j//3)*625))
sheet.save(D/'template_style_contact.jpg')
