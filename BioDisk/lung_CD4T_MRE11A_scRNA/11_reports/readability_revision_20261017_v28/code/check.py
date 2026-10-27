from pathlib import Path
import fitz,json
O=Path(r'C:\Users\31937\lung_readability_revision_20261017_v28');p=json.loads((O/'report_delivery.json').read_text(encoding='utf-8'))['pdf'];d=fitz.open(p)
for i,pg in enumerate(d):
 t=pg.get_text();print(i+1,len(t),t.splitlines()[2:3])
print('HEADINGS',json.loads((O/'report_page_index.json').read_text(encoding='utf-8'))[-12:])
bad=[]
for i,pg in enumerate(d):
 for b in pg.get_text('blocks'):
  if b[0]<39 or b[2]>pg.rect.width-37 or b[1]<24 or b[3]>pg.rect.height-15:bad.append([i+1,b[:4]])
print('BOUNDS',bad)
(O/'extracted.txt').write_text('\n'.join(pg.get_text() for pg in d),encoding='utf-8')