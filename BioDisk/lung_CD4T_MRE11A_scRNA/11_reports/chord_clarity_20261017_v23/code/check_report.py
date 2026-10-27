from pathlib import Path
import json,pymupdf as fitz,pandas as pd
O=Path(__file__).parent;d=json.loads((O/'report_delivery.json').read_text(encoding='utf-8'));doc=fitz.open(d['pdf']);print('PAGES',len(doc));print([(i+1,len(p.get_text()),len(p.get_images())) for i,p in enumerate(doc)]);flags=[]
for i,p in enumerate(doc):
 assert p.rect.height>p.rect.width
 for b in p.get_text('blocks'):
  if b[0]<39 or b[2]>p.rect.width-37 or b[1]<24 or b[3]>p.rect.height-15:flags.append([i+1,list(b[:4])])
print('BOUNDS',flags);print('PAGE1_END',doc[0].get_text()[-400:]);(O/'basic_QA.json').write_text(json.dumps({'pages':len(doc),'figures':sum(len(p.get_images()) for p in doc),'bounds':flags}),encoding='utf-8')