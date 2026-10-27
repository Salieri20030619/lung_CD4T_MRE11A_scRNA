from pathlib import Path
import fitz,json,re
O=Path(r'C:\Users\31937\lung_selection_joint_revision_20261017_v30');p=json.loads((O/'report_delivery.json').read_text(encoding='utf-8'))['pdf'];d=fitz.open(p);bad=[]
for i,pg in enumerate(d):
 print(i+1,len(pg.get_text()),pg.get_text().splitlines()[2:3])
 for b in pg.get_text('blocks'):
  if b[0]<39 or b[2]>pg.rect.width-37 or b[1]<24 or b[3]>pg.rect.height-15:bad.append([i+1,b[:4],b[4][:150]])
print('BOUNDS',bad,'VISUALS',sum(len(pg.get_images()) for pg in d));(O/'extracted.txt').write_text('\n'.join(pg.get_text() for pg in d),encoding='utf-8')