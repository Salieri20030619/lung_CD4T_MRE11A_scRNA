from pathlib import Path
import json,hashlib,subprocess,pymupdf as fitz,pandas as pd
from PIL import Image,ImageDraw
O=Path(__file__).parent;d=json.loads((O/'report_delivery.json').read_text(encoding='utf-8'));doc=fitz.open(d['pdf']);index=json.loads((O/'report_page_index.json').read_text(encoding='utf-8'));checks=[]
for i,p in enumerate(doc):
 assert p.rect.height>p.rect.width
 bad=[b[:4] for b in p.get_text('blocks') if b[0]<39 or b[2]>p.rect.width-37 or b[1]<24 or b[3]>p.rect.height-15];checks.append({'page':i+1,'images':len(p.get_images()),'bounds':bad,'chars':len(p.get_text())})
assert sum(x['images'] for x in checks)==24
pages=sorted(set([9]+[x['page'] for x in index if '候选筛选' in x['heading']]+[x['page'] for x in index if '候选与验证' in x['heading']]))
R=O/'render';R.mkdir(exist_ok=True);pop=r'C:\Users\31937\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin\pdftoppm.exe'
for n in pages:subprocess.run([pop,'-f',str(n),'-l',str(n),'-scale-to','1600','-png',d['pdf'],str(R/('p'+str(n).zfill(2)))],check=True,creationflags=subprocess.CREATE_NO_WINDOW)
sheet=Image.new('RGB',(1440,1020*((len(pages)+1)//2)),'#dce4e8');draw=ImageDraw.Draw(sheet)
for i,n in enumerate(pages):
 f=list(R.glob('p'+str(n).zfill(2)+'-*.png'))[0];im=Image.open(f).convert('RGB');im.thumbnail((700,990));x=i%2*720;y=i//2*1020;sheet.paste(im,(x,y+20));draw.text((x+4,y+4),'Page '+str(n),fill='black')
sheet.save(R/'changed_pages.jpg',quality=92)
q={'pages':len(doc),'figures':24,'all_portrait':True,'page_checks':checks,'changed_pages_rendered':pages,'screening':json.loads((O/'screening_QA.json').read_text()),'labels':json.loads((O/'F05_label_QA.json').read_text())};(O/'QA.json').write_text(json.dumps(q,ensure_ascii=False,indent=2),encoding='utf-8');print('PAGES',len(doc),'CHANGED',pages,'BOUNDS',[(x['page'],x['bounds']) for x in checks if x['bounds']]);print('HEADINGS',[(x['page'],x['heading']) for x in index if '候选筛选' in x['heading']])
supp=pd.read_csv(O.parent/'lung_three_cohort_20261017_v16/final_tables/历史候选逐项回溯.tsv',sep='\t');print('SUPP_COLS',supp.columns.tolist());supp.to_csv(O/'tables/supplemental42_evidence.tsv',sep='\t',index=False,encoding='utf-8')