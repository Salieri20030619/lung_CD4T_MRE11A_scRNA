from pathlib import Path
import json,subprocess,pymupdf as fitz
from PIL import Image,ImageDraw
O=Path(__file__).parent;B=O.parent/'lung_client_rebuild_20261017_v20';d=json.loads((O/'report_delivery.json').read_text(encoding='utf-8'));doc=fitz.open(d['pdf']);pages=[];checks=[]
for i,p in enumerate(doc):
 assert p.rect.height>p.rect.width
 text=p.get_text();bad=[b[:4] for b in p.get_text('blocks') if b[0]<39 or b[2]>p.rect.width-37 or b[1]<24 or b[3]>p.rect.height-15];checks.append({'page':i+1,'images':len(p.get_images()),'bounds':bad})
 if any(s in text for s in ['示意图A。','示意图B。','示意图C。','图5。','图6。','图1。','图3。']):pages.append(i+1)
assert len(doc)==37 and sum(x['images'] for x in checks)==27
for id in ['F05','F06']:assert (O/'plot_data'/(id+'.tsv')).read_bytes()==(B/'plot_data'/(id+'.tsv')).read_bytes()
for id in ['F05','F06']:
 tx=''.join(p.get_text() for p in fitz.open(O/'figures'/(id+'.pdf')));assert 'Joint ten' not in tx and 'GSE314943' in tx and 'GSE183682' in tx
pages=sorted(set(pages+[2]));R=O/'render';R.mkdir(exist_ok=True);pop=r'C:\Users\31937\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin\pdftoppm.exe'
for n in pages:subprocess.run([pop,'-f',str(n),'-l',str(n),'-scale-to','1500','-png',d['pdf'],str(R/('p'+str(n).zfill(2)))],check=True,creationflags=subprocess.CREATE_NO_WINDOW)
for k in range(0,len(pages),4):
 subset=pages[k:k+4];sheet=Image.new('RGB',(1440,1020*((len(subset)+1)//2)),'#dce4e8');dr=ImageDraw.Draw(sheet)
 for j,n in enumerate(subset):
  f=list(R.glob('p'+str(n).zfill(2)+'-*.png'))[0];im=Image.open(f).convert('RGB');im.thumbnail((700,990));x=j%2*720;y=j//2*1020;sheet.paste(im,(x,y+20));dr.text((x+4,y+4),'Page '+str(n),fill='black')
 sheet.save(R/('contact_'+str(k//4+1)+'.jpg'),quality=92)
q={'pages':len(doc),'numbered_figures':24,'explanatory_displays':3,'all_portrait':True,'total_images':27,'page_checks':checks,'changed_pages':pages,'differential_inputs_unchanged':True,'cohort_counts':json.loads((O/'cohort_QA.json').read_text(encoding='utf-8'))};(O/'QA.json').write_text(json.dumps(q,ensure_ascii=False,indent=2),encoding='utf-8');print('PAGES',len(doc),'RENDERED',pages,'BOUNDS',[(x['page'],x['bounds']) for x in checks if x['bounds']])