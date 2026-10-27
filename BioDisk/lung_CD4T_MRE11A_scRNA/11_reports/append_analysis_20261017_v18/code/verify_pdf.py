from pathlib import Path
import json,hashlib,subprocess,fitz
from PIL import Image,ImageDraw
O=Path(r'C:\Users\31937\lung_append_analysis_20261017_v18');j=json.loads((O/'delivery.json').read_text(encoding='utf-8'));old=fitz.open(j['baseline_pdf']);doc=fitz.open(j['pdf']);R=O/'render';R.mkdir(exist_ok=True);same=[]
assert len(old)==9 and len(doc)==9+j['appendix_pages']
for i in range(9):
 assert old[i].get_text()==doc[i].get_text();a=old[i].get_pixmap(matrix=fitz.Matrix(1,1),alpha=False);b=doc[i].get_pixmap(matrix=fitz.Matrix(1,1),alpha=False);assert a.samples==b.samples;same.append(i+1)
for i in range(9,len(doc)):
 txt=doc[i].get_text();assert len(txt)>100 and '\ufffd' not in txt
 for b in doc[i].get_text('blocks'):assert b[0]>=35 and b[2]<=doc[i].rect.width-35 and b[1]>=25 and b[3]<=doc[i].rect.height-15,(i+1,b[:4])
exe=r'C:\Users\31937\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin\pdftoppm.exe';subprocess.run([exe,'-f','10','-l',str(len(doc)),'-r','90','-png',j['pdf'],str(R/'page')],check=True,capture_output=True,creationflags=subprocess.CREATE_NO_WINDOW);files=sorted(R.glob('page-*.png'));assert len(files)==j['appendix_pages']
for k in range(0,len(files),4):
 canvas=Image.new('RGB',(1000,1500),'#dfe7e9');draw=ImageDraw.Draw(canvas)
 for n,f in enumerate(files[k:k+4]):
  im=Image.open(f);im.thumbnail((480,710));xx=(n%2)*500+10;yy=(n//2)*750+25;canvas.paste(im,(xx,yy));draw.text((xx,yy-18),f.stem,fill='black')
 canvas.save(R/f'contact_{k//4+1}.png')
out={'combined_pages':len(doc),'appendix_pages':j['appendix_pages'],'baseline_text_and_pixels_identical_pages':same,'new_pages_bounds_ok':True,'all_new_pages_rendered':True};(O/'pdf_QA.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out))