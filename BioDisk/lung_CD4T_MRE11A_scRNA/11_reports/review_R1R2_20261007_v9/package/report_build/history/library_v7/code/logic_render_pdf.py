from pathlib import Path
import subprocess, json, re, hashlib
from lxml import html
import pymupdf as fitz
from PIL import Image, ImageOps, ImageDraw

B=Path(__file__).resolve().parent
D=B/'logic_revision_20261007'
source=D/'provenance/logic_render_source.html'
candidate=D/'logic_pdf_candidate.pdf'
exe=r'C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe'
args=[exe,'--headless','--disable-gpu','--no-first-run','--no-default-browser-check',
      '--no-pdf-header-footer','--user-data-dir='+str(B/'logic_edge_private'),
      '--print-to-pdf='+str(candidate),source.resolve().as_uri()]
print('HEADLESS_PDF_START',flush=True)
r=subprocess.run(args,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=120,creationflags=subprocess.CREATE_NO_WINDOW)
(D/'QA/renderer_log.txt').write_text(r.stderr,encoding='utf-8')
assert r.returncode==0 and candidate.exists(),r.returncode
print('RENDER_COMPLETE',candidate.stat().st_size,flush=True)
pdf=fitz.open(candidate)
tree=html.fromstring(source.read_text(encoding='utf-8'))
chunks=[]
for e in tree.xpath('//body/*'):
    if e.tag=='table':chunks.extend(c.text_content() for c in e.xpath('.//th|.//td'))
    elif e.text_content().strip():chunks.append(e.text_content())
normalize=lambda s:re.sub(r'\s+','',s).replace('\u00ad','')
full=normalize('\n'.join(p.get_text() for p in pdf))
missing=[s for s in chunks if normalize(s) not in full]
(D/'QA/missing_pdf_text.json').write_text(json.dumps(missing,ensure_ascii=False,indent=2),encoding='utf-8')
assert not missing,[s[:120] for s in missing]
page_records=[]
for i,p in enumerate(pdf):
    blocks=p.get_text('blocks')
    outside=[list(b[:4])+[b[4][:100]] for b in blocks if b[0]<45 or b[1]<45 or b[2]>p.rect.width-45 or b[3]>p.rect.height-45]
    page_records.append({'page':i+1,'text_blocks':len(blocks),'outside_main_margins':outside})
    p.insert_textbox(fitz.Rect(200,p.rect.height-31,p.rect.width-200,p.rect.height-17),str(i+1),fontsize=9,fontname='tiro',align=1,color=(.1,.1,.1))
pdf.set_metadata({'title':tree.xpath('//title')[0].text,'author':'',
 'subject':'CD4 T-cell phenotype associations in CWP patients and mouse silica models',
 'creator':'Local HTML rendering with Microsoft Edge; content and pagination verified',
 'keywords':'CD4;UMAP;GSEA;IKZF2;STAT1'})
out=D/'矽肺研究中巨噬细胞与CD4T细胞表型关联分析_20261007_最终逻辑核查版.pdf'
pdf.save(out,garbage=4,deflate=True)
proof=D/'QA/rendered_pages';proof.mkdir(exist_ok=True)
thumbs=[]
for i,p in enumerate(pdf):
    pix=p.get_pixmap(matrix=fitz.Matrix(1.3,1.3),alpha=False)
    pix.save(proof/f'page_{i+1:02}.png')
    img=Image.open(proof/f'page_{i+1:02}.png').convert('RGB')
    img.thumbnail((300,425))
    canvas=Image.new('RGB',(320,455),'#dddddd');canvas.paste(img,((320-img.width)//2,22))
    ImageDraw.Draw(canvas).text((10,5),f'Page {i+1}',fill='black')
    thumbs.append(canvas)
for j in range(0,len(thumbs),8):
    sheet=Image.new('RGB',(1280,910),'white')
    for k,img in enumerate(thumbs[j:j+8]):sheet.paste(img,((k%4)*320,(k//4)*455))
    sheet.save(D/'QA'/f'contact_{j//8+1}.jpg',quality=90)
result={'renderer':'Local headless Microsoft Edge from verified HTML', 'pages':len(pdf),
 'html_text_blocks_and_cells_verified':len(chunks),'missing_text_blocks':missing,
 'page_numbers_added':True,'page_records':page_records,
 'pdf':{'name':out.name,'bytes':out.stat().st_size,'sha256':hashlib.sha256(out.read_bytes()).hexdigest()}}
(D/'QA/PDF_content_parity.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in result.items() if k!='page_records'},ensure_ascii=True))
