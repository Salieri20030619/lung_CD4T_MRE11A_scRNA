from pathlib import Path
from docx import Document
from docx.oxml.ns import qn
from lxml import etree
import html,base64,subprocess,json,time
B=Path(__file__).parent;D=B/'review_delivery';doc=Document(D/'internal/revised_source.docx')
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','a':'http://schemas.openxmlformats.org/drawingml/2006/main'}
def inline(e):
 out=[]
 for x in e:
  if x.tag==qn('w:r'):
   s=''
   for z in x:
    if z.tag==qn('w:t'):s+=html.escape(z.text or '')
    elif z.tag==qn('w:br') and z.get(qn('w:type'))!='page':s+='<br>'
    elif z.tag==qn('w:tab'):s+='&emsp;'
    elif z.tag==qn('w:drawing'):
     for a in z.findall('.//a:blip',ns):
      part=doc.part.related_parts[a.get(qn('r:embed'))];s+='<img src="data:'+part.content_type+';base64,'+base64.b64encode(part.blob).decode()+'">'
   pr=x.find(qn('w:rPr'));styles=[]
   if pr is not None:
    for tag,css in [('b','font-weight:bold'),('i','font-style:italic')]:
     a=pr.find(qn('w:'+tag))
     if a is not None and a.get(qn('w:val')) not in ['0','false']:styles.append(css)
    a=pr.find(qn('w:sz'))
    if a is not None:styles.append('font-size:'+str(int(a.get(qn('w:val')))/2)+'pt')
   out.append('<span style="'+';'.join(styles)+'">'+s+'</span>')
  elif x.tag in [qn('w:hyperlink'),qn('w:smartTag')]:out.append(inline(x))
 return ''.join(out)
parts=[];body=list(doc.element.body);i=0
while i<len(body):
 e=body[i]
 if e.tag==qn('w:p'):
  br=e.findall('.//w:br',ns)
  if any(x.get(qn('w:type'))=='page' for x in br):parts.append('<div class="pagebreak"></div>')
  st=e.find('w:pPr/w:pStyle',ns);style=st.get(qn('w:val')) if st is not None else 'Normal';s=inline(e)
  if '<img ' in s:
   cap=inline(body[i+1]);parts.append('<div class="figure"><p class="image">'+s+'</p><p class="Caption">'+cap+'</p></div>');i+=1
  elif s:
   if style.startswith('Heading'):parts.append('<h'+style[-1]+'>'+s+'</h'+style[-1]+'>')
   else:parts.append('<p class="'+style+'">'+s+'</p>')
 elif e.tag==qn('w:tbl'):
  rows=e.findall('w:tr',ns);tr=[]
  for j,row in enumerate(rows):
   cells=[]
   for c in row.findall('w:tc',ns):
    texts=['<p>'+inline(p)+'</p>' for p in c.findall('w:p',ns)];w=c.find('w:tcPr/w:tcW',ns)
    width=str(float(w.get(qn('w:w')))/1440)+'in' if w is not None else 'auto'
    cells.append('<'+('th' if j==0 else'td')+' style="width:'+width+'">'+''.join(texts)+'</'+('th' if j==0 else'td')+'>')
   tr.append('<tr>'+''.join(cells)+'</tr>')
  parts.append('<table><thead>'+tr[0]+'</thead><tbody>'+''.join(tr[1:])+'</tbody></table>')
 i+=1
css='''@page{size:A4;margin:18mm 18.5mm 18mm}*{box-sizing:border-box}body{margin:0;color:#111;font-family:"Times New Roman",SimSun,serif;font-size:10.5pt;line-height:1.18}p{margin:0 0 5pt;text-align:justify;text-indent:21pt;orphans:2;widows:2}h1,h2,h3{font-family:SimHei,serif;line-height:1.18;margin:11pt 0 5pt;break-after:avoid;font-weight:bold}h1{font-size:13pt}h2{font-size:11.5pt}h3{font-size:10.5pt}.Title{font-family:SimHei,serif;font-size:18pt;text-align:center;text-indent:0}.Subtitle{font-size:10pt;text-align:center;text-indent:0}.Caption{font-size:8.6pt;text-indent:0;text-align:left;line-height:1.18}.figure{break-inside:avoid;margin:5pt 0 8pt}.image{margin:0 0 4pt;text-indent:0;text-align:center}.image img{display:block;width:171.5mm;max-width:100%;height:auto}.pagebreak{break-before:page}table{border-collapse:collapse;table-layout:fixed;width:172mm;font-size:8.3pt;line-height:1.02;margin:3pt 0 5pt}thead{display:table-header-group}tr{break-inside:avoid}th,td{border:none;padding:2.5pt 3pt;vertical-align:top;overflow-wrap:anywhere}th{background:#e8eff4}td p,th p{font-size:8.3pt;text-indent:0;text-align:left;line-height:1.02;margin:1.7pt 0}.Caption:has(+table){break-after:avoid}'''
page='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>CD4 scientific report</title><style>'+css+'</style><body>'+''.join(parts)+'</body></html>'
f=D/'provenance/review_render_source.html';f.write_text(page,encoding='utf8')
exe=r'C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe';out=D/'review_pdf_candidate.pdf';profile=B/'review_edge_private'
args=[exe,'--headless','--disable-gpu','--no-first-run','--no-default-browser-check','--no-pdf-header-footer','--user-data-dir='+str(profile),'--print-to-pdf='+str(out),f.resolve().as_uri()]
print('HEADLESS_PDF_START',flush=True)
r=subprocess.run(args,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=120,creationflags=subprocess.CREATE_NO_WINDOW)
print('HEADLESS_PDF_RESULT',r.returncode,out.exists(),out.stat().st_size if out.exists() else 0,flush=True)
(D/'QA/renderer_log.txt').write_text(r.stderr,encoding='utf8')
assert r.returncode==0 and out.exists()
