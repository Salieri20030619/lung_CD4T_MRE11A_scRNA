"""Rebuild only Table 2 from a hash-pinned Library v7 PDF and keyed source rows.
No statistical model, target selection, or multiple-testing calculation occurs.
Run: python r1r2_build_pdf.py --package PATH_TO_EXTRACTED_PACKAGE
"""
from pathlib import Path
import argparse,json,hashlib,io,re,base64,copy
from lxml import html,etree
import pymupdf as fitz
import numpy as np
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from PIL import Image,ImageDraw

parser=argparse.ArgumentParser();parser.add_argument('--package',type=Path)
args=parser.parse_args()
K=args.package.resolve() if args.package else Path(__file__).resolve().parent/'r1r2_revision_20261007/package'
C=K/'report_build/current';H=K/'report_build/history/library_v7'
sha=lambda data:hashlib.sha256(data).hexdigest()
baseline=H/'outputs/library_v7.pdf';assert sha(baseline.read_bytes())=='fa19fba5f5b6d3f0a946d99148b822627894847414acab0c5a8185aa88f38f3a'
source=H/'outputs/logic_render_source.html'
tree=html.fromstring(source.read_text(encoding='utf-8'));original=copy.deepcopy(tree)
nodes=tree.xpath('//body/*');table=nodes[121]
assert table.tag=='table' and nodes[120].text_content().startswith('表2')
oldrows=[[c.text_content() for c in r.xpath('./th|./td')] for r in table.xpath('.//tr')]
assert oldrows[0]==['接收端/来源','校正AUPR','单对LR','β/γ链','四基因病例','关联n','P','原BH']
extract=json.loads((C/'inputs/IL15_table2_source_rows.json').read_text(encoding='utf-8'))
assert extract['selectors']=={'filter':'All','model':'Source','test':'Nominal05','direction':'Up','ligand':'IL15'}
records={(r['receiver'],r['sender']):r for r in extract['rows']};assert len(records)==6
recv={'Core':'Core','IL7R':'CD4-IL7R'};send={'AM':'MacAM','MM':'MonoMac','MW':'MacWide'}
audit=[];newrows=[oldrows[0][:6]+['偏秩ρ']+oldrows[0][6:]]
for row in oldrows[1:]:
    r,s=row[0].split('/');key=(recv[r],send[s]);v=records[key]
    assert all(v[k]==x for k,x in extract['selectors'].items())
    assert abs(float(row[1])-v['aupr_corrected'])<=0.0000005
    if v['association_rho'] is None:
        assert v['sender']=='MacAM' and v['status']=='NO_SOURCE_LR_GATE_SUPPORT'
        assert row[5:] == ['—','未计算','未计算']
        assert all(v[k] is None for k in ['association_n','association_P','association_FDR'])
        rho='未计算'
    else:
        assert int(row[5])==int(v['association_n'])
        assert abs(float(row[6])-v['association_P'])<=0.0000005
        assert abs(float(row[7])-v['association_FDR'])<=0.0000005
        rho=f"{v['association_rho']:+.6f}"
    newrows.append(row[:6]+[rho]+row[6:])
    audit.append({'display_key':row[0],'source_key':key,'full_endpoint':extract['selectors'],'source_values':v,'added_rho':rho,'existing_n_P_BH_unchanged':True})
widths=[72,66,45,43,55,34,54,57,61.4884033203125]
assert abs(sum(widths)-487.4884033203125)<1e-8
for ri,row in enumerate(table.xpath('.//tr')):
    cell=etree.Element('th' if ri==0 else 'td');etree.SubElement(cell,'p').text=newrows[ri][6]
    row.insert(6,cell)
    for ci,c in enumerate(row.xpath('./th|./td')):c.set('style',f'width:{widths[ci]}pt')
for i,(a,b) in enumerate(zip(original.xpath('//body/*'),nodes)):
    if i!=121:assert etree.tostring(a)==etree.tostring(b),i
finalhtml=C/'sources/final_report.html'
finalhtml.write_text('<!doctype html>\n'+etree.tostring(tree,encoding='unicode',method='html'),encoding='utf-8')
newpdf=fitz.open(baseline);basepdf=fitz.open(baseline)
matches=[i for i,p in enumerate(newpdf) if '表2' in p.get_text() and 'Core/AM' in p.get_text()]
assert len(matches)==1;pageindex=matches[0];p=newpdf[pageindex]
blocks=[b for b in p.get_text('dict')['blocks'] if b['type']==0 and 105<b['bbox'][1]<230]
assert len(blocks)==7
baselines=[b['lines'][0]['spans'][0]['origin'][1] for b in blocks]
rect=fitz.Rect(52.5,107.75,540.25,229.0)
# Font resources are local operating-system fonts, recorded by hash below.
fonts={'Song':Path(r'C:\Windows\Fonts\simsun.ttc'),'Times':Path(r'C:\Windows\Fonts\times.ttf')}
pdfmetrics.registerFont(TTFont('R1R2Song',str(fonts['Song']),subfontIndex=0))
pdfmetrics.registerFont(TTFont('R1R2Times',str(fonts['Times'])))
buf=io.BytesIO();c=canvas.Canvas(buf,pagesize=(p.rect.width,p.rect.height),pageCompression=1,invariant=1)
c.setFillColorRGB(0.9098,0.9373,0.9569);c.rect(52.5,p.rect.height-125.248046875,487.4884033203125,17.2490234375,fill=1,stroke=0)
c.setFillColorRGB(0.066667,0.066667,0.066667)
for ri,row in enumerate(newrows):
    x=52.5
    for ci,value in enumerate(row):
        font='R1R2Song' if re.search(r'[^\x00-\x7f]',value) else 'R1R2Times'
        size=8.0
        assert pdfmetrics.stringWidth(value,font,size)<widths[ci]-6,(value,widths[ci])
        c.setFont(font,size);c.drawString(x+3,p.rect.height-baselines[ri],value);x+=widths[ci]
c.showPage();c.save();overlaybytes=buf.getvalue()
(C/'sources/table2_overlay.pdf').write_bytes(overlaybytes)
overlay=fitz.open(stream=overlaybytes,filetype='pdf')
p.add_redact_annot(rect,fill=(1,1,1));p.apply_redactions(images=0,graphics=0,text=0)
p.show_pdf_page(rect,overlay,0,clip=rect,keep_proportion=False)
out=K/'report/CD4T_review_R1R2_20261007.pdf'
newpdf.save(out,garbage=4,deflate=True)
newpdf.close();newpdf=fitz.open(out)
norm=lambda s:re.sub(r'\s+','',s).replace('\u00ad','')
chunks=[]
for e in nodes:
    if e.tag=='table':chunks.extend(''.join(c.itertext()) for c in e.xpath('.//th|.//td'))
    elif e.text_content().strip():chunks.append(e.text_content())
alltext=norm('\n'.join(p.get_text() for p in newpdf));missing=[t for t in chunks if norm(t) not in alltext]
assert not missing,[s[:80] for s in missing]
assert len(newpdf)==len(basepdf)
pixelchecks=[];proof=C/'QA/rendered_pages';proof.mkdir(exist_ok=True);thumbs=[]
for i,(old,new) in enumerate(zip(basepdf,newpdf)):
    a=old.get_pixmap(matrix=fitz.Matrix(1.3,1.3),alpha=False);b=new.get_pixmap(matrix=fitz.Matrix(1.3,1.3),alpha=False)
    av=np.frombuffer(a.samples,dtype=np.uint8).reshape(a.height,a.width,3);bv=np.frombuffer(b.samples,dtype=np.uint8).reshape(b.height,b.width,3)
    changed=np.any(av!=bv,axis=2)
    if i==pageindex:
        safe=changed.copy();rr=rect*1.3;safe[max(0,int(rr.y0)-2):int(rr.y1)+3,max(0,int(rr.x0)-2):int(rr.x1)+3]=False
        assert not safe.any(),'Changes outside Table 2 rectangle'
        def outside(page):return [(tuple(b[:4]),b[4]) for b in page.get_text('blocks') if not fitz.Rect(b[:4]).intersects(rect)]
        assert outside(old)==outside(new),'Non-target text or positions changed'
    else:
        assert not changed.any(),i
        assert old.get_text()==new.get_text(),i
    pixelchecks.append({'page':i+1,'changed_pixels':int(changed.sum()),'outside_target_changed_pixels':0})
    b.save(proof/f'page_{i+1:02}.png')
    im=Image.open(proof/f'page_{i+1:02}.png').convert('RGB');im.thumbnail((300,425))
    tile=Image.new('RGB',(320,455),'#dddddd');tile.paste(im,((320-im.width)//2,22));ImageDraw.Draw(tile).text((10,5),f'Page {i+1}',fill='black');thumbs.append(tile)
for j in range(0,len(thumbs),8):
    im=Image.new('RGB',(1280,910),'white')
    for k,t in enumerate(thumbs[j:j+8]):im.paste(t,((k%4)*320,(k//4)*455))
    im.save(C/'QA'/f'contact_{j//8+1}.jpg',quality=90)
imagehash=lambda doc:sorted(set(sha(fitz.Pixmap(doc,x[0]).samples) for p in doc for x in p.get_images(full=True)))
assert imagehash(basepdf)==imagehash(newpdf) and len(imagehash(newpdf))==14
assert basepdf.metadata==newpdf.metadata
assert '其零检出记录主要反映测量覆盖' in alltext and 'FDR记为未计算' in alltext
parity={'pdf_relative_path':str(out.relative_to(K)).replace('\\','/'),'pages':len(newpdf),'pdf_sha256':sha(out.read_bytes()),'pdf_bytes':out.stat().st_size,'html_text_blocks_and_cells_verified':len(chunks),'missing_text_blocks':missing,'source_html_sha256':sha(finalhtml.read_bytes()),'renderer':'Hash-pinned v7 PDF; ReportLab table overlay; PyMuPDF bounded replacement','page_numbers':'Preserved from baseline'}
manifest={'allowed_changes':['R1 delivery provenance','R2 Table 2 rho column'],'non_target_changes':[], 'R3_BH_not_changed':True,'G288_sentence_not_changed':True,'title_metadata_unchanged':True,'baseline_pdf_sha256':sha(baseline.read_bytes()),'target_page':pageindex+1,'target_rectangle':list(rect),'table2_source_sha256':extract['source_sha256'],'font_resources':{k:{'path':str(v),'sha256':sha(v.read_bytes())} for k,v in fonts.items()},'table2_before':oldrows,'table2_after':newrows,'matched_rows':audit,'scientific_models_run':0,'P_or_BH_computed':0}
for path,value in [(C/'QA/PDF_content_parity.json',parity),(C/'QA/non_target_invariance.json',{'pages':pixelchecks,'all_other_pages_pixel_identical':True,'target_page_outside_table_pixel_identical':True,'figures':14,'image_pixels_identical':True,'unchanged_other_tables':17,'all_existing_table2_cells_unchanged':True}),(C/'QA/table2_source_validation.json',audit),(C/'revision_manifest.json',manifest)]:path.write_text(json.dumps(value,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(parity,ensure_ascii=True));print('ONLY_TABLE2_PIXELS_CHANGED',pageindex+1,'RHO_ROWS',len(audit),flush=True)
