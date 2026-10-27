from pathlib import Path
from collections import Counter
import hashlib, io, json, math, re, zipfile
import numpy as np
import pymupdf as fitz

B = Path(__file__).resolve().parent
D = B / 'title_revision_20261007'
for part in ['QA', 'provenance']:
    (D / part).mkdir(parents=True, exist_ok=True)
source = B / 'review_delivery/review_report.pdf'
old_zip = B / 'review_delivery/CWP_Silica_CD4_IKZF2_STAT1_20261006_review_PDF_code.zip'
sha = lambda raw: hashlib.sha256(raw).hexdigest()
assert sha(source.read_bytes()) == '43fca1d334a898fa160e5136daf9f81ac153dcd0c5ec0c7508b02aa674acd52e'
assert sha(old_zip.read_bytes()) == '052af0f292c2fadf939af44ef03dc7abc5385816181b1c92f79e396aac27bd16'
title = '矽肺研究中巨噬细胞与 CD4⁺T 细胞表型关联分析'
subtitle = '基于煤工尘肺患者与小鼠硅尘模型的整合证据'
pdf_name = '矽肺研究中巨噬细胞与CD4T细胞表型关联分析_20261007_标题修订版.pdf'
zip_name = '矽肺研究中巨噬细胞与CD4T细胞表型关联分析_PDF及代码_20261007_标题修订版.zip'
dest = D / pdf_name
original = fitz.open(source)
assert len(original) == 29
doc = fitz.open(source)
page = doc[0]
region = fitz.Rect(50, 50, 547, 111)
old_title_blocks = [b for b in page.get_text('blocks') if fitz.Rect(b[:4]).intersects(region)]
assert len(old_title_blocks) == 2
page.add_redact_annot(region, fill=(1,1,1))
page.apply_redactions(images=0, graphics=0, text=0)
hei = 'C:/Windows/Fonts/simhei.ttf'
song = 'C:/Windows/Fonts/simsun.ttc'
symbol = 'C:/Windows/Fonts/times.ttf'
fh, fs, fp = fitz.Font(fontfile=hei), fitz.Font(fontfile=song), fitz.Font(fontfile=symbol)
assert fp.has_glyph(0x207a)
page.insert_font(fontname='TitleHei', fontfile=hei)
page.insert_font(fontname='TitleSong', fontfile=song)
page.insert_font(fontname='TitlePlus', fontfile=symbol)
segments = [('矽肺研究中巨噬细胞与 CD4', fh, 'TitleHei'), ('⁺', fp, 'TitlePlus'), ('T 细胞表型关联分析', fh, 'TitleHei')]
size = 18
width = sum(font.text_length(text, fontsize=size) for text,font,_ in segments)
assert width < 490
x = (page.rect.width-width)/2
for text,font,name in segments:
    page.insert_text((x,80), text, fontsize=size, fontname=name, color=(17/255,)*3)
    x += font.text_length(text, fontsize=size)
subsize=11
subwidth=fs.text_length(subtitle,fontsize=subsize)
page.insert_text(((page.rect.width-subwidth)/2,104),subtitle,fontsize=subsize,fontname='TitleSong',color=(17/255,)*3)
metadata = dict(doc.metadata)
metadata['title'] = title + '——' + subtitle
doc.set_metadata(metadata)
doc.subset_fonts()
doc.save(dest, garbage=4, deflate=True)
doc.close()
revised = fitz.open(dest)
assert len(revised) == 29
assert revised.metadata['title'] == title+'——'+subtitle
assert {k:v for k,v in original.metadata.items() if k!='title'} == {k:v for k,v in revised.metadata.items() if k!='title'}
norm=lambda s:re.sub(r'\s+','',s)
new_title_text = revised[0].get_text(clip=region)
assert norm(title) in norm(new_title_text), repr(new_title_text)
assert norm(subtitle) in norm(new_title_text)
checks=[]
for i,(old,new) in enumerate(zip(original,revised)):
    assert old.rect == new.rect
    if i==0:
        clip=fitz.Rect(0,111,old.rect.width,old.rect.height)
        ot,nt=old.get_text(clip=clip),new.get_text(clip=clip)
    else:
        ot,nt=old.get_text(),new.get_text()
    assert ot==nt, ('body text changed',i+1)
    assert re.findall(r'\d+(?:\.\d+)?',ot) == re.findall(r'\d+(?:\.\d+)?',nt)
    # Preserve every image stream and its placement, including all 14 report figures.
    def imagedata(p):
        values=[]
        for info in p.get_image_info(xrefs=True):
            item=dict(info)
            item.pop('xref',None);item.pop('number',None)
            values.append(item)
        return values
    assert imagedata(old)==imagedata(new), ('image content/placement changed',i+1)
    op=old.get_pixmap(matrix=fitz.Matrix(2,2),alpha=False)
    npix=new.get_pixmap(matrix=fitz.Matrix(2,2),alpha=False)
    oa=np.frombuffer(op.samples,dtype=np.uint8).reshape(op.height,op.width,op.n)
    na=np.frombuffer(npix.samples,dtype=np.uint8).reshape(npix.height,npix.width,npix.n)
    different=np.any(oa!=na,axis=2)
    changed_total=int(different.sum())
    if i==0:
        different[math.floor(region.y0*2):math.ceil(region.y1*2),math.floor(region.x0*2):math.ceil(region.x1*2)]=False
        new.get_pixmap(matrix=fitz.Matrix(1.5,1.5)).save(D/'QA/title_after.png')
    outside=int(different.sum())
    assert outside==0, ('rendering outside title differs',i+1,outside)
    checks.append(dict(page=i+1,body_text_identical=True,numeric_tokens_identical=True,image_content_and_positions_identical=True,pixels_changed_outside_title=outside,total_pixels_changed=changed_total))
(D/'QA/page_comparison.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
raw_pdf=dest.read_bytes()
pdf_hash=sha(raw_pdf)
with zipfile.ZipFile(old_zip) as old:
    entries={x.filename:old.read(x.filename) for x in old.infolist()}
    infos=old.infolist()
pdf_member=[n for n in entries if n.endswith('.pdf')]
assert pdf_member==['report/CWP_Silica_CD4_IKZF2_STAT1_20261006_review.pdf']
assert entries[pdf_member[0]]==source.read_bytes()
new_entries=dict(entries)
new_entries[pdf_member[0]]=raw_pdf
validation='provenance/package_validation.json'
old_validation=json.loads(entries[validation])
new_validation=dict(old_validation)
new_validation.update(pdf_sha256=pdf_hash,pdf_bytes=len(raw_pdf))
assert {k:v for k,v in new_validation.items() if k not in ['pdf_sha256','pdf_bytes']}=={k:v for k,v in old_validation.items() if k not in ['pdf_sha256','pdf_bytes']}
new_entries[validation]=json.dumps(new_validation,ensure_ascii=False,indent=2).encode('utf-8')
manifest='FILE_MANIFEST_SHA256.tsv'
lines=entries[manifest].decode('utf-8').splitlines()
assert lines[0]=='relative_path\tbytes\tsha256'
for i,line in enumerate(lines[1:],1):
    name,old_bytes,old_hash=line.split('\t')
    assert name in entries and int(old_bytes)==len(entries[name]) and old_hash==sha(entries[name]),name
    if name in [pdf_member[0],validation]:
        lines[i]=f'{name}\t{len(new_entries[name])}\t{sha(new_entries[name])}'
new_entries[manifest]=('\n'.join(lines)+'\n').encode('utf-8')
changed=[name for name in entries if entries[name]!=new_entries[name]]
assert set(changed)=={pdf_member[0],validation,manifest},changed
assert not any('264182' in n or 'silicosis_bridge' in n for n in new_entries)
code=[n for n in entries if n.startswith('code/') and Path(n).suffix.lower() in ['.py','.r','.sh','.ps1']]
assert len(code)==285,len(code)
assert all(entries[n]==new_entries[n] for n in code)
with zipfile.ZipFile(D/zip_name,'w') as z:
    for info in infos:
        z.writestr(info,new_entries[info.filename])
with zipfile.ZipFile(D/zip_name) as z:
    assert z.testzip() is None
    assert z.namelist()==list(entries)
    assert len(z.namelist())==375
    assert all(not n.startswith(('/', '\\')) and '..' not in Path(n).parts for n in z.namelist())
    assert z.read(pdf_member[0])==raw_pdf
    for line in z.read(manifest).decode('utf-8').splitlines()[1:]:
        name,count,digest=line.split('\t')
        raw=z.read(name)
        assert len(raw)==int(count) and sha(raw)==digest,name
word=B/'graphics_delivery/CD4_IKZF2_STAT1_integrated_manuscript_zh.docx'
status=dict(source_library_id='libfile_51bec083de2c81919ea20e968043ae82',source_library_version=5,
    source_pdf_sha256=sha(source.read_bytes()),source_zip_sha256=sha(old_zip.read_bytes()),
    title=title,subtitle=subtitle,pages=29,pdf=str(dest),pdf_sha256=pdf_hash,pdf_bytes=len(raw_pdf),
    zip=str(D/zip_name),zip_sha256=sha((D/zip_name).read_bytes()),zip_bytes=(D/zip_name).stat().st_size,
    title_region_points=list(region),all_29_pages_body_text_identical=True,all_29_pages_rendering_outside_title_identical=True,
    all_image_contents_and_positions_identical=True,numeric_tokens_outside_title_identical=True,
    report_figures_unchanged=14,report_tables_unchanged=18,
    code_files_unchanged=285,zip_members=375,changed_zip_members=changed,
    all_other_zip_members_byte_identical=True,zip_manifest_records_verified=len(lines)-1,
    new_validation_not_added=True,new_scientific_models=0,word_modified=False,word_sha256=sha(word.read_bytes()),
    old_files_preserved=True,library_update_status='pending')
(D/'provenance/TITLE_REVISION_STATUS.json').write_text(json.dumps(status,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(status,indent=2,ensure_ascii=True))
