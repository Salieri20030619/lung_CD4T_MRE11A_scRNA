from pathlib import Path
import hashlib,json,shutil
import pymupdf
import numpy as np

authored=Path('technical_delivery/output/pdf/IKZF2_STAT1_technical_delivery_zh.pdf')
saved=Path('library_verify/IKZF2_STAT1_technical_delivery_zh.pdf')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
a=pymupdf.open(authored);b=pymupdf.open(saved)
assert len(a)==len(b)==14
pages=[]
for i in range(14):
 pa=a[i].get_pixmap(matrix=pymupdf.Matrix(1.5,1.5));pb=b[i].get_pixmap(matrix=pymupdf.Matrix(1.5,1.5))
 same=pa.width==pb.width and pa.height==pb.height and pa.samples==pb.samples
 text_same=a[i].get_text()==b[i].get_text()
 pages.append({'page':i+1,'raster_identical':same,'text_identical':text_same})
 assert same and text_same,('Saved PDF page differs',i+1)
result={'authored_sha256':sha(authored),'saved_library_pdf_sha256':sha(saved),'authored_bytes':authored.stat().st_size,'saved_library_pdf_bytes':saved.stat().st_size,'all_14_pages_raster_and_text_identical':True,'pages':pages,'action':'Use exact saved Library PDF bytes inside technical ZIP, preserving original authoring output separately.'}
a.close();b.close()
(Path('library_verify')/'pdf_revision_verification.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
shutil.copy2(authored,'library_verify/authored_before_library.pdf')
shutil.copy2(saved,authored)
shutil.copy2('library_verify/pdf_revision_verification.json','technical_delivery/provenance/pdf_revision_verification.json')
print(json.dumps(result,indent=2))
