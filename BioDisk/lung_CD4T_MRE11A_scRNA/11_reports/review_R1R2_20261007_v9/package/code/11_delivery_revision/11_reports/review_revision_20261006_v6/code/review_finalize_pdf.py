from pathlib import Path
import pymupdf as fitz,json,re
from docx import Document
B=Path(__file__).parent;D=B/'review_delivery'
pdf=fitz.open(D/'review_pdf_candidate.pdf');doc=Document(D/'internal/revised_source.docx')
chunks=[p.text for p in doc.paragraphs if p.text.strip()]+[c.text for t in doc.tables for row in t.rows for c in row.cells]
normalize=lambda s:re.sub(r'\s+','',s).replace('\u00ad','')
full=normalize('\n'.join(p.get_text() for p in pdf))
missing=[s for s in chunks if normalize(s) not in full]
assert not missing,[(s[:120]) for s in missing]
for i,p in enumerate(pdf):
 p.insert_textbox(fitz.Rect(200,p.rect.height-31,p.rect.width-200,p.rect.height-17),str(i+1),fontsize=9,fontname='tiro',align=1,color=(.1,.1,.1))
pdf.set_metadata({'title':'煤工尘肺患者与小鼠硅尘模型：CD4 T细胞状态、巨噬信号与IKZF2、STAT1关联','author':'','subject':'Review revision: endpoint definitions, sender evidence, source coverage and existing composition sensitivity','creator':'Validated DOCX content rendered with local headless Edge; pagination verified','keywords':'CD4;UMAP;GSEA;IKZF2;STAT1'})
pdf.save(D/'review_report.pdf',garbage=4,deflate=True)
result={'renderer':'Local headless Microsoft Edge from exact DOCX paragraphs, tables and embedded images','reason':'Reuse validated local HTML/Edge rendering workflow; no Word export in this revision','pages':len(pdf),'docx_text_blocks_verified':len(chunks),'missing_text_blocks':missing,'page_numbers_added':True}
(D/'QA/PDF_content_parity.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8');print(json.dumps(result,ensure_ascii=True))
