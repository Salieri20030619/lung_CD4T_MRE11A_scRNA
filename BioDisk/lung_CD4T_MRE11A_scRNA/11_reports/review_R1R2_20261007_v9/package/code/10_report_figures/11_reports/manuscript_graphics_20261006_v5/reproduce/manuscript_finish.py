doc.add_page_break();h('参考文献',1)
refs=json.loads((B/'revision_delivery/document_manifest.json').read_text(encoding='utf8'))['references']
refs += [
['Subramanian A, et al. Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. PNAS. 2005;102:15545–15550. DOI: 10.1073/pnas.0506580102.','https://www.pnas.org/doi/10.1073/pnas.0506580102'],
['Hänzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis for microarray and RNA-Seq data. BMC Bioinformatics. 2013;14:7. DOI: 10.1186/1471-2105-14-7.','https://link.springer.com/article/10.1186/1471-2105-14-7'],
['Andreatta M, Carmona SJ. UCell: robust and scalable single-cell gene signature scoring. Computational and Structural Biotechnology Journal. 2021;19:3796–3798. DOI: 10.1016/j.csbj.2021.06.043. 方法与引用信息见作者维护资源。','https://github.com/carmonalab/UCell']]
for i,(citation,url) in enumerate(refs,1):
 pp=p(f'[{i}] {citation}\n{url}');pp.paragraph_format.first_line_indent=Pt(0);pp.paragraph_format.space_after=Pt(7);pp.paragraph_format.alignment=WD_ALIGN_PARAGRAPH.LEFT
 for r in pp.runs:r.font.size=Pt(9)
# Fix font inheritance at run level, preserving meaningful bold, italics and scientific names.
for run in doc.element.iter(qn('w:r')):
 rp=run.find(qn('w:rPr'))
 if rp is None:rp=OxmlElement('w:rPr');run.insert(0,rp)
 fonts=rp.find(qn('w:rFonts'))
 if fonts is None:fonts=OxmlElement('w:rFonts');rp.insert(0,fonts)
 for key in ['ascii','hAnsi','cs']:fonts.set(qn('w:'+key),'Times New Roman')
 for key in list(fonts.attrib):
  if 'Theme' in key or 'theme' in key:del fonts.attrib[key]
out=D/'CD4_IKZF2_STAT1_integrated_manuscript_zh.docx';doc.save(out)
(D/'document_manifest.json').write_text(json.dumps({'tables':tables,'figures':figures,'references':refs,'new_scientific_runs':0,'reference_style_only':True},ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({'docx':str(out),'tables':len(tables),'figures':len(figures),'references':len(refs)},ensure_ascii=True))
