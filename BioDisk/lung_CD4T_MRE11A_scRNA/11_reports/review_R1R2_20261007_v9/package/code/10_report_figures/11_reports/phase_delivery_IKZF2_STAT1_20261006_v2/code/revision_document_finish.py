h('参考文献',1)
refs=[
('Cassel SL, et al. The Nalp3 inflammasome is essential for the development of silicosis. PNAS. 2008. DOI: 10.1073/pnas.0803933105. PMID: 18577586.','https://pmc.ncbi.nlm.nih.gov/articles/PMC2449360/'),
('Peng Z, et al. Impaired interferon-γ signaling promotes the development of silicosis. iScience. 2022;25:104647. DOI: 10.1016/j.isci.2022.104647. PMID: 35800765.','https://pubmed.ncbi.nlm.nih.gov/35800765/'),
('You Y, et al. Crystalline silica-induced recruitment and immuno-imbalance of CD4+ tissue resident memory T cells promote silicosis progression. Communications Biology. 2024;7:971. DOI: 10.1038/s42003-024-06662-z. PMID: 39122899.','https://www.nature.com/articles/s42003-024-06662-z'),
('Zhang Y, et al. Dihydrotanshinone I Alleviates Crystalline Silica-Induced Pulmonary Inflammation by Regulation of the Th Immune Response and Inhibition of STAT1/STAT3. Mediators of Inflammation. 2019:3427053. DOI: 10.1155/2019/3427053. PMID: 31379467.','https://pmc.ncbi.nlm.nih.gov/articles/PMC6652093/'),
('Ding M, et al. Silica-exposed macrophages-secreted exosomal miR125a-5p induces Th1/Th2 and Treg/Th17 cell imbalance and promotes fibroblast transdifferentiation. Ecotoxicology and Environmental Safety. 2023;267:115647. DOI: 10.1016/j.ecoenv.2023.115647. PMID: 37918332.','https://pubmed.ncbi.nlm.nih.gov/37918332/'),
('Akimova T, et al. Helios Expression Is a Marker of T Cell Activation and Proliferation. PLoS ONE. 2011;6:e24226. DOI: 10.1371/journal.pone.0024226. PMID: 21918685.','https://pmc.ncbi.nlm.nih.gov/articles/PMC3168881/'),
('Szabo PA, et al. Distinct transcription factors control tissue adaptation and effector function in infant and adult memory T cells. Nature Immunology. 2026;27:1517–1527. DOI: 10.1038/s41590-026-02535-1. PMID: 42230957.','https://www.nature.com/articles/s41590-026-02535-1'),
('Afkarian M, et al. T-bet is a STAT1-induced regulator of IL-12R expression in naïve CD4+ T cells. Nature Immunology. 2002;3:549–557. DOI: 10.1038/ni794. PMID: 12006974.','https://pubmed.ncbi.nlm.nih.gov/12006974/'),
('Aibar S, et al. SCENIC: single-cell regulatory network inference and clustering. Nature Methods. 2017;14:1083–1086. DOI: 10.1038/nmeth.4463. PMID: 28991892.','https://pubmed.ncbi.nlm.nih.gov/28991892/'),
('Browaeys R, et al. NicheNet: modeling intercellular communication by linking ligands to target genes. Nature Methods. 2020 (online 2019). DOI: 10.1038/s41592-019-0667-5. PMID: 31819264.','https://pubmed.ncbi.nlm.nih.gov/31819264/'),
('Squair JW, et al. Confronting false discoveries in single-cell differential expression. Nature Communications. 2021;12:5692. DOI: 10.1038/s41467-021-25960-2.','https://www.nature.com/articles/s41467-021-25960-2'),
('Saravanan V, Berman GJ, Sober SJ. Application of hierarchical bootstrap to multi-level data in neuroscience. Neurons, Behavior, Data Analysis, and Theory. 2020;3(5). PMID: 33644783; PMCID: PMC7906290.','https://pmc.ncbi.nlm.nih.gov/articles/PMC7906290/')]
for i,(citation,url) in enumerate(refs,1):
 pp=p(f'[{i}] {citation}\n{url}');pp.paragraph_format.space_after=Pt(7)
 for r in pp.runs:r.font.size=Pt(9)
for st in doc.styles:
 if hasattr(st,'font'):
  st.font.name='Microsoft YaHei';rp=st._element.get_or_add_rPr();fonts=rp.find(qn('w:rFonts'))
  if fonts is None:fonts=OxmlElement('w:rFonts');rp.insert(0,fonts)
  for key in list(fonts.attrib):
   if 'Theme' in key or 'theme' in key:del fonts.attrib[key]
  for key in ['ascii','hAnsi','eastAsia','cs']:fonts.set(qn('w:'+key),'Microsoft YaHei')
for run in doc.element.iter(qn('w:r')):
 rp=run.find(qn('w:rPr'))
 if rp is None:rp=OxmlElement('w:rPr');run.insert(0,rp)
 fonts=rp.find(qn('w:rFonts'))
 if fonts is None:fonts=OxmlElement('w:rFonts');rp.insert(0,fonts)
 for key in ['ascii','hAnsi','eastAsia','cs']:fonts.set(qn('w:'+key),'Microsoft YaHei')
out=D/'CD4_IKZF2_STAT1_stage_delivery_zh.docx';doc.save(out)
(D/'document_manifest.json').write_text(json.dumps({'tables':tables,'figures':figures,'references':refs,'new_scientific_runs':0},ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'docx':str(out),'tables':len(tables),'figures':len(figures),'references':len(refs)},ensure_ascii=False))
