from pathlib import Path
from collections import Counter
import json, re, hashlib, zipfile, copy
from lxml import html
import pymupdf as fitz

B=Path(__file__).resolve().parent;D=B/'logic_revision_20261007'
sha=lambda b:hashlib.sha256(b).hexdigest()
tree=html.fromstring((D/'provenance/logic_render_source.html').read_text(encoding='utf-8'))
nodes=tree.xpath('//body/*')
old=json.loads((D/'provenance/original_blocks.json').read_text(encoding='utf-8'))
text=[e.text_content() for e in nodes]
figures=[re.match(r'图(S?\d+)',e.xpath('./p[@class="Caption"]')[0].text_content()).group(1) for e in tree.xpath('//div[@class="figure"]')]
tables=[re.match(r'表(S?\d+)',text[i]).group(1) for i in range(len(nodes)) if nodes[i].get('class')=='Caption' and i+1<len(nodes) and nodes[i+1].tag=='table']
assert figures==[str(i) for i in range(1,11)]+['S1','S2','S3','S4'],figures
assert tables==[str(i) for i in range(1,5)]+['S'+str(i) for i in range(1,15)],tables
fig_refs=set(re.findall(r'图\s*(S?\d+)', '\n'.join(text)))
tab_refs=set(re.findall(r'表\s*(S?\d+)', '\n'.join(text)))
assert fig_refs<=set(figures) and tab_refs<=set(tables)
assert len(nodes[232:])==15 and all(re.match(r'\['+str(i+1)+r'\]',n.text_content()) for i,n in enumerate(nodes[232:]))
refs=set(int(x) for s in text[:232] for x in re.findall(r'\[(\d+)\]',s))
assert refs<=set(range(1,16))
headers={re.match(r'(\d+\.\d+)\s',e.text_content()).group(1) for e in nodes if e.tag in ['h1','h2','h3'] and re.match(r'(\d+\.\d+)\s',e.text_content())}
section_refs=set(re.findall(r'(?:方法|结果|讨论)\s*(\d+\.\d+)', '\n'.join(text)))
assert section_refs<=headers,(section_refs,headers)
numeric=re.compile(r'[−+]?\d+(?:\.\d+)?(?:[eE][−+-]?\d+)?')
table_checks=[]
for i,e in enumerate(nodes):
    if e.tag!='table':continue
    actual=[[c.text_content() for c in row.xpath('./th|./td')] for row in e.xpath('.//tr')]
    assert len(actual)==len(old[i]['text'])
    for ri,(r1,r2) in enumerate(zip(old[i]['text'],actual)):
        assert len(r1)==len(r2)
        for ci,(a,b) in enumerate(zip(r1,r2)):
            assert Counter(numeric.findall(a))==Counter(numeric.findall(b)),(i,ri,ci,a,b)
    table_checks.append({'block':i,'rows':len(actual),'numeric_cells_identical':True})
decimal=json.loads((D/'QA/decimal_diff_review.json').read_text(encoding='utf-8'))
assert decimal==[{'block':17,'cell':None,'removed':{},'added':{'3.4':1}}, {'block':42,'cell':None,'removed':{},'added':{'3.3':1}}]
assert json.loads((D/'QA/history_term_scan.json').read_text(encoding='utf-8'))==[]
assert all(term not in '\n'.join(text) for term in ['GSE264182','silicosis_bridge'])
word=B/'graphics_delivery/CD4_IKZF2_STAT1_integrated_manuscript_zh.docx'
assert sha(word.read_bytes())=='1bd05f6e28a5c60da9a012b5397f237e0999919d43e0890b893074314882476c'
pdf=next(D.glob('*最终逻辑核查版.pdf'));pbytes=pdf.read_bytes();doc=fitz.open(pdf)
parity=json.loads((D/'QA/PDF_content_parity.json').read_text(encoding='utf-8'))
assert len(doc)==parity['pages'] and not any(r['outside_main_margins'] for r in parity['page_records'])
checks=[
 {'topic':'CD4 cohort and source confounding','blocks':[15,16,115,116],'finding':'RNA19=10 control/9 CWP; source counts 2/0,8/2,0/7; donor unit; source-confounding and CWP applicability consolidated in discussion. AUC and subtype rosters remain distinct.'},
 {'topic':'Models and estimability','blocks':[24,25,26,167,187],'finding':'Outcome-specific transformations stated; AUC reference disease BH and association P output settings attributed to those methods only; missing is not negative. Source code stats.py and align_compute.py inspected; no model rerun.'},
 {'topic':'IL7R excess and Bootstrap','blocks':[38,40,41,82,91,140,141,142],'finding':'88-gene matched-background endpoint separated from PRF1/IFNG/NKG7 endpoint. Fixed16, RNA19 and matchedRNA17 stated. Bootstrap CIs and paired differences interpreted separately.'},
 {'topic':'IL15 evidence layers','blocks':[29,31,32,69,119,121,122,125],'finding':'Receiver-common AUPR separated from sender LR/complex coverage and same-donor tests. MacAM one qualifying case, MonoMac/MacWide two, strict four-gene cases zero; incomplete activity target set can have independent program correlation. Source stage2.R inspected.'},
 {'topic':'GZMA source and IKZF2 coverage','blocks':[81,127,128,130,131],'finding':'n12 comprises G227 7 controls/2 cases and G288 3 cases; G135 absent. 29/2554 detected cells and 31 UMI retained.'},
 {'topic':'Composition sensitivity','blocks':[36,106,107,134,135,137,138],'finding':'Fixed Core17 composition model and background-filtered15 analysis are distinct; composition df8 vs primary df11. Nominal and full-family BH results retained.'},
 {'topic':'Independent candidates and species','blocks':[78,84,86,87,133,143,149,153,154,155],'finding':'IKZF2/STAT1 interpreted independently. Mouse AUC differences nominal with uncomputed BH; mouse unit is library; GSE183682 single-library time points direction only.'},
 {'topic':'Pathways and scores','blocks':[44,45,49,50,51,52,98,104,145,146,147,148,150,151,222],'finding':'GSEA, GSVA and UCell meanings distinct; apoptosis regulation is not overall apoptosis. Top5 curves cross-reference corrected; PyroResp/Inflam/ApoExec caption definitions checked against module members.'},
 {'topic':'UMAP and manuscript structure','blocks':[54,61,65],'finding':'Full-lung and CD4 embeddings and lineage labels distinct; source attribution remains precise. All figure/table/section citations resolve; bibliography 1-15 retained.'}
]
qa={'pages_reviewed':len(doc),'all_pages_rendered':True,'visual_review':'All 30 pages reviewed in four contact sheets; dense methods/table pages inspected separately. No clipped blocks or isolated captions; continued tables repeat headers.',
 'figures':figures,'tables':tables,'figure_references':sorted(fig_refs),'table_references':sorted(tab_refs),'section_references':sorted(section_refs),'unresolved_references':[],
 'table_numeric_checks':table_checks,'decimal_differences':'Only two added section references (3.4 and 3.3); scientific decimal values unchanged',
 'word_sha256':sha(word.read_bytes()),'word_unchanged':True,'new_scientific_models':0,'logic_checks':checks}
(D/'QA/final_logic_audit.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')

baseline=next((B/'title_revision_20261007').glob('*.zip'))
assert sha(baseline.read_bytes())=='87e789a2351eb70a1f5c887941290b6a2e8f694fd8d980a580ed602688e5ae97'
oldzip=zipfile.ZipFile(baseline)
pdfpath='report/CWP_Silica_CD4_IKZF2_STAT1_20261006_review.pdf'
validationpath='provenance/package_validation.json'
manifestpath='FILE_MANIFEST_SHA256.tsv'
validation=json.loads(oldzip.read(validationpath));validation.update(pdf_pages=len(doc),pdf_sha256=sha(pbytes),pdf_bytes=len(pbytes),all_pdf_pages_reviewed=len(doc))
replacement={pdfpath:pbytes,validationpath:json.dumps(validation,ensure_ascii=False,indent=2).encode('utf-8')}
manifest=oldzip.read(manifestpath).decode('utf-8').splitlines()
newmanifest=[manifest[0]]
for line in manifest[1:]:
    path,size,digest=line.split('\t')
    if path in replacement:size=str(len(replacement[path]));digest=sha(replacement[path])
    newmanifest.append('\t'.join([path,size,digest]))
replacement[manifestpath]=('\n'.join(newmanifest)+'\n').encode('utf-8')
out=D/'矽肺研究中巨噬细胞与CD4T细胞表型关联分析_PDF及代码_20261007_最终逻辑核查版.zip'
with zipfile.ZipFile(out,'w') as z:
    for info in oldzip.infolist():z.writestr(copy.copy(info),replacement.get(info.filename,oldzip.read(info.filename)))
z=zipfile.ZipFile(out)
assert z.testzip() is None and z.namelist()==oldzip.namelist() and len(z.namelist())==375
changed=[n for n in z.namelist() if z.read(n)!=oldzip.read(n)]
assert set(changed)==set(replacement),changed
scripts=[n for n in z.namelist() if n.startswith('code/') and Path(n).suffix.lower() in ['.py','.r','.sh','.ps1']]
assert len(scripts)==285 and all(z.read(n)==oldzip.read(n) for n in scripts)
assert z.read(pdfpath)==pbytes
for line in z.read(manifestpath).decode('utf-8').splitlines()[1:]:
    n,size,digest=line.split('\t');b=z.read(n);assert len(b)==int(size) and sha(b)==digest,n
assert not any('GSE264182' in n or 'silicosis_bridge' in n for n in z.namelist())
packageqa={'archive':out.name,'bytes':out.stat().st_size,'sha256':sha(out.read_bytes()),'members':len(z.namelist()),'manifest_entries_verified':len(manifest)-1,'script_count':len(scripts),'all_scripts_byte_identical':True,'changed_members':changed,'all_other_members_byte_identical':True,'archive_integrity':'passed','embedded_pdf_equals_standalone':True,'max_member_path_length':max(map(len,z.namelist()))}
(D/'QA/ZIP_integrity.json').write_text(json.dumps(packageqa,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'pdf':parity['pdf'],'pages':len(doc),'zip':packageqa},ensure_ascii=True))
