from pathlib import Path
import json,hashlib,subprocess,pymupdf as fitz,pandas as pd,numpy as np
from PIL import Image,ImageDraw
O=Path(__file__).parent;B=O.parent/'lung_client_rebuild_20261017_v20';d=json.loads((O/'report_delivery.json').read_text(encoding='utf-8'));doc=fitz.open(d['pdf']);pagechecks=[]
for i,p in enumerate(doc):
 assert p.rect.height>p.rect.width,(i+1,'landscape')
 assert '\ufffd' not in p.get_text()
 bad=[]
 for b in p.get_text('blocks'):
  if b[0]<39 or b[2]>p.rect.width-37 or b[1]<24 or b[3]>p.rect.height-15:bad.append(b[:4])
 pagechecks.append(dict(page=i+1,width=p.rect.width,height=p.rect.height,images=len(p.get_images()),text_chars=len(p.get_text()),bounds=bad))
assert sum(x['images'] for x in pagechecks)==22
parity=[]
for p in (O/'plot_data').glob('F*.tsv'):
 if p.stem not in ['F09','F10','F11','F12']:
  assert p.read_bytes()==(B/'plot_data'/p.name).read_bytes();parity.append(p.stem)
q=pd.read_csv(O/'plot_data/GSEA_curve_QA.tsv',sep='\t');assert len(q)==10 and q.error.max()<1e-8
all=pd.read_csv(B/'server_tables/GSEA_all.tsv',sep='\t');tops=pd.read_csv(O/'plot_data/GSEA_top5_selected.tsv',sep='\t');merged=tops.merge(all,left_on=['scope','ID'],right_on=['scope','pathway'],suffixes=('_shown','_original'),validate='one_to_one');assert len(merged)==10
for a,b in [('pvalue','pval'),('p.adjust','padj'),('NES_shown','NES_original'),('enrichmentScore','ES')]:assert np.allclose(merged[a],merged[b],atol=1e-12)
ch=pd.read_csv(O/'plot_data/F11.tsv',sep='\t');old=pd.read_csv(B/'plot_data/F11.tsv',sep='\t');assert len(ch)==12;pd.testing.assert_frame_equal(ch.sort_values(['sender','scope','group']).reset_index(drop=True),old[old.scope.isin(['Naive','ActCD4','Treg'])].sort_values(['sender','scope','group']).reset_index(drop=True))
assert len(pd.read_csv(O/'plot_data/F12.tsv',sep='\t'))==6
R=O/'render';R.mkdir(exist_ok=True)
subprocess.run([r'C:\Users\31937\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin\pdftoppm.exe','-scale-to','1500','-png',d['pdf'],str(R/'page')],check=True,creationflags=subprocess.CREATE_NO_WINDOW)
files=sorted(R.glob('page-*.png'))
for j in range(0,len(files),8):
 sheet=Image.new('RGB',(1440,1020*((min(8,len(files)-j)+1)//2)),'#dce4e8');dr=ImageDraw.Draw(sheet)
 for k,f in enumerate(files[j:j+8]):
  im=Image.open(f).convert('RGB');im.thumbnail((700,990));x=k%2*720;y=k//2*1020;sheet.paste(im,(x,y+20));dr.text((x+4,y+4),f.stem,fill='black')
 sheet.save(R/('contact_'+str(j//8+1)+'.jpg'),quality=90)
qa={'pages':len(doc),'all_portrait':True,'figure_images':22,'page_checks':pagechecks,'unchanged_plot_inputs':parity,'GSEA_curves_reconstructed':10,'GSEA_max_ES_error':float(q.error.max()),'GSEA_NES_P_FDR_match_original':True,'chord_counts_match_original':True,'chord_receptor_edges':6}
(O/'layout_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8');print('PAGES',len(doc),'BOUNDS',sum(len(x['bounds']) for x in pagechecks),'FIGURES',22);print([(x['page'],x['images'],x['text_chars']) for x in pagechecks])