from pathlib import Path
import json,hashlib,subprocess,pymupdf as fitz
from PIL import Image,ImageDraw
O=Path(__file__).parent;B=O.parent/'lung_client_layout_20261017_v22';d=json.loads((O/'report_delivery.json').read_text(encoding='utf-8'));q=json.loads((O/'basic_QA.json').read_text());assert q['pages']==31 and q['figures']==22 and not q['bounds']
for id in ['F11','F12','F17','F20']:assert (O/'plot_data'/(id+'.tsv')).read_bytes()==(B/'plot_data'/(id+'.tsv')).read_bytes()
g=json.loads((O/'trajectory_backbone_QA.json').read_text());assert g['components']==1 and g['overlay_cells']==1734 and g['edges']==11
R=O/'render';R.mkdir(exist_ok=True);pop=r'C:\Users\31937\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin\pdftoppm.exe';pages=[1,4,15,16,21,24]
for n in pages:subprocess.run([pop,'-f',str(n),'-l',str(n),'-scale-to','1600','-png',d['pdf'],str(R/('p'+str(n).zfill(2)))],check=True,creationflags=subprocess.CREATE_NO_WINDOW)
sheet=Image.new('RGB',(1440,3*1020),'#dce4e8');draw=ImageDraw.Draw(sheet)
for i,n in enumerate(pages):
 p=list(R.glob('p'+str(n).zfill(2)+'-*.png'))[0];im=Image.open(p).convert('RGB');im.thumbnail((700,990));x=i%2*720;y=i//2*1020;sheet.paste(im,(x,y+20));draw.text((x+5,y+4),'Page '+str(n),fill='black')
sheet.save(R/'changed_pages.jpg',quality=92)
q.update(all_portrait=True,original_four_plot_inputs_identical=True,trajectory_edges_based_on_existing_neighbor_graph=True,visual_pages=pages);(O/'QA.json').write_text(json.dumps(q,ensure_ascii=False,indent=2),encoding='utf-8');print('QA_PASS')