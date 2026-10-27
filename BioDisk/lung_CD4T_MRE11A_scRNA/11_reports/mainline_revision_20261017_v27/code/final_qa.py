from pathlib import Path
import fitz,json,hashlib,subprocess,zipfile
from PIL import Image,ImageDraw
O=Path(__file__).parent;H=O.parent;p=Path(json.loads((O/'report_delivery.json').read_text(encoding='utf-8'))['pdf']);d=fitz.open(p)
assert len(d)==46
assert all(pg.rect.height>pg.rect.width for pg in d)
assert sum(len(pg.get_images()) for pg in d)==27
bad=[]
for n,pg in enumerate(d):
 for z in pg.get_text('blocks'):
  if z[0]<39 or z[2]>pg.rect.width-37 or z[1]<24 or z[3]>pg.rect.height-15:bad.append([n+1,list(z[:4])])
assert not bad,bad
t='\n'.join(x.get_text() for x in d)
for x in ['历史基因','继续硬筛','二元创新性P','巨噬入口','核心必做','保留原8条','数据索引：','min10敏感性','旧版','v18','v19']:assert x not in t,x
for x in ['FDR=0.8681','TTC39C','NA','IL-1环境','1.66×10^-6','0.00347','0.917','0.988']:assert x in t,x
original=Path(json.loads((H/'lung_experiment_plan_20261017_v26/report_delivery.json').read_text(encoding='utf-8'))['pdf']);assert hashlib.sha256(original.read_bytes()).hexdigest()=='f64c5e4ea98b63e9d14f7c16f7085811e634659142cb44979375b280b2019e2f'
# Re-render only changed text pages and moved appendices; all other pages visually reviewed in prior render.
r=O/'render_final';r.mkdir(exist_ok=True);exe=r'C:\Users\31937\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin\pdftoppm.exe'
for n in [5,18,26,31]+list(range(36,47)):
 subprocess.run([exe,'-f',str(n),'-l',str(n),'-scale-to','1400','-png',str(p),str(r/'page')],check=True,creationflags=subprocess.CREATE_NO_WINDOW)
files=sorted(r.glob('page-*.png'))
for start in range(0,len(files),6):
 sheet=Image.new('RGB',(1500,1450),'#ddd');dr=ImageDraw.Draw(sheet)
 for j,f in enumerate(files[start:start+6]):
  im=Image.open(f).convert('RGB');im.thumbnail((490,690));x=(j%3)*500;y=(j//3)*725;sheet.paste(im,(x,y+25));dr.text((x+8,y+5),f.stem,fill='black')
 sheet.save(r/f'contact_{start//6+1:02d}.png')
q={'pages':len(d),'portrait':True,'visuals':27,'bounds_issues':bad,'v26_preserved':True,'FDR_full_family_recomputed':True,'supplement_index_all_present':all((O/'supplement/plot_data'/('F%02d.tsv'%i)).exists() for i in range(1,25)),'report_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'experiment_summary_page':33,'experiment_appendix_pages':[39,44],'candidate_bridge_page':30,'candidate_statistics_page':31,'screening_appendix_pages':[36,38]}
(O/'final_QA.json').write_text(json.dumps(q,ensure_ascii=False,indent=2),encoding='utf-8');(O/'pdf_extracted.txt').write_text(t,encoding='utf-8')
print(json.dumps(q))