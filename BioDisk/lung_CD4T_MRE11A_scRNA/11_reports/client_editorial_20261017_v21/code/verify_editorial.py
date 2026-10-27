from pathlib import Path
import json,re,collections,hashlib,pymupdf as fitz
O=Path(__file__).parent;B=O.parent/'lung_client_rebuild_20261017_v20';old=json.loads((B/'report_delivery.json').read_text(encoding='utf-8'));new=json.loads((O/'report_delivery.json').read_text(encoding='utf-8'));a=fitz.open(old['pdf']);b=fitz.open(new['pdf']);tx='\n'.join(p.get_text() for p in b);ot='\n'.join(p.get_text() for p in a)
terms=['历史','冻结','旧版','复评','入池','证据池','撤回','补回','复核','审核','审阅','不补造','沿用','v18','v19','v20','版本','发送','客户','本轮']
found={s:[i+1 for i,p in enumerate(b) if s in p.get_text()] for s in terms};print('RESIDUAL',{k:v for k,v in found.items() if v})
pattern=r'(?:P|FDR|rho|logFC)\s*[=≤<>]\s*[-+]?\d+(?:\.\d+)?(?:\s*[×eE]\s*[-+]?\d+(?:\^[-+]?\d+)?)?'
norm=lambda s:re.sub(r'\s+','',s)
c1=collections.Counter(map(norm,re.findall(pattern,ot)));c2=collections.Counter(map(norm,re.findall(pattern,tx)));assert c1==c2,(c1-c2,c2-c1)
assert len(a)==len(b)==31;assert sum(len(p.get_images()) for p in b)==22
bounds=[]
for i,p in enumerate(b):
 assert '\ufffd' not in p.get_text()
 for bb in p.get_text('blocks'):
  if bb[0]<39 or bb[2]>p.rect.width-37 or bb[1]<24 or bb[3]>p.rect.height-15:bounds.append((i+1,bb[:4]))
assert not bounds
R=O/'render';R.mkdir(exist_ok=True)
import subprocess
subprocess.run([r'C:\Users\31937\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin\pdftoppm.exe','-scale-to','1300','-png',new['pdf'],str(R/'page')],check=True,creationflags=subprocess.CREATE_NO_WINDOW)
from PIL import Image,ImageDraw
files=sorted(R.glob('page-*.png'))
for j in range(0,len(files),8):
 sheet=Image.new('RGB',(1400,500*((min(8,len(files)-j)+1)//2)),'#dce4e8');dr=ImageDraw.Draw(sheet)
 for k,f in enumerate(files[j:j+8]):
  im=Image.open(f).convert('RGB');im.thumbnail((690,470));x=k%2*700;y=k//2*500;sheet.paste(im,(x,y+20));dr.text((x+4,y+4),f.stem,fill='black')
 sheet.save(R/('contact_'+str(j//8+1)+'.png'))
qa={'pages':31,'figures':22,'P_FDR_rho_logFC_tokens_identical':True,'statistic_tokens_checked':sum(c1.values()),'residual_terms':{k:v for k,v in found.items() if v},'bounds_flags':bounds,'edited_figures_data_parity':'F01/F12/F17/F20 checked frame-identical before rendering; other18 unchanged source assets','source_v20_preserved':hashlib.sha256(Path(old['pdf']).read_bytes()).hexdigest()==old['sha256']};(O/'editorial_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(qa,ensure_ascii=False))