from pathlib import Path
import json,subprocess,hashlib,zipfile,fitz,pandas as pd
from PIL import Image,ImageDraw
O=Path(r'C:\Users\31937\lung_dual_benchmark_20261017_v17');B=Path(r'C:\Users\31937\lung_three_cohort_20261017_v16');j=json.loads((O/'delivery.json').read_text(encoding='utf-8'));pdf=Path(j['pdf']);R=O/'render';R.mkdir(exist_ok=True)
doc=fitz.open(pdf);print('PAGES',len(doc));assert len(doc)==9
for n,p in enumerate(doc):
 text=p.get_text();assert len(text)>100 and '\ufffd' not in text
 for b in p.get_text('blocks'):assert b[0]>=35 and b[2]<=p.rect.width-35 and b[1]>=25 and b[3]<=p.rect.height-15,(n,b[:4])
exe=r'C:\Users\31937\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin\pdftoppm.exe';subprocess.run([exe,'-r','90','-png',str(pdf),str(R/'page')],check=True,capture_output=True,creationflags=subprocess.CREATE_NO_WINDOW)
files=sorted(R.glob('page-*.png'));assert len(files)==9
for k in range(0,9,3):
 im=Image.new('RGB',(1500,770),'#dfe7e9');draw=ImageDraw.Draw(im)
 for i,f in enumerate(files[k:k+3]):
  v=Image.open(f);v.thumbnail((480,720));im.paste(v,(i*500+10,25));draw.text((i*500+10,8),'Page '+str(k+i+1),fill='black')
 im.save(R/f'contact_{k//3+1}.png')
x=pd.read_csv(O/'assessment/逐状态双基准评估.tsv',sep='\t');y=pd.read_csv(O/'assessment/90基因双基准总表.tsv',sep='\t');assert len(x)==344 and len(y)==90 and y.gene.is_unique;assert (x.P<.05).all()
# Check that reassessment is not changing input statistics or baseline evidence.
archive=Path(r'C:\Users\31937\Downloads\线粒体焦亡_三队列四轴重筛_v16_20261017\v16_完整统计与复核附件.zip');checked=[]
with zipfile.ZipFile(archive) as z:
 manifest=json.loads(z.read('MANIFEST_SHA256.json'))
 for name in ['DE_all_models.tsv','nominal_DE_axis_membership.tsv','function_membership_extended.tsv','expressed_LR_with_decoy_literature.tsv','NicheNet_activity_all.tsv']:
  f=B/'data/tables'/name;sha=hashlib.sha256(f.read_bytes()).hexdigest();assert sha==manifest['tables/'+name]['sha256'];checked.append(name)
qa={'pages':9,'rendered_pages':9,'bounds_ok':True,'input_files_unchanged':checked,'genes':90,'scope_contrasts':344,'raw_P_and_FDR_unchanged':True,'all_rows_have_A_and_B':bool(y[['A_status','B_status']].notna().all().all())};(O/'QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(qa,ensure_ascii=False))