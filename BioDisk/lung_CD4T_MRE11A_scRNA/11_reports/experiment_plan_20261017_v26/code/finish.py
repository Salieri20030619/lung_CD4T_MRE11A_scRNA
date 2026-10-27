from pathlib import Path
import json,subprocess,time,hashlib,fitz
w=Path(r'C:\Users\31937\lung_experiment_plan_20261017_v26');m=json.loads((w/'report_delivery.json').read_text(encoding='utf-8'));p=Path(m['pdf'])
b=Path(r'C:\Users\31937\Downloads\矽肺_CD4候选研究_队列衔接修订_v25_20261017\矽肺免疫全景_CD4候选与验证_报告_v25.pdf')
d=fitz.open(p);o=fitz.open(b)
def body(page):return '\n'.join(z for z in page.get_text().splitlines() if not z.strip().isdigit())
assert body(d[40])==body(o[35])
def imgs(doc):return [hashlib.sha256(doc.extract_image(i[0])['image']).hexdigest() for pg in doc for i in pg.get_images()]
assert imgs(d)==imgs(o)
q=json.loads((w/'qa.json').read_text(encoding='utf-8'));q.update({'visual_review_pages_35_42':'passed','discussion_unchanged':True,'all_image_content_hashes_unchanged':True})
(w/'qa.json').write_text(json.dumps(q,ensure_ascii=False,indent=2),encoding='utf-8')
remote='/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/experiment_plan_20261017_v26'
def run(args):
 for n in range(4):
  x=subprocess.run(args,capture_output=True,text=True,encoding='utf-8',creationflags=subprocess.CREATE_NO_WINDOW)
  if x.returncode==0:return x.stdout
  if n==3:raise RuntimeError(x.stderr)
  time.sleep(2)
ssh=['ssh','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10']
run(ssh+['mkdir -p '+remote+'/delivery '+remote+'/code'])
run(['scp','-o','BatchMode=yes','-o','ConnectTimeout=20',str(p),'sgs10:'+remote+'/delivery/report_v26.pdf'])
run(['scp','-o','BatchMode=yes','-o','ConnectTimeout=20']+[str(x) for x in w.glob('*.py')]+['sgs10:'+remote+'/code/'])
run(['scp','-o','BatchMode=yes','-o','ConnectTimeout=20']+[str(x) for x in w.glob('*.json')]+[str(w/'experiment_chapter_source.txt'),str(w/'report_text.txt'),'sgs10:'+remote+'/'])
result=run(ssh+['sha256sum '+remote+'/delivery/report_v26.pdf'])
assert result.split()[0]==m['sha256']
complete={'status':'complete','pdf':str(p),'pages':42,'experiment_pages':[35,40],'remote':remote,'remote_sha256_verified':True,'analysis_and_figures_unchanged':True}
(w/'completion.json').write_text(json.dumps(complete,ensure_ascii=False,indent=2),encoding='utf-8')
run(['scp','-o','BatchMode=yes','-o','ConnectTimeout=20',str(w/'completion.json'),'sgs10:'+remote+'/'])
print(json.dumps(complete,ensure_ascii=False))