from pathlib import Path
import json,hashlib,shutil,zipfile,subprocess,time,re
import fitz,pandas as pd
O=Path(__file__).parent;H=O.parent;B=H/'lung_readability_revision_20261017_v28';m=json.loads((O/'report_delivery.json').read_text(encoding='utf-8'));p=Path(m['pdf']);d=fitz.open(p);oldp=Path(json.loads((B/'report_delivery.json').read_text(encoding='utf-8'))['pdf']);old=fitz.open(oldp)
assert len(d)==52;assert all(pg.rect.height>pg.rect.width for pg in d)
assert all(d[i].get_text()==old[i].get_text() for i in list(range(0,33))+[34,35] if i not in [1,28,31])
def imgs(doc):return [hashlib.sha256(doc.extract_image(x[0])['image']).hexdigest() for pg in doc for x in pg.get_images()]
assert imgs(d)==imgs(old) and len(imgs(d))==27
assert hashlib.sha256(oldp.read_bytes()).hexdigest()==json.loads((O/'change_scope.json').read_text(encoding='utf-8'))['base_sha256']
text='\n'.join(pg.get_text() for pg in d);assert '待排' not in text;assert '[10]' in d[28].get_text() and '[10]' in d[31].get_text()
for i,pg in enumerate(d):
 for b in pg.get_text('blocks'):assert b[0]>=39 and b[2]<=pg.rect.width-37 and b[1]>=24 and b[3]<=pg.rect.height-15,(i+1,b[:4])
reply='''v29 文献引用修订简要说明
对应v28文献引用核验意见；以下页码为v29。

1. 文献[3]仅加来源注记：第34页说明Results与Discussion关于ENO1遗传降低效应方向不一致；第44页实验路径回指该注记。不替作者认定笔误，不修改候选排序或实验预期方向。
2. 第34页将IL1R2、TTC39C、SKIL各项论断分别对应[1]—[7]；[4]明确为矽尘暴露小鼠中性粒细胞的Il1r2转录表达。
3. 第29、32页补[10]正文回引，明确其为H2-DMb1条目，不据此覆盖H2-DMb2直接功能。
4. 第43页分开“巨噬原研究的LDH/IL-1β结果”和“混合上清不能归属T细胞”的方案推论；第47页区分死亡机制文献与方案应用。
5. 第37—38页统一13项书目：补作者、卷期/文章号及DOI；[9][10]使用正式网页/数据库标题、标识和访问日期。
6. 第52页附录E列出检测试剂的厂家、货号、应用、物种、验证材料和官方链接，区分检测资源与原代T细胞功能验证。

完成相关页面渲染及引文/页码核对。原27幅图示、统计结果和候选顺序保持不变，v28保留。本文献修订不代表13项来源全部原始数据或全文均经独立审计。
'''
rpath=p.parent/'v29_文献引用修订简要回复_20261017.txt';rpath.write_text(reply,encoding='utf-8');assert rpath.read_text(encoding='utf-8')==reply
S=O/'supplement';shutil.copytree(B/'supplement',S,dirs_exist_ok=True)
resources=json.loads((O/'evidence/reagent_resources.json').read_text(encoding='utf-8'));pd.DataFrame(resources).assign(checked_date='2026-10-17').to_csv(S/'tables/report_tables/experimental_feasibility.tsv',sep='\t',index=False)
shutil.copy2(O/'references.json',S/'evidence/references_v29.json');shutil.copy2(O/'evidence/reagent_resources.json',S/'evidence/reagent_resources.json')
mp=pd.read_csv(S/'figure_mapping.tsv',sep='\t');pages={}
for i,pg in enumerate(d):
 for f in re.findall(r'(?m)^图(\d+)。',pg.get_text()):pages[int(f)]=i+1
assert sorted(pages)==list(range(1,25));mp['pdf_page']=mp.figure_number.map(pages);mp.to_csv(S/'figure_mapping.tsv',sep='\t',index=False)
readme=(S/'README.txt').read_text(encoding='utf-8').replace('v28','v29')+'\n文献及试剂：references_v29.json为完整书目；experimental_feasibility.tsv已改为厂家/货号/用途/物种/验证材料分列，报告附录E可直接核对。文献[3]按用户要求仅作原文方向不一致注记，不解决原作者矛盾。\n';(S/'README.txt').write_text(readme,encoding='utf-8')
inv=[{'path':str(f.relative_to(S)),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in S.rglob('*') if f.is_file() and f.name!='file_inventory.json'];(S/'file_inventory.json').write_text(json.dumps(inv,ensure_ascii=False,indent=2),encoding='utf-8')
zpath=O/'v29_补充统计与图值.zip'
with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for f in S.rglob('*'):
  if f.is_file():z.write(f,str(f.relative_to(S)))
with zipfile.ZipFile(zpath) as z:assert z.testzip() is None
qa={'pages':52,'portrait':True,'visuals':27,'all_images_identical_to_v28':True,'unchanged_analysis_pages_text_identical':True,'bounds_issues':[],'visual_review':'changed pages and shifted appendices rendered and checked; final p34 rechecked','ref10_cited_pages':[29,32],'references_pages':[37,38],'reagent_appendix_page':52,'v28_preserved':True,'ref3_annotation_only':True,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()};(O/'final_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
remote='/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/citation_revision_20261017_v29'
def run(args):
 for n in range(4):
  z=subprocess.run(args,capture_output=True,text=True,encoding='utf-8',creationflags=subprocess.CREATE_NO_WINDOW)
  if z.returncode==0:return z.stdout
  if n==3:raise RuntimeError(z.stderr)
  time.sleep(2)
ssh=['ssh','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10'];scp=['scp','-o','BatchMode=yes','-o','ConnectTimeout=20'];run(ssh+['mkdir -p '+remote+'/delivery '+remote+'/code'])
for f,dest in [(p,'delivery/report_v29.pdf'),(rpath,'delivery/reviewer_reply.txt'),(zpath,'supplement_v29.zip')]:run(scp+[str(f),'sgs10:'+remote+'/'+dest])
run(scp+[str(f) for f in O.glob('*.py')]+['sgs10:'+remote+'/code/']);run(scp+[str(f) for f in O.glob('*.json')]+['sgs10:'+remote+'/'])
hashes=run(ssh+['sha256sum '+remote+'/delivery/report_v29.pdf '+remote+'/delivery/reviewer_reply.txt '+remote+'/supplement_v29.zip'])
for f in [p,rpath,zpath]:assert hashlib.sha256(f.read_bytes()).hexdigest() in hashes
c={'status':'complete','pdf':str(p),'reply':str(rpath),'supplement':str(zpath),'pages':52,'remote':remote,'checksums_verified':True};(O/'completion.json').write_text(json.dumps(c,ensure_ascii=False,indent=2),encoding='utf-8');run(scp+[str(O/'completion.json'),'sgs10:'+remote+'/']);print(json.dumps(c,ensure_ascii=False))