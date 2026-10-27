from pathlib import Path
import subprocess,json,hashlib,zipfile,shutil,datetime,pandas as pd,fitz
O=Path(r'C:\Users\31937\lung_three_cohort_20261017_v16');D=Path(json.loads((O/'report_output.json').read_text(encoding='utf-8'))['delivery']);root='/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/three_cohort_four_axis_20261017_v16/'
flags={'creationflags':subprocess.CREATE_NO_WINDOW,'capture_output':True}
for folder in ['logs','code']:
 subprocess.run(['scp','-q','-r','sgs10:'+root+folder,str(O/'data')],check=True,**flags)
for name in ['legacy_v13_evidence.zip']:
 subprocess.run(['scp','-q','sgs10:'+root+'resources/'+name,str(O/'data/resources'/name)],check=True,**flags)
for name in ['*.samples.tsv','receiver_background.tsv']:
 subprocess.run(['scp','-q','sgs10:'+root+'inputs/'+name,str(O/'data/inputs')],check=True,**flags)
# The public resource bytes must match the download manifest.
manifest=json.loads((O/'data/resources/download_manifest.json').read_text(encoding='utf-8'))
original_dir=O/'data/resources/original_downloads';original_dir.mkdir(exist_ok=True)
analysis_manifest={r['file']:r for r in json.loads((O/'data/resources/manifest.json').read_text(encoding='utf-8'))}
resource_integrity=[]
for r in manifest:
 raw=O/'resources'/r['file'];analysis=O/'data/resources'/r['file'];shutil.copy2(raw,original_dir/r['file'])
 assert hashlib.sha256(raw.read_bytes()).hexdigest()==r['sha256']
 assert hashlib.sha256(analysis.read_bytes()).hexdigest()==analysis_manifest[r['file']]['sha256']
 same=raw.read_bytes()==analysis.read_bytes()
 if not same:assert raw.suffix=='.json' and json.loads(raw.read_text(encoding='utf-8'))==json.loads(analysis.read_text(encoding='utf-8'))
 resource_integrity.append({'raw_path':'original_downloads/'+r['file'],'analysis_path':r['file'],'raw_sha256':r['sha256'],'analysis_sha256':analysis_manifest[r['file']]['sha256'],'identical_bytes':same,'identical_JSON_content':True if raw.suffix=='.json' else None,'note':'Reactome analysis copy is indented JSON; raw download retained; no biological content difference' if not same else 'byte-identical'})
(O/'data/resources/resource_manifest_resolved.json').write_text(json.dumps(resource_integrity,ensure_ascii=False,indent=2),encoding='utf-8')
qa=json.loads((O/'pdf_qa.json').read_text(encoding='utf-8'));assert qa['pages']==18 and qa['rendered_pages']==18 and not qa['out_of_bounds']
pdf=next(D.glob('*.pdf'));doc=fitz.open(pdf);txt='\n'.join(p.get_text() for p in doc)
for g in ['IL1R2','DDIT4','SKIL','FOXP1','GSN']:assert g in txt
checks=pd.read_csv(O/'data/tables/QA_statistical_integrity.tsv',sep='\t');assert (checks.status=='PASS').all() and len(checks)==29
finalchecks=pd.read_csv(O/'final_tables/交集逐项复核.tsv',sep='\t');assert len(finalchecks)==155 and (finalchecks.status=='PASS').all()
summary=json.loads((O/'final_tables/summary.json').read_text(encoding='utf-8'));assert summary['novelty_queries']==165
readme='''v16 三队列四轴重筛复核附件 | 2026-10-17

先读18页PDF，再读final_tables/优先验证候选_5基因.tsv和四轴候选总表_55基因.tsv。
结论：三项生物学交集55基因，99基因-状态组合，155项对比。创新性按明确问题分层，优先IL1R2、DDIT4、SKIL、FOXP1，GSN备选。五基因是科学优先级，不是创新性P值或全球首次证明。

统计：原始P<0.05入候选；|log2FC|>=0.5为效应分层；完整FDR按预设家族保存。两个疾病数据集10独立样本采用分层模型；第三队列12样本为独立时间分析。不要将22只小鼠混作同一个病例对照比较。

文件目录：
- final_tables：55基因四轴表、5基因优先表、155项对比、42历史基因回溯、165项创新性查询、命中题名/摘要、交集逐项检查。
- tables：全基因DE、全部BH检验家族、样本/模型、GO+Reactome功能、文献功能补充、表达兼容LR、NicheNet活性/靶排名、完整通讯路径、分样本QC、全肺来源比较。
- inputs：仅复核所需样本名单、接收端背景基因和历史候选提取记录；不含细胞原始矩阵。
- resources：数据库快照/哈希、GO版本（原始下载在original_downloads；分析JSON仅缩进格式不同，resource_manifest_resolved记录两个SHA和完整JSON等价检查）、历史v13六个完整S17表压缩包（归档名为ASCII，原文件路径和逐文件哈希保留，不改变原始内容）。
- literature_success：最终成功检索原始JSON；*.previous.json为别名扩充前记录，最终结果以同名不带previous的JSON及165项汇总为准。题名/摘要命中有限制，零命中不证明新颖性；预印本独立标记。
- code：服务器分析代码。主要运行顺序：export.py -> resources.R -> de.R -> axes.py -> literature_extensions.py -> nichenet.R -> integrate_fast.py -> qc_sources.py -> qa.py。资源需先按download_manifest获取。literature_extensions是明确文献补充，不调整差异统计。
- local_code：统一文献检索、分层表、报告和渲染/打包代码。
- logs：实际执行环境及完成标记；README和summary说明最终验收。

全流程原始输入及edgeR拟合对象在服务器：
/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/three_cohort_four_axis_20261017_v16
旧v13-v15保持原样。本附件用于内部复审，不是已外审通过的客户版本。
文本UTF-8无BOM；TSV为审计用精确表，原始数值不按PDF显示精度舍入。
'''
(O/'README_review.txt').write_text(readme,encoding='utf-8')
files={}
for folder in ['tables','resources','logs','code','inputs']:
 for f in (O/'data'/folder).rglob('*'):
  if f.is_file():files[folder+'/'+f.relative_to(O/'data'/folder).as_posix()]=f
for folder in ['final_tables','literature_success']:
 for f in (O/folder).rglob('*'):
  if f.is_file():files[folder+'/'+f.relative_to(O/folder).as_posix()]=f
for n in ['download.py','freeze.py','search_novelty.py','finalize_tables.py','report.py','render_check.py','package.py']:
 files['local_code/'+n]=O/n
files['protocol_frozen.json']=O/'data/protocol_frozen.json';files['README_review.txt']=O/'README_review.txt';files['pdf_visual_QA.json']=O/'pdf_qa.json';files['resources/threshold_reference_20261017.json']=Path(r'C:\Users\31937\lung_threshold_cohorts_20261017\threshold_reference_20261017.json')
hashes={n:{'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'bytes':f.stat().st_size} for n,f in sorted(files.items())}
zip_path=D/'v16_完整统计与复核附件.zip'
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for n,f in files.items():z.write(f,n)
 z.writestr('MANIFEST_SHA256.json',json.dumps(hashes,ensure_ascii=False,indent=2).encode('utf-8'))
with zipfile.ZipFile(zip_path) as z:
 assert z.testzip() is None
 for n,rec in hashes.items():assert hashlib.sha256(z.read(n)).hexdigest()==rec['sha256']
completion={'version':'v16','completed_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'analysis_complete':True,'report_pages':18,'biological_intersection_genes':55,'priority_genes':summary['priority_genes'],'biological_contrasts_verified':155,'BH_families_verified':17,'pseudobulk_inputs_verified':12,'novelty_queries_success':165,'visual_QA':'all pages rendered and inspected; final appendix pagination rechecked','outputs':{f.name:{'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'bytes':f.stat().st_size} for f in [pdf,zip_path]},'limitations':'Computational screening and literature appraisal complete; no causal mechanism claimed; third cohort temporal/TLike and second cohort source annotation unavailable as documented'}
(O/'completion.json').write_text(json.dumps(completion,ensure_ascii=False,indent=2),encoding='utf-8')
(D/'阅读说明.txt').write_text('已完成v16重筛。请将PDF与完整统计附件一起交给独立评审。\n55个三项交集基因；优先IL1R2、DDIT4、SKIL、FOXP1，GSN备选。\n创新性为问题分层，不将文献零命中视作首次发现。\n详细方法、完整P/FDR、历史去留和阴性/矛盾结果均在附件。\n',encoding='utf-8')
print(json.dumps(completion,ensure_ascii=False,indent=2))