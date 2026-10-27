from pathlib import Path
import json,hashlib,zipfile,subprocess,datetime,sys,importlib.metadata as md
O=Path(r'C:\Users\31937\lung_append_analysis_20261017_v18');j=json.loads((O/'delivery.json').read_text(encoding='utf-8'));D=Path(j['directory']);pdf=Path(j['pdf']);base=Path(r'C:\Users\31937\Downloads\线粒体焦亡_双基准复评_v17_20261017');files={}
readme='''v18追加分析完整复核包 | 2026-10-17

阅读顺序：
1. 合并PDF第1-9页为原v17，逐页文本及渲染像素保持一致；第10-20页为v18追加。
2. review_tables/相对v17的逐基因证据更新.tsv：新增支持、矛盾、未补齐环节。
3. review_tables/候选表型损伤_主分析与敏感性对照.tsv：全部131项联合主候选检验与调整/背景/组成/时间对照。
4. review_tables/上游连接逐段证据_不跨来源拼接.tsv：来源、状态、配体、近端TF及候选；未预设TF的字段保留NA。
5. tables/all_associations.tsv（2544项、108家族）、programme_group_effects.tsv（221项）、focal_TF_disease_effects.tsv（28项）。所有计划结果均保留；P=预设置换P（时间F效应除外），P_model_t=另列近似模型P，FDR只校正对应P列。
6. common_gene_sensitivity.tsv（500项）及成员记录：覆盖差异诊断后的共同基因敏感性，不能用较小家族FDR重标主分析。
7. baseline_v17：原报告和原完整附件原样保留；旧DE全家族与90基因基准在旧附件中。

统计单位：独立小鼠。两个疾病队列联合n10按队列/时间分层；巨噬和AUCell关联为原6只小鼠；第三队列n12 TLike单独分析，n11是>=10细胞敏感性。GSE183682每时间点n1+1只给方向，无独立P。
候选与全部评分成员已从新TF签名移除；9个候选从所有模块移除。程序为RNA代理，不是死亡执行、蛋白活化或代谢通量。BG1只过滤接收端细胞背景负荷，不是去污染；发送端不因BG1改变。原六样本及独立队列/亚群不是新的多个独立验证。

特别需要核查：
- 联合条件调整只允许块内交换，四个新样本为单例单元；全部置换规则和样本名单见protocol及统计表。
- FOXP1 Core-PyroExec固定组成后减弱。
- IL1R2-ENO1是RNA正相关，不等于抑制酶活；IL1R2-PyroExec当前不支持关联。
- SKIL更支持Treg活化；p53/损伤段仍未闭合。
- TTC39C/MDFIC新增较强ApoExec关联，属于探索性候选线索；TTC39C D7方向对聚合/归一化敏感、D56反向；MDFIC第二队列两个时间点反向。
- AllMac/AMlike和全肺Broad=Mac是不同范围；不能拼接显著段或把来源表达当分泌贡献。
- 时间TLike不是纯ActCD4/Treg，不能用时间趋势充当疾病对照验证。
- 阈值、敏感性和创新性不串成所有候选必须全部显著的硬链；B2湿实验实证不作为筛选阶段交付前提。

复算：包内样本分数、发送端表达、TF均值足以复算样本级P/FDR。全细胞级原始输入、score.loom和细胞AUCell对象在服务器；其路径和SHA在inputs。若需核查原始细胞计数，须使用这些服务器输入，不能仅凭本包完成原始矩阵审计。
代码顺序：resource.R -> export_scores.py -> aucell_clean.py（原pySCENIC环境） -> analyze.py -> source_checks.py -> 下载tables -> tf_effect_and_time.py -> common_gene_sensitivity.py -> verify_numbers.py -> review_tables.py -> normalization_audit.py -> report.py -> verify_pdf.py。
本补充不修改v16/v17统计结果；所有文本UTF-8无BOM。补充是内部复审版，尚未获得独立AI评审通过，也未向外部接收方发送。
'''
(O/'README_for_review.txt').write_text(readme,encoding='utf-8')
for folder in ['tables','review_tables','inputs','logs','figures']:
 for f in (O/folder).rglob('*'):
  if f.is_file():files[folder+'/'+f.relative_to(O/folder).as_posix()]=f
for n in ['protocol_frozen.json','coverage_sensitivity_record.json','QA.json','pdf_QA.json','README_for_review.txt']:files[n]=O/n
for n in ['resource.R','export_scores.py','aucell_clean.py','analyze.py','source_checks.py','tf_effect_and_time.py','common_gene_sensitivity.py','verify_numbers.py','review_tables.py','normalization_audit.py','report.py','verify_pdf.py','final_provenance.py','package_sync.py']:files['code/'+n]=O/n
files['code/prepare_v18.py']=Path(r'C:\Users\31937\prepare_v18.py');files['baseline_v17/v17_report_zh.pdf']=base/'v17_双基准二次评估报告.pdf';files['baseline_v17/v17_evidence.zip']=base/'v17_双基准复评完整证据.zip'
meta={'python':sys.version,'numpy':md.version('numpy'),'pandas':md.version('pandas'),'scipy':md.version('scipy'),'reportlab':md.version('reportlab'),'pypdf':md.version('pypdf')};(O/'local_environment.json').write_text(json.dumps(meta,indent=2),encoding='utf-8');files['logs/local_environment.json']=O/'local_environment.json'
manifest={n:{'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'bytes':f.stat().st_size} for n,f in files.items()};zp=D/'v18追加分析_完整复核附件.zip'
with zipfile.ZipFile(zp,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for n,f in files.items():z.write(f,n)
 z.writestr('MANIFEST_SHA256.json',json.dumps(manifest,ensure_ascii=False,indent=2).encode('utf-8'))
with zipfile.ZipFile(zp) as z:
 assert z.testzip() is None
 for n,m in manifest.items():assert hashlib.sha256(z.read(n)).hexdigest()==m['sha256']
complete={'version':'v18_append_after_v17','completed_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'complete','focus_genes':9,'gene_state_pairs':10,'routes':18,'associations':2544,'programme_effects':221,'focal_TF_effects':28,'common_gene_sensitivity':500,'BH_families_verified':176,'key_permutation_P_independently_verified':4,'PDF_total_pages':20,'baseline_pages_unchanged':9,'added_pages':11,'visual_QA':'all added pages inspected; changed result pages rechecked; original9pages text and pixels identical','output_files':{f.name:{'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'bytes':f.stat().st_size} for f in [pdf,zp]}}
(O/'completion.json').write_text(json.dumps(complete,ensure_ascii=False,indent=2),encoding='utf-8');(D/'交给独立AI的阅读说明.txt').write_text(readme,encoding='utf-8')
root='/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/append_analysis_20261017_v18';flags={'creationflags':subprocess.CREATE_NO_WINDOW,'capture_output':True};subprocess.run(['ssh','-T','sgs10','mkdir -p '+root+'/delivery'],check=True,**flags)
for f,n in [(pdf,'delivery/v17_with_v18_appendix_zh.pdf'),(zp,'delivery/v18_review_evidence.zip'),(O/'completion.json','completion_local.json')]:subprocess.run(['scp','-q',str(f),'sgs10:'+root+'/'+n],check=True,**flags)
remote="""from pathlib import Path
import json,hashlib,zipfile
O=Path('ROOT');j=json.loads((O/'completion_local.json').read_text());result={}
for n,v in j['output_files'].items():
 f=O/'delivery'/('v17_with_v18_appendix_zh.pdf' if n.endswith('.pdf') else 'v18_review_evidence.zip');sha=hashlib.sha256(f.read_bytes()).hexdigest();assert sha==v['sha256'];result[f.name]=sha
with zipfile.ZipFile(O/'delivery/v18_review_evidence.zip') as z:
 for n in z.namelist():
  if n.startswith(('tables/','review_tables/','code/')) or n in ['QA.json','pdf_QA.json','README_for_review.txt','coverage_sensitivity_record.json']:
   dest=(O/n).resolve();assert dest.is_relative_to(O.resolve());dest.parent.mkdir(parents=True,exist_ok=True);data=z.read(n)
   if dest.exists():assert dest.read_bytes()==data,('existing_file_differs',n)
   else:dest.write_bytes(data)
for f in ['export_complete.txt','AUC_complete.txt','source_complete.txt']:assert (O/'logs'/f).read_text().strip()=='success'
assert json.loads((O/'logs/association_complete.json').read_text())['status']=='success'
j['server_sha256_verified']=result;(O/'FINAL_COMPLETION.json').write_text(json.dumps(j,ensure_ascii=False,indent=2),encoding='utf-8');print('SERVER_VERIFIED',result)
""".replace('ROOT',root)
r=subprocess.run(['ssh','-T','sgs10','/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python -'],input=remote.encode('utf-8'),check=True,**flags);print(r.stdout.decode('utf-8'));print(json.dumps(complete,ensure_ascii=False,indent=2))