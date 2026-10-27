from pathlib import Path
import json,hashlib,zipfile,subprocess,datetime
O=Path(r'C:\Users\31937\lung_dual_benchmark_20261017_v17');B=Path(r'C:\Users\31937\lung_three_cohort_20261017_v16');j=json.loads((O/'delivery.json').read_text(encoding='utf-8'));D=Path(j['directory']);pdf=Path(j['pdf']);files={}
for folder in ['assessment','tables']:
 for f in (O/folder).rglob('*'):
  if f.is_file():files[folder+'/'+f.relative_to(O/folder).as_posix()]=f
for n in ['protocol_v17.json','genes.json','QA.json']:files[n]=O/n
for n in ['expression.py','network.R','assess.py','report.py','verify.py','package_sync.py']:files['code/'+n]=O/n
for n in ['DE_all_models.tsv','nominal_DE_axis_membership.tsv','function_membership_extended.tsv','expressed_LR_with_decoy_literature.tsv','NicheNet_activity_all.tsv','QA_statistical_integrity.tsv']:
 files['baseline_v16/'+n]=B/'data/tables'/n
for n in ['优先验证候选_5基因.tsv','四轴候选总表_55基因.tsv','历史候选逐项回溯.tsv']:files['baseline_v16/'+n]=B/'final_tables'/n
readme='''v17双基准二次评估 | 2026-10-17
A=候选入池；B=机制主候选。90基因=旧55交集与42历史基因的并集。
A1 52个：计算支持的候选，不代表已经证明巨噬导致T细胞损伤。
A2 37个：条件保留待补，不是删除。A3 1个：本轮未达名义差异或不可评估。
B 0个：目前缺目标T细胞自身损伤/死亡终点及连续连接，不是FDR原因，也不是无生物学作用。

先读PDF，再读assessment/90基因双基准总表.tsv与逐状态双基准评估.tsv。
功能先验与T细胞自身终点分开；给出的作用方向/来源缺口逐项保留。
NicheNet目标排名及配体排名不作硬排除；非零权重只表示先验覆盖，不能认为显著、因果或调控方向已知。
全部P/FDR复用v16，完整检验家族在baseline_v16，未新增合并P值。
参考表内best_target_rank和best_activity_rank是各自最小值，不能拼成同一条边；确切配对见network_links_no_rank_exclusion.tsv.gz。
v14完整模块/近端/复合物/去重叠表在assessment。轨迹与人类旧结果不得外推为新基因已验证。
源代码路径冻结到本项目。服务器v16保留原始矩阵和先验；v17为独立复评目录。
新增论文与角色复核的具体来源见assessment/新增文献及细胞作用角色复核.tsv及PDF。
全部文本为UTF-8无BOM。完整源文件逐项哈希见MANIFEST_SHA256.json。
'''
(O/'README.txt').write_text(readme,encoding='utf-8');files['README.txt']=O/'README.txt'
manifest={n:{'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'bytes':f.stat().st_size} for n,f in files.items()};zp=D/'v17_双基准复评完整证据.zip'
with zipfile.ZipFile(zp,'w',zipfile.ZIP_DEFLATED) as z:
 for n,f in files.items():z.write(f,n)
 z.writestr('MANIFEST_SHA256.json',json.dumps(manifest,ensure_ascii=False,indent=2).encode('utf-8'))
with zipfile.ZipFile(zp) as z:
 assert z.testzip() is None
 for n,m in manifest.items():assert hashlib.sha256(z.read(n)).hexdigest()==m['sha256']
complete={'version':'v17','finished_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'assessment_complete','genes':90,'scope_contrasts':344,'A1':52,'A2':37,'A3':1,'B_qualified':0,'pdf_pages':9,'visual_QA':'all pages rendered and inspected','outputs':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in [pdf,zp]}}
(O/'completion.json').write_text(json.dumps(complete,ensure_ascii=False,indent=2),encoding='utf-8')
root='/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/dual_benchmark_20261017_v17';flags={'creationflags':subprocess.CREATE_NO_WINDOW,'capture_output':True};subprocess.run(['ssh','-T','sgs10','mkdir -p '+root+'/delivery'],check=True,**flags)
for f,n in [(pdf,'delivery/v17_report_zh.pdf'),(zp,'delivery/v17_evidence.zip'),(O/'completion.json','completion_local.json')]:subprocess.run(['scp','-q',str(f),'sgs10:'+root+'/'+n],check=True,**flags)
remote="""from pathlib import Path
import json,zipfile,hashlib
O=Path('ROOT');j=json.loads((O/'completion_local.json').read_text());checks={}
for name,sha in j['outputs'].items():
 p=O/'delivery'/('v17_report_zh.pdf' if name.endswith('.pdf') else 'v17_evidence.zip');actual=hashlib.sha256(p.read_bytes()).hexdigest();assert actual==sha;checks[p.name]=actual
with zipfile.ZipFile(O/'delivery/v17_evidence.zip') as z:
 for n in z.namelist():
  if n.startswith(('assessment/','code/')) or n in ['QA.json','README.txt']:
   dest=(O/n).resolve();assert dest.is_relative_to(O.resolve());dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(z.read(n))
j['server_hash_verification']=checks;(O/'FINAL_COMPLETION.json').write_text(json.dumps(j,ensure_ascii=False,indent=2),encoding='utf-8');print('SERVER_VERIFIED',checks)
""".replace('ROOT',root)
r=subprocess.run(['ssh','-T','sgs10','/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python -'],input=remote.encode('utf-8'),check=True,**flags);print(r.stdout.decode('utf-8'));print(json.dumps(complete,ensure_ascii=False,indent=2))