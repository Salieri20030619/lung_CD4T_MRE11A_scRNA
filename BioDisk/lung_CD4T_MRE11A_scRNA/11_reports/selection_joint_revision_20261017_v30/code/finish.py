from pathlib import Path
import json,hashlib,shutil,zipfile,re,subprocess,time
import pandas as pd,fitz
O=Path(__file__).parent;H=O.parent;B=H/'lung_citation_revision_20261017_v29';m=json.loads((O/'report_delivery.json').read_text(encoding='utf-8'));p=Path(m['pdf']);d=fitz.open(p);assert len(d)==56
q=json.loads((O/'layout_QA.json').read_text(encoding='utf-8'));assert not q['bounds'] and q['visuals']==27
oldp=Path(json.loads((B/'report_delivery.json').read_text(encoding='utf-8'))['pdf']);oldhash=hashlib.sha256(oldp.read_bytes()).hexdigest();assert oldhash==json.loads((B/'report_delivery.json').read_text(encoding='utf-8'))['sha256']
t='\n'.join(pg.get_text() for pg in d)
for x in ['待排','历史基因','继续硬筛','核心必做']:assert x not in t,x
for x in ['0.0375','0.675','0.8681','候选比较对象11','高效应核查信号27','定向补核对象4']:assert x in t,x
reply='''v30 筛选过程与联合样本分析修订说明
报告56页；以下页码为v30。

一、42→90→9→3的实际过程
第3、43页说明42个来源：11个候选比较对象 + 27个高效应/细胞来源核查信号 + 4个定向补核对象（DUSP5、GSDMD、CXCR6、MRE11）。已与来源表逐项匹配，无组间重叠。
与55个三项交集重叠7个，共90个。资格分层为52/37/1。
9个是追加关联分析前确定的定向面板：5个功能机制候选 + 4个补充探索对象，不是90个按统一P阈值筛剩9个。第33—34页列九者的提名问题、后续结果、保留/推荐理由；第45页流程图把资格分层与提名分开。
最终三个按结果与功能可测性对应不同问题：IL1R2的Treg/IL-1环境，TTC39C的存活/损伤探索，SKIL的状态调控。保留MDFIC等较强关联及其他六个的主要限制，未宣称三者所有指标均最强。各分析段加入其对候选判断的具体贡献。

二、3+3与5+5
这里3+3双侧组标签精确检验有20种分配，最小P=0.10；分层5+5有80种分配，最小P=0.025。此限制不适用于所有检验，edgeR差异和秩关联另有模型/排列空间。
原联合程序评分已存在但未充分展示。本轮复算主队列和联合模型共104项程序效应、P及各自完整52项BH家族，与原统计一致；第25—26页并列两模型。联合焦亡及线粒体应激没有普遍疾病增强支持；Treg凋亡P=0.050、FDR=0.144，未达到既定P<0.05。
本轮新增3项联合CD4比例补充分析：按队列/时间分层、保持层内组别数置换80种，Naive/ActCD4/Treg三项独立构成BH家族。Naive下降38.52个百分点、Treg上升25.61个百分点，均P=0.025/FDR=0.0375；ActCD4上升12.91个百分点，P=0.150/FDR=0.150。图4展示全部10个样本并保留3+3对照。

三、复核范围
新增统计仅为上述联合细胞比例补充；既有DE、候选关联、程序评分原值不改。阈值、候选集合与推荐顺序不变，文献[3]来源注记保留，v29原版保留。
已检查56页渲染、27幅图示、连续图号和索引。请重点复查第3、12—13、25—26、33—36、43—45页及新比例分析的分层模型适用性。
'''
rpath=p.parent/'v30_筛选与联合分析_简要复核说明_20261017.txt';rpath.write_text(reply,encoding='utf-8')
S=O/'supplement';shutil.copytree(B/'supplement',S,dirs_exist_ok=True)
for f in (O/'plot_data').glob('*'):shutil.copy2(f,S/'plot_data'/f.name)
(S/'tables/v30').mkdir(exist_ok=True)
for f in (O/'tables').glob('*'):shutil.copy2(f,S/'tables/v30'/f.name)
for f in (O/'evidence').glob('*'):shutil.copy2(f,S/'evidence'/f.name)
shutil.copy2(H/'lung_append_analysis_20261017_v18/protocol_frozen.json',S/'evidence/nine_panel_protocol_before_associations.json')
shutil.copy2(H/'lung_three_cohort_20261017_v16/finalize_tables.py',S/'evidence/source_nomination_rules.py')
shutil.copy2(H/'lung_dual_benchmark_20261017_v17/report.py',S/'evidence/source_five_plus_four_context.py')
shutil.copy2(H/'lung_three_cohort_20261017_v16/final_tables/优先验证候选_5基因.tsv',S/'evidence/five_nominees_source.tsv')
# A self-contained replication of the newly added proportions, using archived10-sample inputs.
code="""from pathlib import Path
import pandas as pd,numpy as np,itertools
from scipy.stats import false_discovery_control
P=Path(__file__).resolve().parent
w=pd.read_csv(P/'plot_data/F04.tsv',sep='\\t');out=[]
for co in ['314','Joint']:
 z=(w if co=='Joint' else w[w.dataset=='GSE314943']).reset_index(drop=True)
 labels=(z.dataset.astype(str)+'_'+z.day.astype(str)) if co=='Joint' else pd.Series(['all']*len(z))
 C=np.column_stack([np.ones(len(z)),pd.get_dummies(labels,drop_first=True,dtype=float).to_numpy()]);M=np.eye(len(z))-C@np.linalg.pinv(C);g=z.group.eq('Silica').to_numpy(float);gx=M@g
 groups=[np.flatnonzero(labels.eq(k)) for k in sorted(labels.unique())];choices=[list(itertools.combinations(a,int(g[a].sum()))) for a in groups]
 for sc in ['Naive','ActCD4','Treg']:
  y=(z[sc]/z.Core).to_numpy();ry=M@y;effect=float(gx@ry/(gx@gx));null=[]
  for parts in itertools.product(*choices):
   a=np.zeros(len(z))
   for ids in parts:a[list(ids)]=1
   x=M@a;null.append(float(x@ry/(x@x)))
  out.append(dict(cohort=co,scope=sc,effect=effect,P=np.mean(np.abs(null)>=abs(effect)-1e-12),allocations=len(null)))
a=pd.DataFrame(out)
for k,z in a.groupby('cohort'):a.loc[z.index,'FDR']=false_discovery_control(z.P)
ref=pd.read_csv(P/'tables/v30/proportion_joint_supplement.tsv',sep='\\t');joined=a.merge(ref,on=['cohort','scope'],suffixes=('_recomputed','_saved'),validate='one_to_one')
for col in ['effect','P','FDR','allocations']:assert np.allclose(joined[col+'_recomputed'],joined[col+'_saved'])
print(a.to_string(index=False));print('MATCHED_SAVED_RESULTS')
"""
(S/'reproduce_proportions.py').write_text(code,encoding='utf-8')
r=subprocess.run([str(H/'anaconda3/python.exe'),str(S/'reproduce_proportions.py')],capture_output=True,text=True,encoding='utf-8',creationflags=subprocess.CREATE_NO_WINDOW);assert r.returncode==0,r.stderr;(O/'evidence/proportion_replication.log').write_text(r.stdout,encoding='utf-8')
manifest=json.loads((O/'figure_manifest.json').read_text(encoding='utf-8'));pages={}
for i,pg in enumerate(d):
 for f in re.findall(r'(?m)^图(\d+)。',pg.get_text()):pages[int(f)]=i+1
mapping=[{'figure_number':x['number'],'pdf_page':pages[x['number']],'data_file':'plot_data/'+x['id']+'.tsv','title':x['title']} for x in sorted(manifest,key=lambda x:x['number'])];pd.DataFrame(mapping).to_csv(S/'figure_mapping.tsv',sep='\t',index=False)
(S/'README.txt').write_text('v30配套复核材料\n\nfigure_mapping.tsv映射连续图号、PDF页码及稳定数据ID。\n新增：tables/v30/proportion_joint_supplement.tsv为3项联合细胞比例检验及3+3对照；plot_data/F04.tsv包含全部10样本原始细胞数与比例。运行reproduce_proportions.py可独立复算。\n原有联合程序结果：group_scores_both_models.tsv；104项复算见evidence/group_score_recalculation.tsv。图15—16改为3+3与分层5+5并列。\n选择来源：supplement42_provenance.tsv逐项记录11+27+4；nine_gene_comparison_statistics.tsv为九候选所述关联原值。source_nomination_rules.py、five_nominees_source.tsv及nine_panel_protocol_before_associations.json保留定向提名和追加分析先后依据。\n候选来源证据文件保留原始字段与表述，仅供复核；正式报告使用当前定义与证据定位。\n所有原统计表与阴性结果继续保留；本包不含全部原始测序数据，不等同FASTQ全流程重跑。\n',encoding='utf-8')
inv=[{'path':str(f.relative_to(S)),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in S.rglob('*') if f.is_file() and f.name!='file_inventory.json'];(S/'file_inventory.json').write_text(json.dumps(inv,ensure_ascii=False,indent=2),encoding='utf-8')
zpath=O/'v30_补充统计与筛选来源.zip'
with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for f in S.rglob('*'):
  if f.is_file():z.write(f,str(f.relative_to(S)))
with zipfile.ZipFile(zpath) as z:assert z.testzip() is None
qa={'pages':56,'visuals':27,'layout_passed':True,'visual_review':'all56 pages reviewed; no clipping/overlap after layout repairs','old_v29_preserved':True,'old_v29_sha256':oldhash,'source42_exact_match':True,'existing_group_score_tests_reproduced':104,'new_joint_proportion_tests':3,'new_proportion_replication_passed':True,'nominal_threshold_unchanged':True,'candidate_set_and_order_unchanged':True,'pdf_sha256':hashlib.sha256(p.read_bytes()).hexdigest()};(O/'final_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
remote='/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/selection_joint_revision_20261017_v30'
def run(args):
 for n in range(4):
  r=subprocess.run(args,capture_output=True,text=True,encoding='utf-8',creationflags=subprocess.CREATE_NO_WINDOW)
  if r.returncode==0:return r.stdout
  if n==3:raise RuntimeError(r.stderr)
  time.sleep(2)
ssh=['ssh','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10'];scp=['scp','-o','BatchMode=yes','-o','ConnectTimeout=20'];run(ssh+['mkdir -p '+remote+'/delivery '+remote+'/code'])
for f,tail in [(p,'delivery/report_v30.pdf'),(rpath,'delivery/review_response.txt'),(zpath,'supplement_v30.zip')]:run(scp+[str(f),'sgs10:'+remote+'/'+tail])
run(scp+[str(f) for f in O.glob('*.py')]+['sgs10:'+remote+'/code/']);run(scp+[str(f) for f in O.glob('*.json')]+['sgs10:'+remote+'/'])
checks=run(ssh+['sha256sum '+remote+'/delivery/report_v30.pdf '+remote+'/delivery/review_response.txt '+remote+'/supplement_v30.zip'])
for f in [p,rpath,zpath]:assert hashlib.sha256(f.read_bytes()).hexdigest() in checks
c={'status':'complete','pdf':str(p),'review_txt':str(rpath),'supplement':str(zpath),'remote':remote,'sha256_verified':True,'pages':56};(O/'completion.json').write_text(json.dumps(c,ensure_ascii=False,indent=2),encoding='utf-8');run(scp+[str(O/'completion.json'),'sgs10:'+remote+'/']);print(json.dumps(c,ensure_ascii=False))