from pathlib import Path
import json,hashlib,re,shutil,zipfile,subprocess,time
import fitz,pandas as pd
O=Path(__file__).parent;H=O.parent;B=H/'lung_mainline_revision_20261017_v27';m=json.loads((O/'report_delivery.json').read_text(encoding='utf-8'));p=Path(m['pdf']);d=fitz.open(p)
assert len(d)==50 and sum(len(x.get_images()) for x in d)==27
text='\n'.join(x.get_text() for x in d);assert '待排' not in text
bounds=[]
for i,pg in enumerate(d):
 for b in pg.get_text('blocks'):
  if b[0]<39 or b[2]>pg.rect.width-37 or b[1]<24 or b[3]>pg.rect.height-15:bounds.append([i+1,b[:4]])
assert not bounds,bounds
old=Path(json.loads((B/'report_delivery.json').read_text(encoding='utf-8'))['pdf']);assert hashlib.sha256(old.read_bytes()).hexdigest()=='2f984747d3db9df763dca3347db41ef241cc54eb51cba267b6e6a9b37ec05655'
figs=[]
for i,pg in enumerate(d):
 for f in re.findall(r'(?m)^图(\d+)。',pg.get_text()):figs.append((int(f),i+1))
assert [f for f,n in figs]==list(range(1,25))
for n in ['F19','F24']:pd.testing.assert_frame_equal(pd.read_csv(O/'plot_data'/(n+'.tsv'),sep='\t'),pd.read_csv(B/'supplement/plot_data'/(n+'.tsv'),sep='\t'))
reply='''致审核者：v27可读性与逻辑意见修订回复
修订版：v28，共50页。以下均为v28页码。

1. 主线与导航：第2页前置问题—证据—阅读位置总览；第26页说明轨迹仅承担Naive/ActCD4状态背景，不验证Treg路线。
2. 候选入选解释：第3页集中列来源、分层及推荐依据。已核对成员表：IL1R2属于55交集及42补充两集合，SKIL来自55交集，二者为计算支持；TTC39C来自补充集合，功能证据待补。未将三者写成同一硬筛路径产物。
3. 线粒体结果：第1、25页补充主队列MitoISR方向及统计结果，明确其为转录应激程序，并与候选关联区分。
4A. NicheNet：第28页选取、排名、颜色统一为校正AUPR，格内标排名；第48页方法同步。沿用原排序与数值。
4B. 调整模型：第31页图22并列M1（队列/时间层）和M2（进一步疾病调整）的rho、P、FDR，醒目标示原始坐标；第33页使用同一命名。
5. 术语与状态：第7页集中解释程序及模型名称；第10页说明Core/Naive/ActCD4/Treg/TLike，并据已有标志表列描述性表达特征，不新增注释阈值。
6. 冗余与方法：详细公式、参数及备选PROS1路线说明移至第47—48页方法附录；压缩讨论中的重复候选介绍，保留正文与实验附录各自用途。
7. 完成度命名：第30—31页改称“候选表达与关联证据”；第35页明确为“实验验证建议”。
8. 版式与导航：表格字号由8.5提高至10磅；图号按阅读顺序连续编为1—24；第49—50页索引增加页码及稳定数据文件映射。已逐页检查并修复关联图说明重叠。

本轮为叙述、已有结果补充与图形呈现修订；样本、统计阈值、原P/FDR及候选推荐顺序保持不变。未新增湿实验结论，也不将此次自检表述为独立统计审计通过。原v27保留。
请重点复查第2—3、25、28、31—33、49—50页及其衔接。
'''
reply_path=p.parent/'v28_给独立审核的简要回复_20261017.txt';reply_path.write_text(reply,encoding='utf-8');assert reply_path.read_text(encoding='utf-8')==reply
# Supplement retains stable data IDs despite renumbered display figures.
S=O/'supplement';shutil.copytree(B/'supplement',S,dirs_exist_ok=True)
for f in (O/'evidence').glob('*'):shutil.copy2(f,S/'evidence'/f.name)
for f in (O/'plot_data').glob('*'):shutil.copy2(f,S/'plot_data'/f.name)
manifest=json.loads((O/'figure_manifest.json').read_text(encoding='utf-8'));mapping=[]
for x in sorted(manifest,key=lambda x:x['number']):mapping.append({'figure_number':x['number'],'pdf_page':dict(figs)[x['number']],'data_file':'plot_data/'+x['id']+'.tsv','title':x['title']})
pd.DataFrame(mapping).to_csv(S/'figure_mapping.tsv',sep='\t',index=False)
(S/'README.txt').write_text('v28配套统计与图值\n\n报告图号按阅读顺序1—24排列；数据文件保留稳定ID。请以figure_mapping.tsv核对图号、报告页码、数据路径。\n图19统一使用校正AUPR；图22（数据F24.tsv）原始点对应M1/M2统计见evidence/two_model_plot_statistics.tsv。完整P/FDR及家族表沿用原统计。\n新增候选来源和状态标志说明的依据见evidence/candidate_origin_and_status.tsv及state_marker_sample_mean.tsv；原始标志表见tables/enrichment_and_figure_sources/cd4_markers.tsv。\n其他目录保留完整统计、阴性结果、筛选成员、图值与来源。H2-DM仍单列抗原呈递网络信息；TTC39C网络NA表示未纳入评价。\n本包不包含全部原始测序文件，不等同从FASTQ重跑的全流程包。\n',encoding='utf-8')
inv=[{'path':str(f.relative_to(S)),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in S.rglob('*') if f.is_file() and f.name!='file_inventory.json'];(S/'file_inventory.json').write_text(json.dumps(inv,ensure_ascii=False,indent=2),encoding='utf-8')
zpath=O/'v28_补充统计与图值.zip'
with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for f in S.rglob('*'):
  if f.is_file():z.write(f,str(f.relative_to(S)))
with zipfile.ZipFile(zpath) as z:assert z.testzip() is None
qa={'pages':50,'visuals':27,'bounds':bounds,'figures_sequential':True,'table_font_pt':10,'source_data_F19_F24_identical':True,'v27_preserved':True,'visual_review':'50 pages reviewed; revised association page rerendered and passed','pdf_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'reply':str(reply_path),'supplement':str(zpath)};(O/'final_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
remote='/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/readability_revision_20261017_v28'
def run(args):
 for n in range(4):
  r=subprocess.run(args,capture_output=True,text=True,encoding='utf-8',creationflags=subprocess.CREATE_NO_WINDOW)
  if r.returncode==0:return r.stdout
  if n==3:raise RuntimeError(r.stderr)
  time.sleep(2)
ssh=['ssh','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10'];scp=['scp','-o','BatchMode=yes','-o','ConnectTimeout=20']
run(ssh+['mkdir -p '+remote+'/delivery '+remote+'/code'])
run(scp+[str(p),'sgs10:'+remote+'/delivery/report_v28.pdf']);run(scp+[str(reply_path),'sgs10:'+remote+'/delivery/reviewer_reply.txt']);run(scp+[str(zpath),'sgs10:'+remote+'/supplement_v28.zip'])
run(scp+[str(f) for f in O.glob('*.py')]+['sgs10:'+remote+'/code/']);run(scp+[str(f) for f in O.glob('*.json')]+['sgs10:'+remote+'/'])
checks=run(ssh+['sha256sum '+remote+'/delivery/report_v28.pdf '+remote+'/delivery/reviewer_reply.txt '+remote+'/supplement_v28.zip'])
for f in [p,reply_path,zpath]:assert hashlib.sha256(f.read_bytes()).hexdigest() in checks
complete={'status':'complete','pdf':str(p),'reply_txt':str(reply_path),'supplement':str(zpath),'remote':remote,'checksums_verified':True,'pages':50};(O/'completion.json').write_text(json.dumps(complete,ensure_ascii=False,indent=2),encoding='utf-8');run(scp+[str(O/'completion.json'),'sgs10:'+remote+'/']);print(json.dumps(complete,ensure_ascii=False))