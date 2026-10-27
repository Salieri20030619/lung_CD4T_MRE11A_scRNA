from pathlib import Path
import json,hashlib,zipfile,subprocess,time,shutil
O=Path(__file__).parent;H=O.parent;m=json.loads((O/'report_delivery.json').read_text(encoding='utf-8'));p=Path(m['pdf'])
q=json.loads((O/'final_QA.json').read_text(encoding='utf-8'));q['visual_review']='all 46 pages reviewed; revised figures and final changed pages rechecked; no clipping or overlap';(O/'final_QA.json').write_text(json.dumps(q,ensure_ascii=False,indent=2),encoding='utf-8')
res='v27修订及复核记录\n\nV26-R01 已解决：摘要、30页连接表、33页研究路径统一主线；细胞表型、线粒体与焦亡分开判定。\nV26-R02 已解决：31页SKIL补充FDR=0.8681；完整50项BH复算一致，无不同分析类型混入。\nV26-R03 呈现问题已解决，直接作用未获证实：27页H2-DM单列抗原呈递信息；追溯CellChatDB小鼠原条目、KEGG注释及内体/溶酶体定位；不作为普通表面信号推荐，不宣称直接结合已验证。原网络统计未重算。\nV26-R04 已解决：27页TTC39C标NA，区别于显示范围没有支持路线。\nV26-R05 已解决：正文与筛选图清理内部措辞，背景指标按预定义标志UMI比例说明。\nV26-R06 已采纳：主文保留交集示意，36—38页保留详细阈值和完整集合，45—46页统一索引；TF与轨迹保留以满足已要求分析范围。\nE01 已解决：39—44页实验附录为建议研究问题与可选路径；先模型/可测性，再选一个候选功能干预。\nE02 已解决：41页IL1R2不以常规抑制试验为唯一判定，加入IL-1环境、状态与数量，并区分洗涤后功能和持续共培养效应。\nE03 已解决：42页TTC39C优先蛋白/递送/基础存活，作用方向按实测；不指定p53。\nE04 已解决：43页SKIL成熟Treg与诱导分化、双基因模型分开。\nE05 已解决：模型阶段有限线粒体/焦亡可测性核查；执行轴依据目标T细胞表达，混合上清不能归属T细胞。\nE06 已解决：40页报告各组、主效应及交互，说明结局模型和Δ尺度；交互阴性不等同基因无功能。\nE07 已解决：平台、细胞量、递送、周期、预算与正式样本量在执行范围确定后落实；扩展机制、动物和患者验证不作全套承诺。\n\n验收范围：报告文字/数字对应、变更图示、完整50项FDR复算、101项补充材料文件完整性、46页渲染与27图计数。未重跑原始测序分析或全部细胞网络；未执行湿实验。\n'
(O/'修订复核记录.txt').write_text(res,encoding='utf-8')
zip_path=O/'v27_补充统计与图值.zip'
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for f in (O/'supplement').rglob('*'):
  if f.is_file():z.write(f,str(f.relative_to(O/'supplement')))
with zipfile.ZipFile(zip_path) as z:assert z.testzip() is None
remote='/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/mainline_revision_20261017_v27'
def run(args):
 for n in range(4):
  r=subprocess.run(args,capture_output=True,text=True,encoding='utf-8',creationflags=subprocess.CREATE_NO_WINDOW)
  if r.returncode==0:return r.stdout
  if n==3:raise RuntimeError(r.stderr)
  time.sleep(2)
ssh=['ssh','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10'];scp=['scp','-o','BatchMode=yes','-o','ConnectTimeout=20']
run(ssh+['mkdir -p '+remote+'/delivery '+remote+'/code '+remote+'/evidence'])
run(scp+[str(p),'sgs10:'+remote+'/delivery/report_v27.pdf'])
run(scp+[str(zip_path),'sgs10:'+remote+'/supplement_v27.zip'])
run(scp+[str(f) for f in O.glob('*.py')]+['sgs10:'+remote+'/code/'])
run(scp+[str(f) for f in O.glob('*.json')]+[str(O/'修订复核记录.txt'),'sgs10:'+remote+'/'])
run(scp+[str(f) for f in (O/'evidence').glob('*') if f.is_file()]+['sgs10:'+remote+'/evidence/'])
remote_hash=run(ssh+['sha256sum '+remote+'/delivery/report_v27.pdf '+remote+'/supplement_v27.zip']);assert m['sha256'] in remote_hash;assert hashlib.sha256(zip_path.read_bytes()).hexdigest() in remote_hash
c={'status':'complete','pdf':str(p),'pages':46,'supplement_zip':str(zip_path),'supplement_zip_bytes':zip_path.stat().st_size,'remote':remote,'remote_hash_verified':True,'old_version_preserved':True,'statistics_recomputed':'BH audit only; no changed inference values','wet_lab_experiments_performed':False}
(O/'completion.json').write_text(json.dumps(c,ensure_ascii=False,indent=2),encoding='utf-8');run(scp+[str(O/'completion.json'),'sgs10:'+remote+'/']);print(json.dumps(c,ensure_ascii=False))