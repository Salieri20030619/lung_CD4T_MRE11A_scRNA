from pathlib import Path
import json,shutil,hashlib,os
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/death_priority_20261018_v3';Z=P/'00_project'
d=json.loads((O/'delivery_complete.json').read_text(encoding='utf-8'))
assert d['status']=='complete' and d['remote_sha256_verified']
assert hashlib.sha256((O/'delivery/report_v3.pdf').read_bytes()).hexdigest()==d['sha256_pdf']
assert json.loads((O/'network_loo/matched_reference_complete.json').read_text())['status']=='success'
latest=Z/'LATEST_DEATH_PRIORITY.json';old=json.loads(latest.read_text(encoding='utf-8'))
assert old['version'] in ['v2','v3']
backup=Z/'LATEST_DEATH_PRIORITY_before_v3_20261018.json'
if old['version']=='v2' and not backup.exists():shutil.copy2(latest,backup)
new={'version':'v3','status':'complete','completion':str(O/'delivery_complete.json'),'report':str(O/'delivery/report_v3.pdf'),'baseline':'v1/v2 preserved','data_root':'/mnt/BioData/bio/wangmingyuan/lung_CD4T_MRE11A_scRNA','root_available_bytes_at_delivery':d['storage']['root_available_bytes_at_delivery']}
latest.write_text(json.dumps(new,ensure_ascii=False,indent=2),encoding='utf-8')
text='''线粒体/焦亡项目：凋亡焦亡优先独立路线 v3 完成
完整报告75页、31图，全部竖版，正文移除旧TF投射；v1/v2保留。
原完整5+5框架继续有效。追加：All/BG1各10次逐鼠全网络重建、活性和关联统计、匹配完整seed17参考；公开IL2RA/HIF1A真实干预、独立人肺、BAL及关键成员复验；当前TF靶程序在真实干预中的响应。
351项新增数值/完整性检查通过，含42个程序效应独立复算、40项留一关键关联和4项匹配参考关键关联的独立复算。20个网络两阶段退出码均为0。
153项临界GSEA统一提高精度，初次值保留。结构统计和展示排除软件TF自身成员，RNA-自身TF活性关联另保留。
IL2RA外周Treg敲除支持存活/抗凋亡功能。人肺13例多来源模型支持代谢关联，但G227同来源疾病上调未复现。NFKB2统计强度依赖签名构建；All留一9/10同向、0/10名义，BG1 10/10同向、7/10名义；结构边9/10与7/10保留。不能提升为矽肺内已确立的死亡机制。
ARL6IP1未找到匹配CD4/Treg直接敲除转录组；可选虚拟敲除经模型适用性检查未执行。
项目下载迁至/mnt/BioData/bio/wangmingyuan/lung_CD4T_MRE11A_scRNA，约9.03GB，经SHA256逐文件校验；BioDisk原输入路径以软链接兼容。
本账号180项历史系统盘文件/目录约3.07GB迁至BioData，保留兼容链接。交付时根分区可用字节：'''+str(d['storage']['root_available_bytes_at_delivery'])+'''；剩余系统盘占用需管理员处理，不涉及其他用户文件。
后续计算使用TMPDIR=/mnt/BioData/bio/wangmingyuan/lung_CD4T_MRE11A_scRNA/tmp/death_priority_v3。
正式完成记录：'''+str(O/'delivery_complete.json')+'\n服务器报告：'+str(O/'delivery/report_v3.pdf')+'\n'
(Z/'PROJECT_STATUS_death_priority_v3_20261018.txt').write_text(text,encoding='utf-8')
(O/'analysis_complete.json').write_text(json.dumps({'status':'success','report_status':'complete','network_folds':20,'matched_seed17_reference':True,'public_perturbation_series':3,'numerical_checks':351,'original_versions_preserved':True},indent=2),encoding='utf-8')
print('LATEST_V3_UPDATED_AND_REPORT_HASH_VERIFIED',flush=True)
