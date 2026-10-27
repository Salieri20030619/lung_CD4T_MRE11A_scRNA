"""Delivery-only checks, readable numerical brief and archive; never refit."""
from pathlib import Path
import json,hashlib,zipfile,datetime
from html.parser import HTMLParser
from collections import Counter
import pandas as pd
B=Path(__file__).resolve().parents[1];D=B/'deliverables';R=B/'results';P=B/'provenance'
def rd(n):return pd.read_csv(R/n,sep='\t',low_memory=False)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
class Inspect(HTMLParser):
 def __init__(self):super().__init__();self.tags=Counter();self.external=[]
 def handle_starttag(self,tag,attrs):
  self.tags[tag]+=1
  for k,v in attrs:
   if k in ['src','href'] and v and v.startswith(('http://','https://')):self.external.append(v)
h=D/'TNF_IL18_donor_bootstrap_zh.html';txt=h.read_text(encoding='utf-8');a=Inspect();a.feed(txt)
assert a.tags['svg']==4 and a.tags['table']>=10 and not a.external
assert '\ufffd' not in txt and '旧单次 P/FDR 完整保留' in txt
frozen=json.loads((P/'FROZEN.json').read_text());assert sha(B/'FROZEN_PLAN.txt')==frozen['plan_sha256'];assert sha(B/'code/bootstrap.py')==frozen['code_sha256']
for i in frozen['inputs']:assert sha(B/i['path'])==i['sha256']
s=rd('effect_summary.tsv');x=rd('paired_difference_summary.tsv');c=rd('chain_summary.tsv');alls=rd('all_seeds_chain_summary.tsv')
E=['ligand_TF','TF_phenotype','ligand_phenotype'];EN={'ligand_TF':'配体—TF','TF_phenotype':'TF—效应','ligand_phenotype':'配体—效应'}
lines=['TNF/IL18 供体配对 bootstrap：最终数值摘要','',
'分层单元顺序：G135-Control、G227-Control、G227-CWP、G288-CWP。RNA19=2/8/2/7；17人=1/7/2/7。17人单例CTL068固定。主区间条件于当前来源/疾病构成，不是独立来源验证。',
'每方案/名单5000次；主分层+预定不分层敏感性。RNA17/AUC三seed、两条线、三边和调整共享索引；RNA19另抽。seed17沿用旧主定义，29/43并列，不按结果选择。',
'分层主分析：75个展示规格；原效应/中位数/区间/符号保持见完整effect_summary.tsv。以下列全部调整规格。',
'边顺序：配体—TF；TF—效应；配体—效应。每个区间为95%百分位区间；有效率均以5000次为分母。','']
brief=[]
for axis,adj in [('IL18','B'),('TNF','B'),('TNF','B_Treg')]:
 for m in ['RNA19','RNA17','AUC17','AUC29','AUC43']:
  q=s[(s.scheme=='stratified')&(s.axis==axis)&(s.adjustment==adj)&(s.modality==m)].set_index('edge');cc=c[(c.scheme=='stratified')&(c.axis==axis)&(c.adjustment==adj)&(c.modality==m)].iloc[0]
  out={'axis':axis,'modality':m,'adjustment':adj}
  values=[]
  for e in E:
   z=q.loc[e];val=f'{z.original_effect:.3f} [{z.CI95_low:.3f}, {z.CI95_high:.3f}]';values.append(val);out[e+'_rho_CI95']=val
  out['all_edges_valid_rate']=cc.valid/5000;out['joint_compatible_fraction']=cc.direction_compatible_fraction;out['all_original_sign_fraction']=cc.all_three_original_sign_fraction;brief.append(out)
  lines.append(f'{axis} {m} {adj}：'+'；'.join(values)+f'；三边均有效 {int(cc.valid)}/5000 ({cc.valid/5000:.2%})；同次方向相容 {cc.direction_compatible_fraction:.2%}；三边均保原符号 {cc.all_three_original_sign_fraction:.2%}。')
pd.DataFrame(brief).to_csv(D/'adjusted_summary_matrix.tsv',sep='\t',index=False)
lines+=['','15项Treg增强：Δρ=ρ(B+Treg)−ρ(B)，同次共同有效配对，全部95%区间跨零。']
q=x[(x.scheme=='stratified')&(x.comparison_type=='Treg_enhancement')];assert len(q)==15 and ((q.CI95_low<0)&(q.CI95_high>0)).all()
for m in ['RNA19','RNA17','AUC17','AUC29','AUC43']:
 z=q[q.modality==m].set_index('edge');lines.append(m+'：'+'；'.join(f'{EN[e]} Δ={z.loc[e].original_effect:.3f} [{z.loc[e].CI95_low:.3f}, {z.loc[e].CI95_high:.3f}]，有效{int(z.loc[e].valid)}/5000，增加比例{z.loc[e].positive_fraction:.2%}' for e in E)+'。')
lines+=['','IL18基础调整seed配对差（同次索引，后者减前者；末边恒等为0）']
for comp in ['AUC29-AUC17','AUC43-AUC17','AUC43-AUC29']:
 z=x[(x.scheme=='stratified')&(x.axis=='IL18')&(x.adjustment=='B')&(x.comparison==comp)].set_index('edge')
 lines.append(comp+'：'+'；'.join(f'{EN[e]} Δ={z.loc[e].original_effect:.3f} [{z.loc[e].CI95_low:.3f}, {z.loc[e].CI95_high:.3f}]，有效{int(z.loc[e].valid)}/5000' for e in E)+'。')
lines+=['','三个seed同次均方向相容：']
for z in alls[(alls.scheme=='stratified')&(alls.adjustment!='unadjusted')].itertuples():lines.append(f'{z.axis} {z.adjustment}: {z.all_seeds_compatible_count}/{z.valid} = {z.all_seeds_compatible_fraction:.2%}；三个seed各自三边都保原符号 {z.all_seeds_original_sign_fraction:.2%}。')
lines+=['','失败与解释边界：',
'分层无秩亏；17人TNF B+Treg有2/5000次零残差失败，其他主规格有效5000/5000。分层独特供体中位数RNA19=13、17人=12；重复权重有效数中位数10.94/9.97。17人B+Treg有16.44%的尝试支持诊断df<3，未按此筛选区间。',
'不分层：RNA19 B有效3793/5000(75.86%)，TNF B+Treg有效3792/5000(75.84%)；17人B有效2807/5000(56.14%)，B+Treg有效2805/5000(56.10%)。高失败率意味着有效重复是条件子集，不能只据其区间作无条件稳健结论；全部失败原因和全尝试方向上下界均保留。',
'方向比例、增加比例不是P值或真实/因果概率；三边代数相容不是中介证明。bootstrap不增加样本量，不单独解释regulon生成随机性。小n、离散秩、单人层和近退化支持度限制区间覆盖。',
'旧单次P/FDR完整保存在inputs/aligned_chain_models.tsv；本次对应旧统计另存results/legacy_statistics_REUSED.tsv，保持旧family和计算方法。bootstrap不取代旧P/FDR，也不新算P显著比例。',
'统计核验PASS。4张PNG已目视检查；HTML结构/自包含检查通过。HTML整页浏览器截图受本机Edge进程/GPU限制未完成，未宣称完整页面像素级核查。',
'完整报告：TNF_IL18_donor_bootstrap_zh.html；完整代码、输入、索引和逐次结果：TNF_IL18_bootstrap_reproducible.zip。']
(D/'SUMMARY_zh.txt').write_text('\n'.join(lines),encoding='utf-8')
qa=dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),statistics='PASS_PREVIOUSLY_COMPLETED_NO_RERUN',frozen_inputs_plan_estimation_code_unchanged=True,html_utf8_valid=True,inline_svg_count=a.tags['svg'],html_table_count=a.tags['table'],external_resources=len(a.external),four_exported_PNGs_visually_inspected=True,full_page_browser_screenshot='NOT_COMPLETED: local Edge GPU/process restrictions; figures and HTML structure verified',legacy_P_FDR_preserved=True)
(P/'DELIVERY_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
archive=D/'TNF_IL18_bootstrap_reproducible.zip'
files=[p for p in sorted(B.rglob('*')) if p.is_file() and p!=archive and p.name!='manifest.tsv' and '__pycache__' not in p.parts and p.name!='LIBRARY_DELIVERY_RECEIPT.json']
manifest=[dict(path=p.relative_to(B).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in files]
pd.DataFrame(manifest).to_csv(B/'manifest.tsv',sep='\t',index=False)
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in files+[B/'manifest.tsv']:z.write(p,p.relative_to(B).as_posix())
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None
 for r in manifest:assert hashlib.sha256(z.read(r['path'])).hexdigest()==r['sha256']
print(json.dumps(dict(qa=qa,archive_bytes=archive.stat().st_size,archive_sha256=sha(archive),archived_files=len(files)+1,report_sha256=sha(h)),indent=2))
