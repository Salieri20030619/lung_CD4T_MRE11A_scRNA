from pathlib import Path
import json,hashlib,zipfile,shutil,math
import numpy as np,pandas as pd,openpyxl,fitz
from PIL import Image,ImageOps,ImageDraw
R=Path('C:/codex_test');T=R/'tmp/silicosis_v11';O=R/'output/silicosis_v11';OLD=R/'output/silicosis_v10'
payload=json.loads((T/'workbook_data.json').read_text(encoding='utf-8'))
xlsx=O/'supplement/S17_candidate_selection_rationale.xlsx';wb=openpyxl.load_workbook(xlsx,read_only=False,data_only=True)
verified=0
for d in payload['sheets']:
 sh=wb[d['name']]
 assert [sh.cell(4,c+1).value for c in range(len(d['columns']))]==d['columns']
 for i,row in enumerate(d['rows'],5):
  for j,val in enumerate(row,1):
   got=sh.cell(i,j).value
   if isinstance(val,(int,float)) and not isinstance(val,bool):assert isinstance(got,(int,float)) and math.isclose(val,got,rel_tol=1e-14,abs_tol=0),(d['name'],i,j,val,got)
   elif val=='':assert got in ['',None]
   else:assert got==val,(d['name'],i,j,val,got)
   verified+=1
 assert sh.freeze_panes is not None
# Read source statistics independently, rather than trusting the workbook intermediate.
full=pd.read_csv(O/'supplement/S17_完整效应比较.tsv',sep='\t');base=pd.read_csv(O/'supplement/S17_候选基准.tsv',sep='\t')
for scope,d in full.groupby('状态'):
 raw=pd.read_csv(OLD/f'source_tables/03_receiver/{scope}_DE_all_genes.tsv.gz',sep='\t').set_index('canonical')
 for _,r in d.iterrows():
  q=raw.loc[r['基因']];assert np.allclose([r.log2FC,r['原始P'],r['原BH']],[q.logFC,q.PValue,q.FDR],rtol=1e-13,atol=0)
  assert q.PValue<.05
  if pd.notna(r['同方向且效应更低的候选']):
   for g in r['同方向且效应更低的候选'].split(';'):
    b=raw.loc[g];assert b.PValue<.05 and b.logFC*q.logFC>0 and abs(b.logFC)<abs(q.logFC)
# Original reports and data are preserved; no figures were regenerated in this update.
for folder in ['assets','figures','plot_data','source_tables']:
 for p in (OLD/folder).rglob('*'):
  if p.is_file():assert hashlib.sha256(p.read_bytes()).digest()==hashlib.sha256((O/folder/p.relative_to(OLD/folder)).read_bytes()).digest()
for p in (OLD/'supplement').iterdir():assert p.read_bytes()==(O/'supplement'/p.name).read_bytes()
hist=json.loads((O/'provenance/screening_history_excerpts.json').read_text(encoding='utf-8'));messages=json.loads((T/'shared_messages.json').read_text(encoding='utf-8'));byid={m['id']:m for m in messages}
for h in hist:assert h['text']==byid[h['message_id']]['text']
for name in ['prepare_s17.py','build_s17.mjs','build_reports.py','extract_share.py','finalize_v11.py']:shutil.copy2(T/name,O/'code'/name)
checks={'xlsx_cells_readback_verified':verified,'numeric_roundtrip_relative_tolerance':1e-14,'comparison_rows_verified_against_original_DE':len(full),'same_direction_comparisons_verified':True,'history_excerpts_match_shared_messages':len(hist),'report_pages':len(fitz.open(O/'silicosis_integrated_report_v11.pdf')),'S17_reading_pages':len(fitz.open(O/'supplement/S17_candidate_selection_rationale.pdf')),'original_figures_and_statistics_preserved':True,'new_inferential_tests':0,'xlsx_export_note':'Exporter emitted XLSX EXPORTED but process exited 1; saved ZIP integrity, all cell values and sheet structure independently verified.','visual_review':'S17 PDF all pages and all workbook sheets inspected; main-report changed pages inspected.'}
(O/'provenance/final_validation.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2),encoding='utf-8')
(O/'README.txt').write_text('矽肺候选选择依据补充版 v11\n\n主报告：silicosis_integrated_report_v11.pdf（27页），新增第10页3.3节。\n补充表：supplement/S17_candidate_selection_rationale.xlsx；同名PDF为12页阅读版。\nExcel包含27个高效应基因的逐项说明、历史筛选依据、54个基因-状态详细统计、44个候选-状态基准、1821个完整效应比较条目和使用说明。\n\n已恢复用户提供的“数据集质量比较分析”原始对话，提取H01-H06筛选相关记录；历史决定、本次文献评估和未解决的来源问题分开表述。正文4处筛选叙述已据此修正，变更见provenance/narrative_changes.json。未取得逐基因历史排除决定的条目不追溯编造。\n原图1-15、原始统计值及全部既有补充表保持不变。没有新增统计检验或重新指定主靶点。\n文献及数据库来源共20项，含疾病/细胞场景及链接；完整对话未打包，仅保留任务相关摘录。\n',encoding='utf-8')
rows=[]
for p in sorted(O.rglob('*')):
 if p.is_file() and p.name!='MANIFEST_SHA256.tsv' and not p.name.endswith('.inspect.ndjson'):
  rows.append({'file':p.relative_to(O).as_posix(),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
pd.DataFrame(rows).to_csv(O/'MANIFEST_SHA256.tsv',sep='\t',index=False,encoding='utf-8')
zp=O.parent/'silicosis_candidate_rationale_v11.zip'
with zipfile.ZipFile(zp,'w',zipfile.ZIP_DEFLATED) as z:
 for row in rows:z.write(O/row['file'],row['file'])
 z.write(O/'MANIFEST_SHA256.tsv','MANIFEST_SHA256.tsv')
with zipfile.ZipFile(zp) as z:assert z.testzip() is None
print(json.dumps(checks,ensure_ascii=False));print('ZIP_SHA256',hashlib.sha256(zp.read_bytes()).hexdigest())
# Contact sheet of all changed report pages, excluding unchanged figures.
nums=[1,2,4,7,10,25,26];c=Image.new('RGB',(1160,1680),'#d8e0e5')
for j,n in enumerate(nums):
 p=T/f'report_{n:02d}.png';c.paste(ImageOps.contain(Image.open(p),(560,396)),((j%2)*580,(j//2)*420+20));ImageDraw.Draw(c).text(((j%2)*580+8,(j//2)*420+3),p.stem,fill='black')
c.save(T/'main_changed_pages.png')
