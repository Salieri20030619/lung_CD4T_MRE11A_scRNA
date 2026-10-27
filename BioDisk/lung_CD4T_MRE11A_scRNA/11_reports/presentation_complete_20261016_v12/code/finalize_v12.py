from pathlib import Path
import json,hashlib,shutil,zipfile,re
from collections import Counter
import numpy as np,pandas as pd
from scipy.stats import spearmanr
from bs4 import BeautifulSoup
import fitz
from PIL import Image,ImageOps,ImageDraw
R=Path('C:/codex_test');O=R/'output/silicosis_v12';T=R/'tmp/silicosis_v12';V4=R/'tmp/silicosis_v9/v4';OLD=R/'output/silicosis_v11'
def rd(p):return pd.read_csv(p,sep='\t')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
checks={}
assert rd(O/'supplement/S18_sample_program_scores.tsv').equals(rd(V4/'01_phenotype/sample_program_scores.tsv'))
assert rd(O/'supplement/S18_death_program_statistics.tsv').equals(rd(V4/'01_phenotype/01_CD4_death_group_tests.tsv'))
checks['score_and_test_tables_unchanged']=True
coords=rd(O/'supplement/S18_composition_scatter_coordinates.tsv');act=rd(V4/'03_receiver/TF_all_regulon_sample_activity.tsv.gz');act=act[act.seed==17];sc=rd(V4/'01_phenotype/sample_program_scores.tsv');sc=sc[sc['filter']=='All'];w=rd(V4/'06_support/Core_fixed_composition_weights.tsv').set_index('state').fixed_weight
for _,r in coords.iterrows():
 if r.composition=='实际亚群比例':
  x=act[(act.scope=='Core')&(act['sample']==r['sample'])&(act.regulon==r.tf)].mean_auc.iloc[0];y=sc[(sc.scope=='Core')&(sc['sample']==r['sample'])][r.module].iloc[0]
 else:
  x=sum(float(act[(act.scope==s)&(act['sample']==r['sample'])&(act.regulon==r.tf)].mean_auc.iloc[0])*weight for s,weight in w.items());y=sum(float(sc[(sc.scope==s)&(sc['sample']==r['sample'])][r.module].iloc[0])*weight for s,weight in w.items())
 assert np.allclose([r.tf_auc,r.program_score],[x,y],rtol=1e-12,atol=1e-14)
checks['sample_scatter_coordinates_verified']=len(coords)
rank=rd(V4/'02_sender/source_specific_activity_ranks.tsv');shown=rd(O/'supplement/S19_NicheNet_activity_figure_values.tsv')
for scope,d in shown.groupby('scope'):
 expected=rank[(rank.sender=='AMlike')&(rank.scope==scope)&(rank.direction=='Up')].sort_values('source_specific_rank').head(8)
 assert set(d.ligand)==set(expected.ligand)
 for _,r in d.iterrows():assert np.isclose(r.aupr_corrected,expected.set_index('ligand').loc[r.ligand,'aupr_corrected'],rtol=1e-13,atol=0)
activity=rd(O/'supplement/S19_NicheNet_all_activities.tsv');assert np.allclose(activity.aupr-activity.target_n/activity.background_n,activity.aupr_corrected,rtol=1e-10,atol=1e-14)
checks['nichenet_metric_and_rank_values_verified']=len(shown)
prior=rd(T/'inputs_R2/nichenet_death_candidate_prior_weights.tsv.gz');plot=rd(O/'supplement/S19_NicheNet_death_weight_figure_values.tsv');joined=plot.merge(prior,on=['target','ligand'],suffixes=('_plot','_source'),validate='one_to_one');assert len(joined)==len(plot) and np.allclose(joined.weight_plot,joined.weight_source,rtol=1e-13,atol=0)
checks['unchanged_prior_weights_verified']=len(plot)
genes=rd(O/'supplement/S19_ActCD4_death_genes_DE_and_target_membership.tsv');assert set(genes.loc[genes.target_membership=='Up_target','gene'])=={'GSDMD'}
de=rd(V4/'03_receiver/ActCD4_DE_all_genes.tsv.gz').set_index('canonical')
for _,r in genes.iterrows():
 if r.gene in de.index:assert np.allclose([r.logFC,r.P,r.BH],de.loc[r.gene,['logFC','PValue','FDR']].astype(float),rtol=1e-13,atol=0)
 else:assert pd.isna(r.logFC) and pd.isna(r.P) and pd.isna(r.BH)
checks['death_gene_membership_and_DE_verified']=len(genes)
for p in (OLD/'supplement').iterdir():assert sha(p)==sha(O/'supplement'/p.name)
checks['all_prior_supplement_files_preserved']=True
hs=BeautifulSoup((O/'report_v12.html').read_text(encoding='utf-8'),'html.parser');old=BeautifulSoup((OLD/'report_v11.html').read_text(encoding='utf-8'),'html.parser')
oldcells=Counter()
for s in old.find_all('section'):
 h=s.find(['h1','h2'])
 if h and h.get_text().startswith('4.2'):continue
 oldcells.update(c.get_text() for c in s.find_all(['td','th']))
newcells=Counter(c.get_text() for c in hs.find_all(['td','th']));assert all(newcells[k]>=v for k,v in oldcells.items());checks['unaffected_main_tables_preserved']=True
figures=[]
for i,img in enumerate(hs.find_all('img'),1):
 p=O/img['src'];assert p.is_file();cap=img.find_next_sibling('p',class_='caption').get_text();assert re.match(r'图'+str(i)+r'\s',cap)
 figures.append({'figure':i,'file':img['src'],'sha256':sha(p),'caption':cap})
pd.DataFrame(figures).to_csv(O/'provenance/figure_index.tsv',sep='\t',index=False,encoding='utf-8')
pdf=fitz.open(O/'silicosis_integrated_report_v12.pdf');assert len(pdf)==31;checks.update(pages=31,figures=len(figures),new_inferential_tests=0,visual_review='Changed pages and replacement figures inspected; unchanged figure layouts retained.',R01='presentation completed',R02='completed',R03='4.2 scatter format completed')
shutil.copytree(OLD/'provenance',O/'provenance/baseline_v11',dirs_exist_ok=True)
for f in ['ligand_activity_v4.R','composition_v4.py']:shutil.copy2(V4/'code'/f,O/'code'/f)
shutil.copy2(__file__,O/'code/finalize_v12.py')
(O/'provenance/final_validation.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2),encoding='utf-8')
(O/'README.txt').write_text('矽肺报告 v12：评分与NicheNet展示完善版\n\n主文件：silicosis_integrated_report_v12.pdf（31页、18幅图）。\n第5页：死亡程序评分公式、完整基因集、匹配对照及NicheNet输入定义。\n第6-7页：焦亡与凋亡分开展示，样本编号、效应量、P和BH均明确。\n第15页：4.2节改为四幅样本散点，左右分别为实际比例与统一比例。\n第19-21页：NicheNet来源排名、corrected AUPR以及配体-死亡程序基因先验权重。\n\nS18补充表：样本评分、组间统计、完整模块/对照成员、组成散点坐标和精确统计；炎症小体单独图为figures/S18_inflammasome_scores。\nS19补充表：完整NicheNet上/下调预测、来源排名、目标/背景成员、死亡基因差异及先验权重。\nS17候选选择依据及全部既有补充文件原样保留；历史来源见provenance/baseline_v11。\n\n本轮复用冻结统计，未新增推断检验。配体-死亡基因图为未改变的NicheNet先验映射，没有将死亡程序重新作为目标集合拟合。R01/R02及4.2格式修改均已完成对应展示范围验收。\n图号与文件名以provenance/figure_index.tsv为准。原报告版本继续保留。\n',encoding='utf-8')
rows=[]
for p in sorted(O.rglob('*')):
 if p.is_file() and p.name!='MANIFEST_SHA256.tsv' and not p.name.endswith('.inspect.ndjson'):rows.append({'file':p.relative_to(O).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)})
pd.DataFrame(rows).to_csv(O/'MANIFEST_SHA256.tsv',sep='\t',index=False,encoding='utf-8');zp=O.parent/'silicosis_presentation_complete_v12.zip'
with zipfile.ZipFile(zp,'w',zipfile.ZIP_DEFLATED) as z:
 for row in rows:z.write(O/row['file'],row['file'])
 z.write(O/'MANIFEST_SHA256.tsv','MANIFEST_SHA256.tsv')
with zipfile.ZipFile(zp) as z:assert z.testzip() is None
print(json.dumps(checks,ensure_ascii=False));print('ZIP_SHA256',sha(zp));print('PDF_SHA256',sha(O/'silicosis_integrated_report_v12.pdf'))
nums=[1,2,5,6,7,15,20,21,30];canvas=Image.new('RGB',(1160,2100),'#d8e0e5')
for j,n in enumerate(nums):
 p=T/f'page_{n:02d}.png';canvas.paste(ImageOps.contain(Image.open(p),(560,396)),((j%2)*580,(j//2)*420+20));ImageDraw.Draw(canvas).text(((j%2)*580+8,(j//2)*420+3),p.stem,fill='black')
canvas.save(T/'changed_pages.png')
