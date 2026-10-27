from pathlib import Path
import pandas as pd,numpy as np,json,shutil,requests,concurrent.futures,hashlib
O=Path(__file__).parent;H=O.parent;V=H/'lung_append_analysis_20261017_v18';W=H/'lung_feedback_revision_20261017_v19';T=O/'tables';T.mkdir(exist_ok=True);(O/'literature').mkdir(exist_ok=True);(O/'inputs').mkdir(exist_ok=True)
for n in ['lung_candidate_experiment_literature_20261017.txt','v19_independent_review_20261017.txt']:shutil.copy2(H/'Downloads'/n,O/'inputs'/n)
plan={'version':'v20','purpose':'客户完整报告重建，交付前限定核查','precedence':'候选实验文献优先于v19独立复核的实验排序建议','recommendation':['IL1R2主选/Treg功能','TTC39C探索备选/ActCD4存活与凋亡','SKIL状态备选/Treg与SnoN-SKI'],'modules':['免疫全景','CD4亚群','差异分析','GO/KEGG','GSEA','细胞通讯','TF网络','焦亡/凋亡评分','状态拟时序','NicheNet','候选与实验'],'figures_per_module':'1-2幅，主要模块目标2幅','statistics':'保留v19通过范围；新增模块另行数值核查，不宣称独立审阅已覆盖','trajectory':'复用冻结DPT，只在Naive/ActCD4内解释状态排序，不推断Treg谱系、不以细胞/bin充当重复','new_enrichment_if_needed':'使用已冻结v16完整DE，明确映射背景/原始P候选/方向/完整家族FDR；不改变候选入池'}
(O/'plan.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2),encoding='utf-8')
x=pd.read_csv(V/'tables/sample_selected_gene_expression.tsv.gz',sep='\t');x=x[x.gene.isin(['IL1R2','TTC39C','SKIL'])].copy();x['positive_fraction']=x.positive_cells/x.cells;x['multiUMI_fraction']=x.multiUMI_cells/x.cells;x['log1p_CP10k']=np.log1p(x.CP10k);x.to_csv(T/'candidate_sample_detection.tsv',sep='\t',index=False)
rows=[]
for (ds,var,sc,g),z in x.groupby(['dataset','variant','scope','gene']):
 z=z.copy(); c=z[z.group=='Saline'];d=z[z.group=='Silica'];r=dict(dataset=ds,variant=var,scope=sc,gene=g,samples=len(z),positive_samples=int((z.positive_cells>0).sum()),min_fraction=z.positive_fraction.min(),median_fraction=z.positive_fraction.median(),max_fraction=z.positive_fraction.max(),max_sample_UMI_share=z.UMI.max()/z.UMI.sum() if z.UMI.sum() else np.nan)
 if len(c) and len(d):
  effect=d.log1p_CP10k.mean()-c.log1p_CP10k.mean();loo=[d.log1p_CP10k.mean()-c.drop(j).log1p_CP10k.mean() for j in c.index if len(c)>1]+[d.drop(j).log1p_CP10k.mean()-c.log1p_CP10k.mean() for j in d.index if len(d)>1];r.update(descriptive_logCP10k_difference=effect,LOO_direction_fraction=np.mean(np.sign(loo)==np.sign(effect)) if loo else np.nan,LOO_min=min(loo) if loo else np.nan,LOO_max=max(loo) if loo else np.nan)
 rows.append(r)
pd.DataFrame(rows).to_csv(T/'candidate_detection_summary.tsv',sep='\t',index=False)
# Pre-delivery evidence register: conservative descriptions, figure/section locator and exact scope.
lit=[('IL1R2','PMC9438093','2022 JITC','小鼠Foxp3特异Il1r2缺失；MC38肿瘤','Fig4E-F；Fig2；Methods','减少肿瘤Treg克隆扩增及效应特征；不是矽肺Treg损伤实验','条件缺失模型有直接Treg先例；普通体外敲低效应仍待测'),('IL1R2','174005','2024 JCI Insight','Foxp3相关缺失，生发中心Tfr','Fig3D-L；IL-1阻断结果','初免GC反应增强并可被IL-1阻断逆转；不能等同全部Treg','特定Tfr功能及IL-1轴证据'),('IL1R2','PMC12533375','2025 Advanced Science','巨噬细胞/脓毒症；缺失和ENO1抑制','Fig2-4；Fig6','胞内IL1R2-ENO1与糖酵解/焦亡；Treg作用需另证','区分胞内、表面、可溶形式；RNA关联非酶活'),('IL1R2','PMC13214385','2026 Particle and Fibre Toxicology','矽肺上皮NLRP3/中性粒细胞','Il1r2表达相关结果','已有矽肺表达报道；不是CD4/Treg功能证明','不宣称IL1R2首次涉矽肺'),('TTC39C','PMC9615393','2022 BMC Pulmonary Medicine','人A549/NCI-H1299，shRNA','Fig3A-B；Cell transfection','敲低增加Annexin V阳性；与项目正向程序关联不直接同义','人肿瘤shRNA可行性不等于原代小鼠CD4已验证'),('SKIL','PMC8421263','2021 JEM','小鼠初始CD4，Arkadia缺失/CRISPR','Fig5A-E','双靶向Ski/SnoN可救援，单靶向不能；特定背景冗余','测SnoN蛋白与SKI；不把RNA上调等同蛋白增加'),('SKIL','40154177','2025 Int Immunopharmacol','小鼠及大鼠RLE-6TN；VC/MG132','原始论文摘要，未取得全文图级验证','ROS/SnoN蛋白降解及肺纤维化已有联系','广义肺纤维化联系非全新；不是Treg直接实验')]
urls={'174005':'https://insight.jci.org/articles/view/174005','40154177':'https://pubmed.ncbi.nlm.nih.gov/40154177/'}
ll=pd.DataFrame(lit,columns=['gene','id','source','context','locator','finding_and_boundary','feasibility']);ll['url']=ll.id.map(lambda s:urls.get(s,'https://pmc.ncbi.nlm.nih.gov/articles/'+s+'/'));ll.to_csv(T/'literature_evidence_register.tsv',sep='\t',index=False)
def fetch(row):
 try:
  r=requests.get(row.url,timeout=35);r.raise_for_status();p=O/'literature'/(row.id+'.html');p.write_text(r.text,encoding='utf-8');return {'id':row.id,'url':row.url,'status':r.status_code,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
 except Exception as e:return {'id':row.id,'url':row.url,'error':str(e)}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:res=list(ex.map(fetch,ll.itertuples()))
(O/'literature/download_manifest.json').write_text(json.dumps(res,ensure_ascii=False,indent=2),encoding='utf-8')
print(pd.DataFrame(rows).query("dataset=='GSE314943' and variant=='Raw' and ((gene=='IL1R2' and scope=='Treg') or (gene=='SKIL' and scope=='Treg') or (gene=='TTC39C' and scope=='ActCD4'))").to_string(index=False))