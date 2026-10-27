from pathlib import Path
import sys,json,hashlib,html,csv
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(sys.argv[1]);P=O.parents[1];V=P/'11_reports/customer_sequence_20261015_v2'
# Curated summaries are original paraphrases, not reproduced abstracts.
R=[]
def add(g,year,label,pmid,doi,model,finding,kind,limit,url=None):
 R.append(dict(evidence_id=f'L{len(R)+1:02d}',gene=g,year=year,study=label,PMID=pmid,DOI=doi,model=model,finding_zh=finding,evidence_type=kind,limit_zh=limit,url=url or ('https://pubmed.ncbi.nlm.nih.gov/'+pmid+'/' if pmid else 'https://doi.org/'+doi)))
add('AHNAK',2008,'CD4 T-cell calcium signaling','18191595','10.1016/j.immuni.2007.11.020','Ahnak1缺失小鼠及CD4 T细胞','缺失导致TCR刺激后的钙流、NFAT活化、IL-2和增殖降低；已建立CD4钙信号功能。','direct_gene_function','非矽肺；不能声称首次发现AHNAK调控CD4活化。')
add('AHNAK',2018,'Whole-lung Ahnak knockout','29724814','10.1158/1541-7786.MCR-17-0726','全身Ahnak缺失小鼠；肺上皮及巨噬细胞','缺失促进II型肺泡上皮增生及促肿瘤M2巨噬细胞状态。','direct_gene_function','非CD4特异干预；不能由CD4中上调推导全身抑制有益。')
add('AHNAK',2021,'Gm40600-Ahnak in CD4 responses','33987830','10.1111/imm.13365','Gm40600转基因/敲除小鼠及CD4细胞','Gm40600与Ahnak互作、影响其定位并促进钙信号和CD4反应。','mechanistic_interaction','验证重点为Gm40600；并非矽肺、自噬或焦亡机制。')
add('MDFIC',2000,'Molecular cloning and viral transcription','10671520','','HIC/MDFIC克隆及HTLV-I/HIV表达体系','早期已研究MDFIC对病毒转录表达的调节。','molecular_function','病毒转录体系不能替代生理CD4疾病机制。')
add('MDFIC',2017,'Glucocorticoid receptor cross-talk','28223352','10.1074/jbc.M116.758888','GR表达细胞与MDFIC过表达/耗减','MDFIC与GR相互作用并改变其磷酸化及转录反应；糖皮质激素反向抑制MDFIC表达。','direct_gene_function','组织和刺激依赖；不是矽肺CD4验证。')
add('MDFIC',2019,'miR-146b-3p/MDFIC in myoblasts','31261950','','成肌细胞；miRNA/MDFIC/PI3K-AKT轴','研究涉及增殖、分化及凋亡，不能把MDFIC与细胞死亡的联系称为全新。','direct_gene_function','非免疫细胞；不证明CD4焦亡或自噬。')
add('MDFIC',2020,'MDFI versus MDFIC in colorectal cancer','32457453','','HCT116结直肠癌细胞及JMJD1A互作','MDFI和MDFIC呈不同作用；MDFIC抑制HCT116生长。','direct_gene_function','MDFI与MDFIC不是同一基因；效应不能跨细胞类型照搬。')
add('MDFIC',2022,'Human lymphatic disease and mouse model','35235341','10.1126/scitranslmed.abm4869','六个家系、MDFIC截断突变小鼠及淋巴内皮细胞','双等位致病变异联系淋巴管异常；实验涉及整合素β1活化、细胞迁移和瓣膜形成。','human_genetics_and_function','淋巴内皮并非T淋巴细胞。')
add('MDFIC',2024,'Standigm ASK IPF target study','38349059','10.1093/bib/bbae035','计算筛选；L132及HEK293T相关TGFβ/EMT实验','MDFIC已被提出为IPF候选靶点；TGFβ诱导表达，敲低逆转EMT相关读数及迁移。','fibrosis_related_function','L132被Cellosaurus标记为HeLa衍生污染细胞；降低肺上皮外推强度，但不消除既往发表事实。')
add('MDFIC',2026,'IFNalpha8-memory CD4 HIV study','','10.21203/rs.3.rs-10466271/v1','静息记忆CD4与HIV再激活；蛋白组学','2026-08-05预印本报告IFNα8响应中MDFIC下降，处于STAT1相关病毒转录背景。','preprint_expression','未同行评审；不是MDFIC特异干预确证，更非矽肺证据。')
add('TTC39C',2019,'Muscle atrophy and signaling','31188487','10.1002/jcp.28950','小鼠神经源性萎缩及培养肌细胞','已研究表达、定位及增减表达；异常表达影响肌分化、ERK1/2和Hedgehog相关读数。','direct_gene_function','过表达和敲低均可损害分化；不宜假设单向通用作用。')
add('TTC39C',2022,'Lung adenocarcinoma perturbation','36303158','10.1186/s12890-022-02173-x','TCGA与A549/NCI-H1299细胞','敲低抑制增殖、迁移和克隆形成并增加凋亡；另有代谢及p53通路富集。','direct_gene_function','富集不等于已证实直接通路；不是CD4或矽肺。')
add('TTC39C',2026,'Double-negative T-cell GC signature','41526926','10.1186/s12935-025-04080-7','胃癌单细胞/机器学习、队列及表达验证','TTC39C进入双阴性T细胞相关四基因诊断模型。','signature_and_expression','双阴性T不是CD4 T；标志物模型不能作为CD4特异机制。')
add('ST6GALNAC3',2005,'Enzyme cloning and substrate specificity','16169874','10.1093/jb/mvi124','人ST6GalNAc III克隆及酶学','已确立唾液酸转移酶相关底物特性，涉及糖脂及部分O-糖链。','biochemical_function','与ST6GAL1、ST6GALNAC1/2区分。')
add('ST6GALNAC3',2018,'Prostate cancer methylation biomarker','29465788','10.1002/1878-0261.12183','前列腺癌组织及液体活检','ST6GALNAC3启动子高甲基化已被作为生物标志物研究。','biomarker','不是免疫细胞功能实验。')
add('ST6GALNAC3',2020,'IPF bulk expression/network analysis','32596346','10.1155/2020/6537462','GEO肺组织表达及网络分析','主文表1将ST6GALNAC3列在IPF下调基因列表。','bulk_expression_association','不是基因特异功能验证；不能声称首次与肺纤维化相关。')
add('ST6GALNAC3',2021,'Rora perturbation in activated CD4','','10.1371/journal.pone.0251233','Rora条件敲除/过表达；活化CD4转录组','已在Rora扰动背景中报道St6galnac3差异表达；存在组织及实验条件依赖。','upstream_perturbation_association','未直接干预St6galnac3；不能声称首次发现其在CD4中表达或变化。',url='https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0251233')
add('ST6GALNAC3',2021,'A549 growth and TFR1 expression','34461439','10.1016/j.bbrc.2021.08.055','A549敲低及TFR1救援','敲低降低TFR1 mRNA稳定性、蛋白量及细胞生长；TFR1过表达部分救援。','direct_gene_function','文中TFR1唾液酸化未显著改变；不能改写成已证明直接糖基化TFR1，更不能外推为CD4自噬机制。')
e=pd.DataFrame(R);e.to_csv(O/'literature/curated_evidence.tsv',sep='\t',index=False,encoding='utf-8')
nov=pd.DataFrame([
 dict(gene='ST6GALNAC3',known='酶学、肿瘤/TFR1功能、IPF表达列表、Rora扰动后CD4表达均已有研究',space='硅尘背景下ActCD4中的基因特异功能及相关糖脂/膜蛋白表达变化仍有明确待检验空间',priority='优先探索性验证',reason='小鼠状态和外部方向较一致；创新应落在特定情境和因果功能，而非表达发现',avoid='不能声称从未在CD4或肺纤维化研究；不预设自噬或直接TFR1糖基化机制'),
 dict(gene='AHNAK',known='CD4钙流/NFAT/活化机制已较成熟；肺上皮及巨噬细胞作用也已有实验',space='长期硅尘背景中特定NaiveCD4状态的适应性/功能作用',priority='机制基础明确的并行备选',reason='鼠内证据可复核，但人类NaiveCD4方向相反；泛化机制创新较小',avoid='不能称首次发现AHNAK调控T细胞，也不由上调推出抑制有益'),
 dict(gene='TTC39C',known='肌细胞信号、肺癌凋亡/代谢、2026年双阴性T细胞相关模型已有研究',space='晚期硅尘ActCD4中的独立功能及病程依赖',priority='阶段依赖候选，先解决可重复性',reason='D7部分口径敏感，但D56下降仍在；20周与早期不能直接合并解释',avoid='不能把肿瘤或双阴性T细胞证据替代CD4机制'),
 dict(gene='MDFIC',known='GR互作、淋巴内皮、人类淋巴管病、2024年IPF/EMT靶点研究及2026年CD4预印本',space='与已发表上皮/内皮工作区分的硅尘CD4功能',priority='当前降低优先级',reason='已有肺纤维化功能主张；外部方向及拟时序敏感性问题并存',avoid='既往L132模型局限不能用来声称该基因从未被研究')])
nov.to_csv(O/'literature/novelty_priority_assessment.tsv',sep='\t',index=False,encoding='utf-8')
exceptions=pd.DataFrame([
 ['AHNAK2','AHNAK','2022 AHNAK2肺纤维化论文不是当前AHNAK候选的直接证据','https://pubmed.ncbi.nlm.nih.gov/35882062/'],
 ['TTC39B/TTC39A','TTC39C','同家族脂代谢等研究不能归给TTC39C','https://www.ncbi.nlm.nih.gov/gene/125488'],
 ['MDFI/MDFIC2/MDFIC-7','MDFIC','相关基因及lncRNA名称不能混用','https://pubmed.ncbi.nlm.nih.gov/32457453/'],
 ['ST6GAL1/ST6GALNAC1/ST6GALNAC2','ST6GALNAC3','唾液酸化酶家族功能不能相互替代','https://www.ncbi.nlm.nih.gov/gene/256435'],
 ['L132','MDFIC prior study','Cellosaurus记录为HeLa衍生污染细胞；作为先前发表事实保留，模型外推降级','https://www.cellosaurus.org/CVCL_1908']
],columns=['item','candidate','interpretation','source'])
exceptions.to_csv(O/'literature/identity_and_model_cautions.tsv',sep='\t',index=False,encoding='utf-8')
# Save explicit cohort-time evidence independently of original day=-1 analysis fields.
time=pd.DataFrame([
 ['GSE314943','20 weeks (~140 days)','3 saline + 3 silica mice/libraries','Original publication Fig.2A legend; GEO series n=3/group','https://www.nature.com/articles/s41467-026-70193-w'],
 ['GSE183682','D7 and D56','one mouse/library per condition/time','GSM5567502 treatment protocol explicitly states four mice for sequencing','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5567502'],
 ['GSE274372','D0/D7/D28/D56','3 libraries/time; baseline contrasts','Measured-time series lacks contemporaneous late saline controls','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274372']
],columns=['dataset','sampling_time','replicates','evidence','source'])
time.to_csv(O/'direction/cohort_time_evidence.tsv',sep='\t',index=False,encoding='utf-8')
# Stable diagnostic issue IDs, bilingual for internal review; no automatic external messages.
issues=[
 ['D01','高／解释边界','已澄清设计差异，生物学原因未证实','GSE314943 Fig.2A; direction/cohort_time_evidence.tsv','主队列20周；外部D7/D56','不同阶段不能当作同条件复制，也不能直接证明阶段反转','解释中分开阶段；同一硅尘模型设置D7/D56/20周及各自对照','需独立生物学重复与相同CD4状态定义后，才能确认阶段依赖','Discovery is 20 weeks; external cohorts are D7/D56. Design mismatch is established, but phase-dependent reversal is unproven. Validate with replicated matched-time controls and harmonized CD4 states.'],
 ['D02','中／候选稳健性','部分解释','direction/direction_by_metric_and_filter.tsv; gene=TTC39C, scope=ActCD4','D7细胞均值差-0.070427 CP10k；伪批量差+4.604780 CPM；D56两者均负','D7方向取决于估计量；D56矛盾仍在','同时呈现检出率、阳性细胞强度及伪批量；不选择有利口径覆盖主结果','同条件重复中检出与表达强度均可复现，或明确证实不同估计量反映不同变化','TTC39C D7 is estimand-sensitive; D56 remains negative. Preserve both summaries and validate detection frequency and expression intensity with independent matched samples.'],
 ['D03','高／选靶风险','未解决生物学方向矛盾','direction/ActCD4_depth_standardization.tsv; gene=MDFIC','D7/D56在均值、伪批量、BG1及身份/深度分层后仍下降；D7病例仅1/95阳性细胞','不能用归一化或选取细胞子集声明复制成功','降低机制实验优先级；先补同阶段独立样本及RNA/蛋白检出','补充有效重复后确认方向，或明确限制为20周模型情境','MDFIC stays negative under normalization, background, identity and depth diagnostics. The D7 case has one detected cell out of 95. Deprioritize until independent matched-stage expression validation.'],
 ['D04','高／人类外推','未解决跨物种/组织差异','direction/human_scope_group_directions.tsv; direction/human_leave_one_donor_out.tsv','AHNAK人类NaiveCD4差-3.767899 CP10k；五次删一供者均负；Core组成变化不能解释其下降','不能写人类方向验证成功；也不能改选MemCD4正向结果替代','分开肺组织小鼠与BALF人类；预先固定匹配状态和对照定义','独立人类同状态样本中重现方向，或明确仅作小鼠情境候选','Human NaiveCD4 AHNAK remains negative in all five leave-one-donor analyses. Do not substitute the positive MemCD4 result. Validate matched tissue/state and retain the mouse-specific interpretation meanwhile.']]
cols=['id','severity','status','file_and_location','evidence','impact','action','acceptance','english_feedback']
iss=pd.DataFrame(issues,columns=cols);iss.to_csv(O/'direction/findings_zh_en.tsv',sep='\t',index=False,encoding='utf-8')
# Diagnostic figures with explicit units and no cell-level inference.
con=pd.read_csv(O/'direction/direction_by_metric_and_filter.tsv',sep='\t');dec=pd.read_csv(O/'direction/detection_intensity_decomposition.tsv',sep='\t')
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axs=plt.subplots(1,2,figsize=(13,5));fig.subplots_adjust(left=.10,right=.85,bottom=.23,top=.84,wspace=.42)
metrics=['mean_CP10k','pseudobulk_CPM','mean_logCP10k','detection_fraction','positive_mean_CP10k']
for ax,gene in zip(axs,['MDFIC','TTC39C']):
 a=con[(con.gene==gene)&(con.scope=='ActCD4')&(con['filter']=='All')]
 mat=np.array([[a[(a.dataset==ds)&(a.day==day)&(a.metric==m)].log2_case_control.iloc[0] for m in metrics] for ds,day in [('GSE314943',-1),('GSE183682',7),('GSE183682',56)]])
 im=ax.imshow(mat,cmap='RdBu_r',vmin=-2,vmax=2,aspect='auto')
 for i in range(3):
  for j in range(5):ax.text(j,i,f'{mat[i,j]:+.2f}',ha='center',va='center',color='white' if abs(mat[i,j])>1.3 else 'black')
 ax.set(title=gene+' / ActCD4',xticks=range(5),xticklabels=['Cell mean\nCP10k','Pseudobulk\nCPM','Mean\nlogCP10k','Detection\nfraction','Positive-cell\nCP10k'],yticks=range(3),yticklabels=['Discovery / 20 wk','External / D7','External / D56'])
cax=fig.add_axes([.90,.26,.02,.5]);fig.colorbar(im,cax=cax,label='log2(case/control); color limited to +/-2')
fig.suptitle('Direction depends on estimand and cohort; no cells used as replicates',fontsize=15);fig.savefig(O/'figures/direction_estimand_grid.png',dpi=170,bbox_inches='tight');plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(12,5))
for ax,gene in zip(axs,['MDFIC','TTC39C']):
 a=dec[(dec.gene==gene)&(dec.scope=='ActCD4')&(dec['filter']=='All')].set_index(['dataset','day']).loc[[('GSE314943',-1),('GSE183682',7),('GSE183682',56)]]
 x=np.arange(3);ax.bar(x-.17,a.detection_component,width=.32,label='Detection-fraction component',color='#357b9b');ax.bar(x+.17,a.positive_intensity_component,width=.32,label='Positive-expression component',color='#d59a45');ax.scatter(x,a.delta,marker='D',c='black',zorder=4,label='Observed total difference');ax.axhline(0,color='grey',lw=.8);ax.set(title=gene+' / ActCD4',xticks=x,xticklabels=['Discovery\n20 weeks','External\nD7','External\nD56'],ylabel='Case - control mean CP10k')
handles,labels=axs[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=3,frameon=False);fig.suptitle('Exact arithmetic decomposition; components do not identify biological causes',fontsize=14);fig.tight_layout(rect=[0,.12,1,.92]);fig.savefig(O/'figures/detection_intensity_components.png',dpi=170,bbox_inches='tight');plt.close(fig)
h=pd.read_csv(V/'05_supporting_evidence/GSE174725_candidate_donor_expression.tsv',sep='\t');fig,axs=plt.subplots(1,4,figsize=(14,4.5))
for ax,scope in zip(axs,['Core','NaiveCD4','MemCD4','Treg']):
 a=h[(h.gene=='AHNAK')&(h.scope==scope)&(h['filter']=='ReviewedAll')].dropna(subset=['CP10k']);ns=[]
 for j,group in enumerate(['Exposed','Silicosis']):
  b=a[a.group==group];ns.append(len(b));ax.scatter(j+np.linspace(-.07,.07,len(b)),b.CP10k,c=['#357b9b','#c55345'][j],s=35);ax.hlines(b.CP10k.mean(),j-.15,j+.15,color='black')
 ax.set(title=f'{scope}: n={ns[0]} vs {ns[1]}',xticks=[0,1],xticklabels=['Exposed','Silicosis'],ylabel='Donor mean CP10k')
fig.suptitle('Human AHNAK: separate scopes, retain every nonempty donor',fontsize=15);fig.tight_layout(rect=[0,0,1,.91]);fig.savefig(O/'figures/human_AHNAK_scope_points.png',dpi=170,bbox_inches='tight');plt.close(fig)
# Validate decompositions and immutable v2 release.
for r in dec.itertuples():
 if np.isfinite(r.detection_component+r.positive_intensity_component):assert abs(r.detection_component+r.positive_intensity_component-r.delta)<1e-8
hc=pd.read_csv(O/'direction/human_Core_composition_decomposition.tsv',sep='\t');hr=pd.read_csv(O/'direction/human_scope_group_directions.tsv',sep='\t')
for (f,g),a in hc.groupby(['filter','gene']):
 original=hr[(hr['filter']==f)&(hr.gene==g)&(hr.scope=='Core')].difference.iloc[0];assert abs(a.composition_component.sum()+a.within_state_component.sum()-original)<1e-8
v2manifest=pd.read_csv(V/'OUTPUT_MANIFEST_SHA256.tsv',sep='\t',dtype=str)
for r in v2manifest.itertuples():assert hashlib.sha256((V/r.path).read_bytes()).hexdigest()==r.sha256
# Internal Chinese research memo; this is new research, not a customer translation of v2.
def table(d):return d.to_html(index=False,escape=True)
source_rows=[]
for r in R:
 source_rows.append('<tr>'+''.join('<td>'+html.escape(str(r[k]))+'</td>' for k in ['evidence_id','gene','year','model','finding_zh','evidence_type','limit_zh'])+'<td><a href="'+html.escape(r['url'],quote=True)+'">'+html.escape(r['study'])+'</a></td></tr>')
parts=['<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>候选创新性与方向矛盾：内部研究备忘录</title><style>body{font:16px "Microsoft YaHei",Arial,sans-serif;line-height:1.65;max-width:1400px;margin:30px auto;padding:0 22px;color:#183044}h1,h2,h3{color:#17496a}table{border-collapse:collapse;font-size:13px;display:block;overflow-x:auto;margin:14px 0}td,th{border:1px solid #c8d5df;padding:8px;vertical-align:top}img{max-width:100%;height:auto}.note{background:#edf3f8;padding:14px;border-left:4px solid #43799a}</style><body>',
'<h1>候选创新性与方向矛盾：内部研究备忘录</h1><p>检索截至 2026-10-15；对象：AHNAK、MDFIC、TTC39C、ST6GALNAC3。v2选靶规则、原P/BH、旧报告和数据均保留，不通过调整口径强行保留候选。</p>',
'<p class="note">当前建议：ST6GALNAC3优先进行情境特异的探索性验证；AHNAK作为机制基础较强的并行备选；TTC39C保留为需要匹配阶段验证的候选；MDFIC降低当前优先级。四基因都不能称为“从未被研究”，本轮也未建立其自噬/焦亡因果机制。</p>',
'<h2>一、创新应落实到“硅尘—具体CD4状态—因果功能”</h2>',table(nov),
'<p>AHNAK的CD4钙流/NFAT机制已有成熟基础。ST6GALNAC3虽更具可探索空间，但已有IPF下调列表、Rora扰动下CD4表达及A549功能研究。TTC39C已有肺癌凋亡/代谢工作和2026年双阴性T细胞相关模型。MDFIC已有2024年IPF靶点/EMT研究，不能把再次得到肺纤维化关联当作主要创新。</p>',
'<p>本轮未确认四候选在“矽肺CD4特定状态中，经基因特异干预证明焦亡/自噬功能”的直接证据。这是可检验研究空白，不是绝对首创证明。ST6GALNAC3—TFR1文献明确报告TFR1唾液酸化未显著改变，不能把表达量及mRNA稳定性证据改写成直接糖基化机制。Rora扰动下的差异表达也不等于已证实直接、单向调控。</p>',
'<h3>已核读的主要研究</h3><table><tr><th>ID</th><th>基因</th><th>年份</th><th>模型</th><th>已做工作</th><th>证据属性</th><th>外推边界</th><th>原始来源</th></tr>'+''.join(source_rows)+'</table>',
'<h3>名称和模型核查</h3>',table(exceptions),
'<p>MDFIC/IPF研究使用的L132被<a href="https://www.cellosaurus.org/CVCL_1908">Cellosaurus</a>列为HeLa衍生污染细胞。这里将该工作作为已发表的细胞实验保留，同时降低其正常肺上皮外推强度；该模型问题不能用于宣称该方向无人研究。</p>',
'<h2>二、方向矛盾：已排查什么，仍缺什么</h2>',
'<p><b>D01：比较的病程并不相同。</b>发现队列原论文图2A明确记载20周采样、每组3只小鼠；外部队列为D7/D56，每条件/时间仅1只小鼠。20周与D7/D56必须分开解释。它们提示阶段依赖这一待检验解释，但没有同模型跨阶段重复，不能证明表达随病程反转。<a href="https://www.nature.com/articles/s41467-026-70193-w">发现队列原文</a>；<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5567502">外部样本原始元数据</a>。</p>',table(time),
'<p><b>D02：TTC39C的D7矛盾部分来自估计量。</b>外部D7细胞均值从0.554209到0.483782 CP10k，差-0.070427；伪批量从55.082458到59.687239 CPM，差+4.604780。检出比例从10/60降到7/95，而阳性细胞归一化强度增加。伪批量按总文库计数归一化，细胞均值对每个细胞等权，不能互换。D56两种主要口径仍下降，不能选D7有利口径宣布一致。</p>',
'<p><b>D03：MDFIC矛盾尚未消除。</b>D7病例仅1/95细胞检出、总目标UMI为1；对照5/60、UMI为5。D56对照18/135阳性，病例10/108阳性。下降主要由检出比例项驱动；阳性细胞强度并未同步明显下降。均值、伪批量、背景过滤、身份分层以及深度分层均没有消除其外部负向结果。稀疏检测和单样本设计使生物学原因仍不能确定。</p>',
'<img src="figures/direction_estimand_grid.png"><img src="figures/detection_intensity_components.png">',
'<p>外部D7病例平均文库为1234.5 UMI/细胞，对照3025.8；D56两组深度较接近。对文库深度作等样本权重的三分位分层标准化后，MDFIC D7/D56仍分别为-0.181963/-0.141121 CP10k；TTC39C仍为-0.054457/-0.087593。四个身份面板的分层结果全部保留，仅TTC39C D7的NaiveIdentity分层转正，不选择此单一规格覆盖其他结果。分层是描述性标准化，不是细胞层面的生物学重复检验。</p>',
'<p><b>D04：AHNAK的人类反向不是仅由Core组成造成。</b>人类NaiveCD4差-3.767899 CP10k，原精确P=0.2、BH=1；五次删一供者均保持负向。Core总差-2.062802，可分解为组成项约+0.214171和亚群内项约-2.276974。组成变化实际上抵消一部分下降。MemCD4为正的结果不能替代预定NaiveCD4验证，且其删一病例结果可变号。鼠肺与人BALF、疾病阶段及暴露对照不同，应保留为未解决的外推问题。</p>',
'<img src="figures/human_AHNAK_scope_points.png">',table(iss[['id','severity','status','evidence','action','acceptance']]),
'<h2>三、下一步如何衔接基础研究</h2><p>先验证ST6GALNAC3在预先定义ActCD4中的表达、蛋白/酶读数与病例方向；功能研究应区分糖脂/膜蛋白处理、增殖存活与自噬通量，TFR1表达/铁稳态仅作为可检验分支。AHNAK可围绕已知钙信号框架做硅尘情境和细胞状态特异验证，但不能预设抑制有益。TTC39C/MDFIC应先在同一模型的D7、D56、20周及对应对照中，用独立生物学重复检验检出比例与表达强度，再决定是否进入机制干预。这里的时间点来自真实队列差异，不是按目标基因结果挑选阈值。</p>',
'<h2>检索与复核范围</h2><p>保留33条Europe PMC可复核检索式，覆盖规范名、别名、题名摘要及全文；另核读PubMed、期刊原文和权威模型记录。每条自动检索最多返回100条；题名摘要零命中、全文偶然提及、同家族基因、综述及预印本均分别处理，命中数不是该基因已验证机制的数量。部分全文接口不可获取的命中未作为已证实功能证据。</p><p>方向诊断覆盖320行样本表达汇总、480行估计量/背景对比、48条身份分层及12条深度分层记录；保留所有非空样本，没有20/10细胞样本门槛，没有新增外部单样本推断P值。检测/强度及组成分解为精确代数分解，不是因果归因。所有v2发布文件哈希已核对未变。</p><p>完整表：literature/curated_evidence.tsv、novelty_priority_assessment.tsv、search_log.tsv、supplemental_search_log.tsv；direction/findings_zh_en.tsv、direction_by_metric_and_filter.tsv及各分解/敏感性表。全文缓存仅供内部核读，不随结果包发布。</p></body></html>']
(O/'research_memo_zh.html').write_text('\n'.join(parts),encoding='utf-8')
qa={'status':'PASS','curated_studies':len(R),'search_queries':33,'v2_manifest_files_unchanged':len(v2manifest),'direction_decomposition_checked':True,'human_composition_checked':True,'diagnostic_ids':['D01','D02','D03','D04'],'candidate_rule_modified':False,'external_cell_level_P_tests':0}
(O/'provenance/report_validation.json').write_text(json.dumps(qa,indent=2),encoding='utf-8');print(json.dumps(qa))