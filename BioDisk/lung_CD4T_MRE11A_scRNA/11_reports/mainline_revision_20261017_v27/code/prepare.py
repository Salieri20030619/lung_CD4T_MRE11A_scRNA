from pathlib import Path
import json,shutil,tokenize,io,ast,pandas as pd,hashlib
H=Path(r'C:\Users\31937');O=H/'lung_mainline_revision_20261017_v27';O.mkdir(exist_ok=True)
for d in ['figures','plot_data','evidence','supplement']: (O/d).mkdir(exist_ok=True)
B=H/'lung_experiment_plan_20261017_v26'
reps=json.loads((B/'editorial_replacements.json').read_text(encoding='utf-8'))
def clean(x):
 for a,b in reps:x=x.replace(a,b)
 for a,b in [('同状态巨噬入口','同一T细胞状态的巨噬来源信号依据'),('巨噬来源入口','巨噬来源信号依据'),('巨噬入口','巨噬来源信号'),('历史匹配对照评分','匹配对照基因评分'),('既有K15','K15'),('min10敏感性为11样本','每样本至少10个TLike细胞的敏感性分析纳入11个样本'),('Raw/BG1','完整细胞集/低背景负担细胞子集'),('上游未补齐','尚无明确上游配体，先检验巨噬刺激下的候选效应'),('p53/凋亡链未补齐','当前数据不支持疾病调整后的凋亡程序关联'),('不作为第四个统一硬门槛','用于确定预实验优先顺序'),('不产生二元创新性P值','结合机制与细胞背景判断研究空间'),('不是52基因的继续硬筛子集','来自完整90基因评价集合'),('全部统计值来自原分析。',''),('图例位于独立区域，颜色对应同名通路；',''),('；该展示子集不用于重新排序候选',''),('保留原8条代表路线；',''),('DPT排序、统计及不包含Treg的范围均保持一致。','分析范围为Naive与ActCD4。')]:x=x.replace(a,b)
 return x
source=(B/'report.py').read_text(encoding='utf-8');tok=[]
for t in tokenize.generate_tokens(io.StringIO(source).readline):
 if t.type==tokenize.STRING:
  try:
   v=ast.literal_eval(t.string)
   if isinstance(v,str):t=t._replace(string=repr(clean(v)))
  except (ValueError,SyntaxError):pass
 tok.append(t)
s=tokenize.untokenize(tok)
s=s.replace('矽肺_CD4候选研究_实验方案修订_v26_20261017','矽肺_CD4候选研究_主线修订_v27_20261017').replace('矽肺免疫全景_CD4候选与验证_报告_v26.pdf','矽肺免疫全景_CD4候选与验证_报告_v27.pdf')
(O/'editorial_replacements.json').write_text('[]',encoding='utf-8')
manifest=json.loads((B/'figure_manifest.json').read_text(encoding='utf-8'))
for m in manifest:
 for k in ['title','module','caption']:m[k]=clean(m[k])
 if m['id'] in ['F09','F10']:m['caption']=m['caption'].replace('采用enrichplot::gseaplot2绘制','展示')
 if m['id'] in ['F11','F12']:m['module']='细胞通讯预测'
 if m['id']=='F13':m['caption']='12个TF去除候选及程序重叠成员后的样本平均AUCell，行标准化仅用于展示。每状态每组3个独立样本；图中展示描述性活性分布。'
 if m['id'] in ['F15','F16']:m['caption']='样本汇总计数评分，移除9个候选与程序的重叠成员。每点为独立样本，3+3；双侧组标签排列20种，FDR基于52个程序-状态检验。PyroExec、ApoExec、MitoISR分别表示焦亡相关、凋亡相关及线粒体应激程序。'
 if m['id']=='F20':
  m['file']=str(O/'figures/F20.png');m['title']='IL1R2诱饵结合与SKIL响应先验提供不同验证线索';m['caption']='蓝色为配体-靶先验及靶排名，橙色为IL-1与IL1R2诱饵结合；灰色为本展示范围未列入支持路线。NA表示TTC39C未达到对应联合差异输入条件，未纳入评价。下方紫色单列H2-DM抗原呈递相关记录，不作为普通表面配体-CD4信号推荐。排名为模型指标。'
 if m['id'] in ['F21','F22']:m['file']=str(O/'figures'/(m['id']+'.png'))
(O/'figure_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8');shutil.copy2(B/'extra_manifest.json',O/'extra_manifest.json')
# Explicit paragraph replacements after materializing previous editorial cleanup.
def rep(a,b):
 global s
 if a not in s:raise RuntimeError('missing replacement '+a[:80])
 s=s.replace(a,b)
rep("fn=H/'lung_cohort_bridge_20261017_v25/figures'/(id+'.png')","fn=(O/'figures' if id=='S_C' else H/'lung_cohort_bridge_20261017_v25/figures')/(id+'.png')")
rep('差异分析前的筛选框架：三项证据交集与后续分层','候选基因筛选策略：三项证据交集与综合评价')
# Rewrite opening to question and result-driven narrative, preserving all key quantitative evidence.
a=s.index("p('研究目的','h2')");b=s.index("page();p('研究设计与数据范围'",a)
opening="""p('研究问题与数据','h2');p('围绕矽肺背景下巨噬细胞如何影响T细胞表型，以及这种变化是否伴随线粒体损伤和焦亡开展研究。分析材料为小鼠公开单细胞数据：主队列GSE314943为3个对照与3个矽尘样本、37,155个全肺细胞，其中CD4细胞2,311个；与GSE183682的4个样本构成分层联合模型，另以GSE274372的12个样本评价时间背景。')
p('矽尘条件下T细胞出现什么变化','h2');p('Naive比例下降、ActCD4及Treg比例上升，组成检验P分别为0.10、0.20、0.10。Treg中IL1R2和SKIL呈疾病相关上调；功能富集显示翻译、免疫响应等程序变化。组成变化方向与表达差异分别构成表型和分子证据。')
p('巨噬信号与候选在哪些环节相联系','h2');p('通讯预测与NicheNet提供巨噬来源信号及T细胞响应候选：IL-1与Treg的IL1R2呈诱饵结合关系；SKIL具备趋化因子等配体的靶先验；TTC39C主要依据凋亡与线粒体应激程序关联。主队列ActCD4、Treg焦亡相关评分P为0.30、0.20，当前未检出显著组间差异。')
tab([['候选与切入环节','主要数据依据','首轮主要终点'],['IL1R2 / Treg表型，优先','联合上调logFC=3.290，P=1.66×10^-6，FDR=0.000662','Treg抑制功能'],['TTC39C / 损伤探索','ApoExec rho=0.943，P=0.00347，FDR=0.0556；联合DE P=0.163','存活及凋亡反应'],['SKIL / 状态调控备选','TAct rho=0.965，P=0.00208，FDR=0.0744','Treg状态稳定性']],[115,250,140])
p('首轮实验连接哪些环节','h2');p('以巨噬刺激×T细胞候选干预的四组设计，检验候选是否改变巨噬刺激引起的表型响应；同步观察膜电位、存活及独立焦亡执行读数。由结果决定候选推进、上游定位与救援，检验表型变化与线粒体、死亡终点之间的联系。','small')
"""
s=s[:a]+opening+s[b:]
rep('疾病调整rho=0.897、P=0.0278','疾病调整rho=0.897、P=0.0278、FDR=0.8681')
s=s.replace(";p('数据索引：'+m['id']+'；图值、来源表及脚本见配套材料索引。','small')",'')
rep("if m['id']=='F21':screening_methods()","if m['id'] in ['F21','F22']:continue")
# Detailed screening moves to appendix; retain user-requested compact overlap before DE.
# Replace results dictionary entries via an explicit update prior to figure rendering.
pos=s.index('for m in manifest:')
updates={
'F10':'ActCD4的抗原加工与呈递集合正向富集（NES=1.903，P=0.000247，FDR=0.0174），核糖体集合负向富集（NES=-1.965，P=5.81×10^-7，FDR=0.000164）。这些方向将活化状态与免疫响应、翻译程序的变化联系起来。',
'F13':'Treg中NFKB1的样本平均AUCell由对照0.0600升至矽尘0.0779，FOXP1由0.1539降至0.1354。这里呈现调控程序的描述性差异，为活化与状态变化提供背景，不新增显著性判定。',
'F14':'代表网络中JUNB连接GADD45B与NFKBIA，FOS连接JUNB与GADD45B，形成应激与活化相关靶网络线索。各边为SCENIC调控先验，后续可据此选择响应读数。',
'F15':'ActCD4和Treg的焦亡相关评分均呈上升方向，精确P分别为0.30和0.20；一般炎症评分的组间P均为0.90。本队列未检出显著组间差异，焦亡作为后续实验的独立观察终点。',
'F16':'ActCD4的凋亡相关程序ApoExec差值为+0.0404、P=0.10、FDR=0.274；Treg差值为+0.0352、P=0.20、FDR=0.433。两状态方向均向上，未达到名义显著；后续结合候选关联检验具体损伤环节。',
'F19':'联合疾病上调目标集的校正AUPR排序中，Treg以ICAM1（0.0552）和TNF（0.0547）居前；ActCD4以WNT4（0.0241）和IL10（0.0237）居前。排序反映各状态响应集合的模型解释度，候选基因的具体连接见下一图。',
'F20':'IL1R2的IL-1诱饵结合与SKIL的CXCL10、CCL2、SPP1靶先验属于不同证据类型。TTC39C未纳入本网络路线评价；H2-DM相关记录单列为抗原呈递信息。',
'F21':'三项证据交集包含55个基因，与42个研究问题相关补充候选合并，重叠7个，得到90个综合评价基因。',
'F22':'90个基因形成52个计算支持、37个证据待补及1个未达到名义差异条件的分层。9个重点基因来自完整90基因集合，其中5个计算支持、4个证据待补；3个预实验建议结合功能关联、文献及可行性确定。'}
s=s[:pos]+'results.update('+repr(updates)+')\n'+s[pos:]
# Methods: explicit definitions, no internal version references.
needle="page();p('评分、关联与状态轨迹方法','h')"
rep(needle,"p('GSEA曲线使用enrichplot::gseaplot2绘制。Raw表示完整细胞集；BG1表示background_fraction≤1%的低背景负担细胞子集。该比例为预定义上皮、髓系及免疫球蛋白标志基因UMI占细胞完整文库UMI的比例，不是环境RNA污染率估计。完整标志集合与计算定义见补充材料。','small')\n"+needle)
# One compact main experiment page; detailed six-page plan retained as appendix.
a=s.index("page();p('基础实验推荐方案｜总体路线与实施前提'");b=s.index("page();p('讨论：证据边界与研究空间'",a);detail=s[a:b]
detail=detail.replace('基础实验推荐方案｜总体路线与实施前提','附录B｜实验实施详案与参数确认')
detail=detail.replace('FOXP3/CD25表型、活细胞绝对数量、存活及凋亡相关读数；区分Treg数量减少与单位活细胞功能改变','FOXP3/CD25、活细胞数、膜电位与死亡读数；平行检测目标Treg的gasdermin裂解及膜完整性，区分细胞数量与单位活细胞功能')
detail=detail.replace('焦亡检测只在裂解性死亡及相应执行分子支持时追加，不能因项目主题而将焦亡设为唯一终点。','首轮平行保留gasdermin裂解与膜完整性等独立焦亡观察；出现执行层与死亡表型共同支持后，再追加通路阻断。凋亡读数不替代焦亡。')
detail=detail.replace('线粒体损伤和细胞死亡作为伴随终点，不替代Treg状态主问题。','首轮同时观察膜电位、存活及gasdermin裂解/膜完整性，判断状态效应是否伴随线粒体或死亡变化；主要终点仍为预设Treg状态。')
main="""page();p('实验验证｜在巨噬—T细胞体系中连接表型与损伤','h')
p('首轮回答三个相连的问题：矽尘相关巨噬刺激是否改变目标T细胞表型；候选干预是否改变这一响应；是否伴随线粒体及死亡终点变化。IL1R2优先回答Treg功能环节，TTC39C更接近损伤问题，SKIL侧重状态调控。')
tab([['组别','巨噬刺激条件','接收端T细胞'],['A','对照巨噬','非靶向对照'],['B','矽尘处理巨噬','非靶向对照'],['C','对照巨噬','候选干预'],['D','矽尘处理巨噬','候选干预']],[50,210,245])
p('主要比较为交互效应Δ=(D−C)−(B−A)，检验候选干预是否改变巨噬刺激效应。优先小鼠体系；采用肺来源巨噬细胞与目标Treg或活化CD4细胞，先导实验确定刺激和递送参数。')
tab([['候选','预设主要终点','共同观察与下一步'],['IL1R2 / Treg','等量活Treg的抑制功能','同步记录膜电位、存活及独立焦亡读数；区分IL-1诱饵结合与胞内作用'],['TTC39C / 活化CD4','存活/凋亡反应及执行层','测线粒体膜电位；焦亡使用gasdermin裂解及膜完整性等独立观察'],['SKIL / Treg','成熟Treg的FOXP3维持率','观察状态变化是否伴随线粒体或死亡表型；并行检测SnoN及SKI']],[100,170,235])
p('共同控制：颗粒残留、递送及基础毒性；目标T细胞内读取终点；独立动物/供体作为重复。先导阶段确认焦亡执行分子可测性与检测对照，阴性结果保留。线粒体损伤与焦亡分别判定，时间与干预实验再检验先后关系。')
p('推进顺序：可重复表型和交互效应→独立干预及回补→上游配体定位与阻断/救援。没有明确上游配体时先验证候选功能。详细实施、参数、统计和暂停条件见附录B。','small')
"""
s=s[:a]+main+s[b:]
# New bridge page.
bridge="""page();p('三候选证据连接｜从巨噬来源到T细胞响应','h')
p('三者对应同一研究主线的不同环节。下表把表达覆盖、网络预测、程序关联和首轮实验问题分开列明；其中AllMac为全体巨噬来源，AMlike为其相关亚群范围，二者不作独立证据相加。')
tab([['环节','IL1R2 / Treg','TTC39C / ActCD4','SKIL / Treg'],['巨噬来源信号','AllMac的IL1A/IL1B与Treg IL1R2在两组满足覆盖','当前联合差异网络未纳入该候选；上游未定位','AllMac的CXCL10、CCL2、SPP1与相应Treg受体在两组满足覆盖'],['受体及网络关系','IL1R2为诱饵受体；结合不等于胞内传导','首轮以巨噬刺激体系评价，不指定配体','CXCR3、CCR2、CD44或ITGAV+ITGB1；SKIL靶排名192、235、155'],['候选/表型依据','疾病上调，Inhib关联及Treg功能文献','ApoExec联合rho=0.943，P=0.00347，FDR=0.0556','TAct联合rho=0.965，P=0.00208，FDR=0.0744'],['损伤依据','PyroExec P=0.988；当前无直接阳性焦亡关联','TLike中MitoISR时间调整rho=0.733，P=0.0185，FDR=0.567','ApoExec疾病调整P=0.917；当前不支持该损伤关联'],['首轮检验','巨噬刺激下的Treg抑制功能是否受IL1R2影响，并伴随线粒体/死亡变化？','候选是否改变巨噬刺激后的存活/凋亡，并伴随膜电位和独立焦亡读数变化？','候选是否改变巨噬刺激后的Treg状态，是否伴随损伤终点变化？']],[75,143,143,144])
p('SKIL的配体-靶排名来自联合Treg上调目标集先验，表达覆盖以主队列发送端和接收端核对；它不是来源特异的处理效应或已验证因果轴。AMlike的配体表达不能替代同一来源的整条实验验证。')
p('H2-DM相关数据库记录另列抗原呈递网络信息，不作为普通表面配体-CD4作用轴推荐。TTC39C的网络NA不否定其凋亡程序关联。完整连接、组别和数据库来源见补充材料。','small')
"""
pos=s.index("page();p('候选比较与外部方向'");s=s[:pos]+bridge+s[pos:]
# Appendices and consolidated index before final build.
pos=s.index('# Layout uses')
append="\npage();p('附录A｜筛选阈值、集合及候选去向','h')\np('候选筛选策略见正文示意图C；本附录保留完整阈值、集合数量和分层去向。55＋42－7＝90为并集扩展，重点候选从完整集合结合功能证据和实验价值选择。')\nscreening_methods()\n"
append+="for m in manifest:\n if m['id'] not in ['F21','F22']:continue\n page();p(m['module']+'｜'+m['title'],'h2');p(results[m['id']]);im=PILImage.open(m['file']);k=min(505/im.width,540/im.height);story.append(Image(m['file'],width=im.width*k,height=im.height*k));p('图'+str(int(m['id'][1:]))+'。'+m['caption'],'small')\n"
append+=detail
append+="""page();p('附录C｜图表与补充材料索引','h')
p('配套补充材料按图号提供图值表、完整统计及来源说明。图号沿用各分析模块的固定标识；图21—22为筛选细节，位于附录A。')
tab([['图号','内容','配套路径']]+[['图'+str(int(m['id'][1:])),m['title'],'plot_data/'+m['id']+'.tsv'] for m in manifest[:12]],[45,320,140])
page();p('附录C｜图表索引（续）','h')
tab([['图号','内容','配套路径']]+[['图'+str(int(m['id'][1:])),m['title'],'plot_data/'+m['id']+'.tsv'] for m in manifest[12:]],[45,320,140])
p('统计表、图值来源及方法','h2');p('补充材料README.txt列明各目录用途；file_inventory.json提供文件哈希。tables包含完整关联、程序效应、差异、GO/KEGG/GSEA及候选检测；evidence包含背景指标定义、H2-DM来源核对及三候选连接来源。原始P、FDR、模型与检验家族分别保留。','small')
"""
s=s[:pos]+append+s[pos:]
# Remove unnecessary generalized caveats and version handling words, while retaining key limitations in discussion.
s=s.replace('PROS1相关功能支持已撤回。','').replace('PROS1路线未纳入。','')
s=s.replace('；不作为统一二元交集','').replace('，不作为统一二元交集','')
s=s.replace('全部DE及原完整家族FDR沿用已复核结果。','')
s=s.replace('旧版六样本平均比例权重','六样本平均比例权重').replace('历史轨迹','状态轨迹').replace('冻结','').replace('v18','')
s=s.replace('主队列为3+3个样本','分析材料为小鼠；主队列为3+3个样本')
s=s.replace('项目支持候选选择，不预先决定湿实验方向','当前证据支持预实验，机制连接由干预检验')
s=s.replace('图谱、细胞通讯、轨迹与功能评分均为计算层证据','图谱、细胞通讯、轨迹与功能评分均为计算层证据')
# Avoid data index promise without concrete package.
a=s.index("p('补充材料','h2')") if "p('补充材料','h2')" in s else -1
if a>=0:
 b=s.index('# Layout uses',a) if '# Layout uses' in s[a:] else len(s)
 # Do not erase appended chapters: only replace original two paragraphs until appendices.
 end=s.index("page();p('附录A",a)
 s=s[:a]+"p('补充材料','h2');p('完整统计、图值及来源索引随本报告提供，见附录C及配套材料。','small')\n"+s[end:]
# Add source reference for corrected antigen-presentation annotation before reference loop.
s=s.replace('for n,title,url in refs:',"refs.append(('[10]','Reactome. H2-DMb1, lysosomal membrane; antigenic peptide loading. R-MMU-2130673.','https://www.reactome.org/content/detail/R-MMU-2130673'))\nfor n,title,url in refs:")
(O/'report.py').write_text(s,encoding='utf-8')
# Record exact statistical evidence and provenance without new inferential tests.
a=pd.read_csv(H/'lung_append_analysis_20261017_v18/tables/all_associations.tsv',sep='\t');z=a[(a.gene=='SKIL')&(a.program=='TAct')&(a.scope=='Treg')&(a.cohort=='Joint')&(a.model=='JointStratumCondition')&(a.variant=='Raw')];assert len(z)==1;z.to_csv(O/'evidence/SKIL_adjusted_source.tsv',sep='\t',index=False)
lr=pd.read_csv(H/'lung_three_cohort_20261017_v16/data/tables/expressed_curated_LR.tsv',sep='\t');lr[(lr.scope=='Treg')&lr.ligand_parts.isin(['CXCL10','CCL2','SPP1','H2-DMB1','H2-DMB2'])].to_csv(O/'evidence/SKIL_route_coverage.tsv',sep='\t',index=False)
db=pd.read_csv(H/'lung_three_cohort_20261017_v16/data/resources/CellChat_interaction.tsv',sep='\t');db[db.astype(str).apply(lambda x:x.str.contains('H2-DMB1_CD4|H2-DMB2_CD4')).any(axis=1)].to_csv(O/'evidence/H2_DM_database_records.tsv',sep='\t',index=False)
info={'resource':'CellChatDB.mouse.rda','sha256':'dd59b19646e546b3ed41f394ad2cf64cf44f3d3fcfbb58ac176f0fe1322f9413','record_version':'CellChatDB v1','record_evidence':'KEGG: mmu04514','record_definition':'H2-DMb1 or H2-DMb2 to Cd4, single-gene entries; not complete H2-DM alpha/beta complex','mapping':'mouse symbols canonicalized to upper case; no human-to-mouse remapping for these CellChat edges','interpretation':'antigen-presentation-related network annotation; not recommended as ordinary cell-surface ligand-CD4 signaling','external_reference':'https://www.reactome.org/content/detail/R-MMU-2130673','boundary':'database records cite pathway annotation, not direct functional proof of these proposed signaling edges; NicheNet rankings retained as model priors'}
(O/'evidence/H2_DM_interpretation.json').write_text(json.dumps(info,ensure_ascii=False,indent=2),encoding='utf-8')
shutil.copy2(H/'lung_client_rebuild_20261017_v20/review_bundle/baseline_v19/code/v16/export.py',O/'evidence/background_definition_source.py')
print('PREPARED',O)