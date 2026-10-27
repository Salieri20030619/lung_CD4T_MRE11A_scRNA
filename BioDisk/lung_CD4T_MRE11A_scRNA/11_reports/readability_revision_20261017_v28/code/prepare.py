from pathlib import Path
import shutil,json,pandas as pd
H=Path(r'C:\Users\31937');O=H/'lung_readability_revision_20261017_v28';O.mkdir(exist_ok=True);B=H/'lung_mainline_revision_20261017_v27'
for d in ['figures','plot_data','evidence']:(O/d).mkdir(exist_ok=True)
for n in ['figure_manifest.json','editorial_replacements.json','extra_manifest.json']:shutil.copy2(B/n,O/n)
shutil.copy2(B/'figures/S_C.png',O/'figures/S_C.png')
s=(B/'report.py').read_text(encoding='utf-8').replace('矽肺_CD4候选研究_主线修订_v27_20261017','矽肺_CD4候选研究_可读性修订_v28_20261017').replace('矽肺免疫全景_CD4候选与验证_报告_v27.pdf','矽肺免疫全景_CD4候选与验证_报告_v28.pdf')
# Reader-facing tables use larger text. Technical captions remain compact.
s=s.replace("REPLACEMENTS=json.loads", "styles['table']=ParagraphStyle('table',fontName='CN',fontSize=10,leading=15,spaceAfter=3,wordWrap='CJK')\nREPLACEMENTS=json.loads")
s=s.replace("html.escape(clean(str(x))),styles['small']","html.escape(clean(str(x))),styles['table']")
# Page references are resolved by two identical layout passes.
pos=s.index('story=[]')
s=s[:pos]+"previous_index=json.loads((O/'report_page_index.json').read_text(encoding='utf-8')) if (O/'report_page_index.json').exists() else []\ndef loc(title):\n return str(next((r['page'] for r in previous_index if r['heading']==title),'待排'))\ndef figloc(m):return loc(m['module']+'｜'+m['title'])\n\n"+s[pos:]
m=json.loads((O/'figure_manifest.json').read_text(encoding='utf-8'))
for x in m:
 x['number']={'F23':21,'F24':22,'F21':23,'F22':24}.get(x['id'],int(x['id'][1:]))
 if x['id'] in ['F23','F24']:x['module']='候选表达与关联证据'
 if x['id']=='F19':
  x['file']=str(O/'figures/F19.png');x['title']='校正AUPR排序提示各CD4状态的优先响应配体';x['caption']='每状态按校正AUPR降序取前6个配体，展示三状态并集；色深同样表示校正AUPR，格内标各状态原排名。空白为该状态无对应活动结果，不填为0。各状态目标集与背景不同，数值不是实测通讯强度。完整Pearson指标保留于补充表。'
 if x['id']=='F24':
  x['file']=str(O/'figures/F24.png');x['caption']='坐标为未作协变量残差化的原始样本表达与评分；圆形为GSE314943，三角为GSE183682。模型M1调整队列/时间层，M2进一步调整疾病状态；各面板并列两模型rho、P和FDR。统计量不是原始散点的简单相关。'
 x['caption']=x['caption'].replace('图23','图21').replace('图24','图22')
(O/'figure_manifest.json').write_text(json.dumps(m,ensure_ascii=False,indent=2),encoding='utf-8')
s=s.replace("str(int(m['id'][1:]))","str(m['number'])")
s=s.replace('图21-22','图23-24').replace('图23-24候选与验证','图21-22候选表达与关联证据').replace('图21—22为筛选细节','图23—24为筛选细节')
# Short opening definitions and mitochondrial results.
s=s.replace('ApoExec rho=0.943','凋亡程序(ApoExec) rho=0.943').replace('TAct rho=0.965','活化程序(TAct) rho=0.965')
s=s.replace('主队列ActCD4、Treg焦亡相关评分P为0.30、0.20，当前未检出显著组间差异。','主队列三类CD4状态的线粒体应激程序均向上（P=0.80、0.80、0.40），ActCD4、Treg焦亡相关评分P为0.30、0.20，均未检出显著组间差异。')
# Navigation early; not another full statistical comparison.
a=s.index("page();p('研究设计与数据范围'")
nav="""page();p('研究证据导航｜四个问题与三个候选切入点','h')
p('本报告依次回答T细胞发生何种变化、有哪些巨噬来源信号、损伤程序是否改变，以及哪些候选值得预实验。各模块提供不同层面的证据，连续排列不表示机制链已得到证明。')
tab([['研究问题','当前答案','阅读位置'],['T细胞有何变化？','Naive比例下降，ActCD4/Treg上升；组成检验未达名义显著。Treg候选存在疾病差异及功能程序变化。','图谱与状态见第'+loc('免疫全景｜GSE314943主队列的全肺细胞图谱')+'页起'],['巨噬信号如何联系T细胞？','配体-受体表达覆盖与靶先验提供候选连接，IL1R2诱饵结合与SKIL响应先验分别解释。','通讯见第'+loc('细胞通讯预测｜巨噬细胞至CD4状态的表达相容通讯候选')+'页；NicheNet见第'+loc('NicheNet｜校正AUPR排序提示各CD4状态的优先响应配体')+'页'],['线粒体/死亡程序改变了吗？','主队列线粒体应激、凋亡与焦亡程序未检出显著组间差异；候选关联提供另一层探索线索。','评分见第'+loc('焦亡/凋亡评分｜焦亡执行与炎症程序分开评价')+'页起'],['为何建议这三个候选？','IL1R2关注Treg与IL-1环境，TTC39C关注损伤反应，SKIL关注状态调控。来源与资格并不相同。','入选解释见第'+loc('三候选为何推荐｜来源、分层与研究价值')+'页']],[105,290,110])
p('如何使用后续实验建议','h2');p('先确认巨噬—T细胞模型及候选可测性，再选一个最可实施的候选进入功能干预。正文说明优先问题，附录B提供四组设计和可选路径。详见第'+loc('实验验证建议｜连接巨噬环境、T细胞表型与损伤')+'页及第'+loc('附录B｜建议研究问题与可选实验路径')+'页起。')
p('轨迹分析承担的任务','h2');p('状态排序检查Naive/ActCD4连续变化及死亡相关程序沿状态的分布。它提供活化CD4背景，不包含Treg，也不作为IL1R2或SKIL的Treg轨迹验证。')
"""
s=s[:a]+nav+s[a:]
# Move detailed formulas/parameters to a methods appendix.
a=s.index("page();p('分析方法','h')");b=s.index('# Results lead',a);methods=s[a:b];methods=methods.replace("p('分析方法','h')","p('附录D｜分析方法与参数','h')").replace("p('评分、关联与状态轨迹方法','h')","p('附录D｜评分、关联与状态排序方法','h')")
methods=methods.replace('活动Pearson及网络权重作为模型排序指标，保留各状态实际目标集和背景。','配体按校正AUPR降序排序，同值按配体名称；图19按同一指标着色并标排名。Pearson为另一个预测指标，在完整表中单列。配体-靶基因权重用于靶排名；保留各状态实际目标集和背景。')
brief="""page();p('读图所需的统计与术语','h')
p('统计以独立生物学样本为单位。差异表达使用样本汇总计数；联合模型控制队列/时间层。候选筛选采用目标效应原始P<0.05，并同时报告完整家族FDR。P<0.05且FDR未达阈值的结果称为名义显著。')
tab([['名称','本报告中的含义'],['TAct','T细胞活化相关转录程序；包含CD69、IL2RA、ICOS等成员'],['Inhib','抑制性受体/免疫调节相关转录程序；包含PDCD1、CTLA4、TIGIT等成员，非实测抑制功能'],['ApoExec','凋亡执行相关基因的转录程序'],['PyroExec','焦亡执行相关基因的转录程序'],['MitoISR','线粒体应激相关转录程序；不等同膜电位或呼吸实测'],['M1 / M2','M1：调整队列/时间层；M2：在M1基础上进一步调整疾病状态'],['校正AUPR','NicheNet配体预测排序指标；用于图19的选取、排名及着色']],[115,390])
p('评分图呈现样本级转录程序，相关图展示原始点并明确统计模型。完整基因集、评分公式、排列数量、背景敏感性和各类FDR家族见方法附录及配套表。')
"""
s=s[:a]+brief+s[b:]
# Current marker-table evidence, not newly invented annotation rules.
pos=s.index('def early_framework():')
statefun="""def state_guide():
 page();p('CD4状态定义｜读懂图谱与候选所在细胞','h')
 p('Core CD4为Naive、ActCD4和Treg三类状态的分析合并范围。下表沿用项目注释，并列出已有标志表中帮助理解状态的表达特征；这些是组合性证据，不是新增单基因分类阈值。')
 tab([['状态','中文含义与范围','已有标志表达特征'],['Naive','初始样CD4状态，主队列933个细胞','CCR7、TCF7较丰富；六样本平均阳性率分别为66.9%、80.0%'],['ActCD4','活化样CD4状态，主队列801个细胞','CTLA4阳性率52.7%，高于Naive的14.2%；CCR7较低（22.3%），FOXP3较少（3.1%）'],['Treg','调节性T细胞状态，主队列577个细胞','FOXP3、IL2RA、CTLA4阳性率分别为57.3%、40.8%、79.1%'],['Core CD4','上述三状态的合并分析，共2,311个细胞','三状态均有CD3D/CD3E/TRAC表达；组成和候选分析按状态区分'],['TLike','时间队列GSE274372的T样细胞范围，494个细胞','用于时间背景，不能直接等同纯CD4、Treg或ActCD4']],[85,200,220])
 p('阳性率为样本内比例再跨6个独立样本平均，仅用于描述标志分布。状态名称沿用细胞注释；转录标签不直接指定湿实验分选门。完整标志表随补充材料提供。','small')

"""
s=s[:pos]+statefun+s[pos:]
s=s.replace("if m['id']=='F05':early_framework()","if m['id']=='F03':state_guide()\n if m['id']=='F05':early_framework()")
# Result lead changes and transitions, all from existing results.
pos=s.index('for m in manifest:')
upd={
'F11':'在CD4组成和分子变化基础上，本节检查巨噬细胞能提供哪些表达相容信号。巨噬与Naive、ActCD4、Treg均存在满足覆盖的配体-受体候选；弦图按来源和状态展示其分布。',
'F13':'转录因子网络补充说明T细胞的响应程序。Treg中NFKB1样本平均AUCell由0.0600升至0.0779，FOXP1由0.1539降至0.1354；这里报告描述性活性差异。',
'F16':'主队列线粒体应激程序MitoISR在Naive、ActCD4、Treg中的差值分别为+0.019、+0.012、+0.030；P分别为0.80、0.80、0.40，FDR为0.945、0.945、0.578，均未达名义显著。ApoExec在ActCD4和Treg的差值为+0.0404与+0.0352，P=0.10与0.20，FDR=0.274与0.433。这里检验疾病组间程序差异；后文候选关联回答另一问题。',
'F17':'本节检查Naive/ActCD4的状态连续性，为活化CD4的损伤研究提供背景，不验证Treg候选路线。四套DPT排序相关为0.904-0.999；描述性邻域主干帮助阅读状态方向。',
'F24':'候选表达与程序的关系在两种模型下分别评价：M1调整队列/时间层，M2进一步调整疾病状态。图中并列效应、P和FDR；原始点用于展示样本分布。疾病调整可检验共同组别差异对关联的解释程度。'}
s=s[:pos]+'results.update('+repr(upd)+')\n'+s[pos:]
# Dedicated candidate origin table; sources verified before authoring.
a=s.index("page();p('三候选证据连接｜")
selection="""page();p('三候选为何推荐｜来源、分层与研究价值','h')
p('先以差异、功能和巨噬来源信号构建55基因交集，再与42个补充候选合并，重叠7个，共90个评价基因。以下三者从完整评价集合中按研究价值推荐，并非都由同一三项交集直接筛出。')
tab([['候选/状态','进入来源与计算层级','推荐依据及待检验问题'],['IL1R2 / Treg','55基因交集与42个补充候选均包含；属于52个计算支持基因','疾病上调证据及Treg功能研究较明确，优先检验Treg状态和IL-1环境响应。当前无直接阳性焦亡关联（P=0.988）。'],['TTC39C / ActCD4','来自42个补充候选，不在55基因交集；属于37个证据待补基因，功能维度待补','ApoExec关联及存活/凋亡功能文献提供损伤方向价值。当前联合DE P=0.163、FDR=0.803，未纳入该联合差异网络，先检查原代CD4可测性及实际干预方向。'],['SKIL / Treg','来自55基因交集，不在42个补充候选；属于52个计算支持基因','疾病差异、TAct关联及调控背景支持状态路线。先确认SnoN蛋白及状态响应；疾病调整ApoExec P=0.917，不支持将其直接作为损伤调控靶点。']],[90,170,245])
p('计算资格与实验价值的关系','h2');p('资格分层记录既定计算条件是否得到支持；文献、候选关联与可操作性决定值得检验的问题。TTC39C的文献先例支持探索性预实验，不能据此把其计算分层改写为已通过。首轮优先核查IL1R2，并依据模型和资源选择一个最可实施候选进入功能干预。')
p('完整成员来源和阶段去向见补充筛选表；三个候选的具体信号连接、调整模型和外部方向见随后两页。','small')
"""
s=s[:a]+selection+s[a:]
s=s.replace('建议验证路径｜连接巨噬环境、T细胞表型与损伤','实验验证建议｜连接巨噬环境、T细胞表型与损伤')
s=s.replace('ApoExec联合rho=0.943','ApoExec M1 rho=0.943').replace('TAct联合rho=0.965','TAct M1 rho=0.965').replace('疾病调整rho=0.922','M2 rho=0.922').replace('疾病调整rho=0.897','M2 rho=0.897')
s=s.replace("p('当前证据支持预实验，机制连接由干预检验','h2');p('疾病差异、程序关联、网络先验和干预功能是不同证据层级。IL1R2具备较清楚的疾病差异和Treg功能研究基础；TTC39C提供损伤相关程序切入点；SKIL具备Treg状态和调控背景。三者均有预实验价值，但不对应已量化的成功概率。')","p('当前证据支持预实验，机制连接由干预检验','h2');p('疾病差异、程序关联与网络先验提供可检验假说；是否存在巨噬刺激特异的候选作用，以及是否连接线粒体损伤或焦亡，需由接收端干预与独立终点检验。预实验推荐不代表已量化的成功概率。')")
pros='PROS1-AXL的胶质瘤文献背景不能直接用于小鼠CD4。TYRO3/MERTK的转录表达未满足所设覆盖条件；RNA低检出不是蛋白不存在的证明。'
s=s.replace(pros,'')
# Insert methods appendix immediately before index; preserve exact technical detail.
a=s.index("page();p('附录C｜")
methods+="p('备选路线的受体资格','h2');p('PROS1路线在候选评价中作为备选网络检查。PROS1-AXL的胶质瘤文献背景不能直接用于小鼠CD4；TYRO3/MERTK转录表达未满足所设覆盖条件，因此不作为本项目CD4功能性受体支持。RNA低检出不证明蛋白不存在。','small')\n"
s=s[:a]+methods+s[a:]
s=s.replace("['图号','内容','配套路径']","['图号','内容','页码','配套路径']").replace("m['title'],'plot_data/'+m['id']+'.tsv'","m['title'],figloc(m),'plot_data/'+m['id']+'.tsv'").replace('[45,320,140]','[40,280,40,145]')
# Index in display order, not internal IDs.
s=s.replace("for m in manifest[:12]","for m in sorted(manifest,key=lambda x:x['number'])[:12]").replace("for m in manifest[12:]","for m in sorted(manifest,key=lambda x:x['number'])[12:]")
s=s.replace('图号沿用各分析模块的固定标识；图23—24为筛选细节，位于附录A。','图号按阅读顺序连续排列；配套文件使用稳定分析编号，表中逐项对应。图23—24为附录A筛选细节。')
s=s.replace('初筛网络分支采用校正AUPR>0的合格配体，按活动排序取前30个','初筛网络分支采用校正AUPR>0的合格配体，按校正AUPR降序取前30个')
(O/'report.py').write_text(s,encoding='utf-8')
z=pd.read_csv(H/'lung_screening_display_20261017_v24/tables/candidate90_stage_membership.tsv',sep='\t');z[z.gene.isin(['IL1R2','TTC39C','SKIL'])].to_csv(O/'evidence/candidate_origin_and_status.tsv',sep='\t',index=False)
a=pd.read_csv(H/'lung_client_rebuild_20261017_v20/server_tables/cd4_markers.tsv',sep='\t');a.groupby(['state','gene']).fraction.mean().reset_index().to_csv(O/'evidence/state_marker_sample_mean.tsv',sep='\t',index=False)
print('PREPARED')