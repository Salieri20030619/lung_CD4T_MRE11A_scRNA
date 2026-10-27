from pathlib import Path
import pandas as pd,numpy as np,json,math,html,hashlib,zipfile,shutil
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak,Image,KeepTogether
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.pagesizes import A4
O=Path(r'C:\Users\31937\lung_three_cohort_20261017_v16');T=O/'data/tables';F=O/'final_tables';V=O/'figures';V.mkdir(exist_ok=True)
D=Path(r'C:\Users\31937\Downloads\线粒体焦亡_三队列四轴重筛_v16_20261017');D.mkdir(exist_ok=True)
rd=lambda n:pd.read_csv(T/n,sep='\t',low_memory=False)
de=rd('DE_all_models.tsv');s=rd('three_axis_candidates_all_contrasts.tsv');q=rd('candidate_sample_expression_QC.tsv');inv=rd('sample_inventory.tsv');routes=rd('candidate_complete_routes.tsv');whole=rd('whole_lung_source_expression.tsv.gz')
summary=pd.read_csv(F/'四轴候选总表_55基因.tsv',sep='\t');pr=pd.read_csv(F/'优先验证候选_5基因.tsv',sep='\t');hist=pd.read_csv(F/'历史候选逐项回溯.tsv',sep='\t');meta=json.loads((F/'summary.json').read_text(encoding='utf-8'))
# Plots use static scientific coordinates; each dot denotes an independent sample.
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
fig,ax=plt.subplots(figsize=(8,3.0));labels=['Nominal DE','DE + function','+ macrophage route','Literature-prioritized'];vals=[2765,111,55,5];ax.barh(np.arange(4),vals,color=['#66849c','#338c92','#17756f','#df9c35']);ax.set_yticks(np.arange(4),labels);ax.invert_yaxis();ax.set_xscale('log');ax.set_xlim(1,6000);ax.set_xlabel('Unique genes (log scale)');[ax.text(v*1.10,i,str(v),va='center') for i,v in enumerate(vals)];fig.tight_layout();fig.savefig(V/'funnel.png');plt.close(fig)
genes=['IL1R2','DDIT4','SKIL','FOXP1','GSN','TNFRSF9','DUSP5','STAT1','PMAIP1','GBP5','ERN1','PPM1K','CD44','TNFRSF1B'];scopes=['Core','Naive','ActCD4','Treg'];z=de[de.model=='joint_stratified'];fc=z.pivot(index='gene',columns='scope',values='logFC').reindex(index=genes,columns=scopes);pv=z.pivot(index='gene',columns='scope',values='PValue').reindex(index=genes,columns=scopes);fd=z.pivot(index='gene',columns='scope',values='FDR').reindex(index=genes,columns=scopes)
fig,ax=plt.subplots(figsize=(6.8,5.4));im=ax.imshow(np.ma.masked_invalid(fc),cmap='RdBu_r',vmin=-4,vmax=4,aspect='auto');ax.set_xticks(range(4),scopes);ax.set_yticks(range(len(genes)),genes)
for i in range(len(genes)):
 for j in range(4):
  v=fc.iloc[i,j];txt='NA' if pd.isna(v) else f'{v:+.2f}'+('**' if fd.iloc[i,j]<.05 else '*' if pv.iloc[i,j]<.05 else '')
  ax.text(j,i,txt,ha='center',va='center',fontsize=9,color='white' if pd.notna(v) and abs(v)>2.8 else 'black')
fig.colorbar(im,ax=ax,label='Joint model log2FC');ax.set_title('* raw P < 0.05; ** FDR < 0.05; NA = not tested');fig.tight_layout();fig.savefig(V/'joint_heatmap.png');plt.close(fig)
fig,axs=plt.subplots(2,3,figsize=(10,5.4));pairs=list(zip(pr.gene,pr.scope))
for ax,(gene,scope) in zip(axs.flat,pairs):
 zz=q[(q.gene==gene)&(q.scope==scope)&q.dataset.eq('GSE314943')]
 for gi,gr in enumerate(['Saline','Silica']):
  zzg=zz[zz.group==gr].sort_values('sample');xx=gi+np.array([-.1,0,.1]);ax.scatter(xx,zzg.CP10k,color=['#357ca5','#c75a3b'][gi],s=35,label='raw' if gi==0 else None);ax.scatter(xx,zzg.BG1_CP10k,facecolors='none',edgecolors=['#357ca5','#c75a3b'][gi],s=70,label='BG1 subset' if gi==0 else None)
 ax.set_xticks([0,1],['Saline','Silica']);ax.set_title(gene+' / '+scope);ax.set_ylabel('Pseudobulk CP10k')
axs.flat[0].legend(fontsize=7);axs.flat[5].axis('off');fig.tight_layout();fig.savefig(V/'expression.png');plt.close(fig)
fig,ax=plt.subplots(figsize=(7.5,3.1));zz=de[(de.model=='274_time')&de.gene.eq('FOXP1')&de.contrast.ne('time_omnibus')];ax.bar(np.arange(3),zz.logFC,color='#327d91');ax.axhline(0,color='black',linewidth=.6);ax.set_xticks(range(3),['D7 / D0','D28 / D0','D56 / D0']);ax.set_ylabel('FOXP1 log2FC (TLike)');ax.set_ylim(-2.2,.25)
for i,r in enumerate(zz.itertuples()):ax.text(i,-.08,f'P={r.PValue:.3g}\nFDR={r.FDR:.3g}',ha='center',va='top',fontsize=9)
fig.tight_layout();fig.savefig(V/'time_foxp1.png');plt.close(fig)
pdfmetrics.registerFont(TTFont('CN',r'C:\Windows\Fonts\ARIALUNI.TTF'))
styles={
 'title':ParagraphStyle('title',fontName='CN',fontSize=22,leading=31,textColor=colors.HexColor('#163c50'),spaceAfter=16),
 'h1':ParagraphStyle('h1',fontName='CN',fontSize=16,leading=23,textColor=colors.HexColor('#156b71'),spaceAfter=12),
 'h2':ParagraphStyle('h2',fontName='CN',fontSize=11.5,leading=17,spaceBefore=10,spaceAfter=5,textColor=colors.HexColor('#163c50')),
 'body':ParagraphStyle('body',fontName='CN',fontSize=9.8,leading=16,spaceAfter=8,wordWrap='CJK'),
 'small':ParagraphStyle('small',fontName='CN',fontSize=8,leading=12,spaceAfter=6,wordWrap='CJK'),
 'cell':ParagraphStyle('cell',fontName='CN',fontSize=8.2,leading=12,wordWrap='CJK'),
 'head':ParagraphStyle('head',fontName='CN',fontSize=8.4,leading=12,textColor=colors.white,wordWrap='CJK')}
story=[]
def p(t,sty='body'):return Paragraph(html.escape(str(t)).replace('\n','<br/>'),styles[sty])
def add(t,sty='body'):story.append(p(t,sty))
def heading(t):add(t,'h1')
def page(t):story.append(PageBreak());heading(t)
def table(headers,rows,widths):
 data=[[p(v,'head') for v in headers]]+[[p(v,'cell') for v in row] for row in rows];tt=Table(data,colWidths=widths,repeatRows=1,hAlign='LEFT');tt.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#17566a')),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),6),('RIGHTPADDING',(0,0),(-1,-1),6),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.HexColor('#f1f6f7'),colors.white]),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#17566a'))]));story.append(tt);story.append(Spacer(1,8))
def img(name,w=490,h=None):
 from PIL import Image as PI
 a=PI.open(V/name);h=h or w*a.height/a.width;story.append(Image(str(V/name),width=w,height=h));story.append(Spacer(1,6))
def num(v):return 'NA' if pd.isna(v) else f'{v:.3g}'
def stat(g,sc,model='joint_stratified',contrast='Silica_vs_Saline'):
 a=de[(de.gene==g)&(de.scope==sc)&(de.model==model)&(de.contrast==contrast)]
 return a.iloc[0] if len(a) else None
add('三队列四轴重筛\n线粒体凋亡 / 焦亡项目','title');add('v16 | 2026-10-17 | 中文内部复审版','h2')
add('结论：获得55个三项生物学交集基因；结合创新问题、细胞内表达及队列一致性，形成4个优先验证方向和1个备选。原v13-v15保留。')
add('四轴分别为：独立样本差异表达、内源性线粒体凋亡/焦亡功能、巨噬细胞作用入口、具体研究问题的创新性。功能轴取“凋亡或焦亡”的并集，不要求一个基因必须同时属于两种死亡机制。')
rows=[]
for r in pr.itertuples():rows.append([r.gene,r.scope,f'{r.logFC:+.2f}',num(r.PValue),num(r.FDR),r.review_priority])
table(['基因','重点接收范围','联合log2FC','原始P','FDR','本轮定位'],rows,[64,88,78,78,78,104])
add('IL1R2：优先核查Treg的抗焦亡/反馈保护功能。DDIT4：优先核查线粒体损伤与焦亡之间的调节节点。SKIL：优先核查Treg应激-p53凋亡调控。FOXP1：优先核查跨队列下降的CD4凋亡易感性；其巨噬来源证据较弱。GSN作为备选，关注胞内线粒体保护。')
add('“5个优先候选”是可追溯的文献和证据分层决定，不是把创新性伪装成P值阈值，也不表示其余50个交集基因被统计学否定。此版本完成候选重筛，未将通讯预测升级为已证实的因果机制。','small')
add('复核入口：附件中“优先验证候选_5基因”“四轴候选总表_55基因”“全部155项候选对比”及“历史候选逐项回溯”。','small')
page('1  新定义恢复了被旧条件过度排除的候选')
img('funnel.png')
add('三队列预设主分析共有2,765个基因在至少一个指定对比中达到原始P<0.05；其中111个有功能成员依据，55个同时具有匹配接收状态的巨噬细胞作用入口。55个均有疾病分析来源，其中4个还进入独立时间分支：CD44、FOXP1、HIF1A、TNFRSF1B；本轮没有仅由时间分支新增的交集基因。')
table(['层面','实际结果','计数解释'],[['功能参考全集','495基因','GO/Reactome 489基因，加6个原始文献补充；并非全部差异'],['通路资源内三项交集','49基因','正式注释足以进入功能轴，不需模块相关性再过关'],['文献补充三项交集','6基因','DUSP5、FOXP1、GSN、IL1R2、STAT1、TNFRSF9'],['接收状态分布','Core 41；Naive 11；ActCD4 17；Treg 26；TLike 4','状态间有重叠，不能相加当独立基因或独立队列'],['完整明细','99个基因-状态组合；155项对比','同一基因可在原队列、联合和时间模型重复出现']],[110,165,215])
add('本次不再把“自定义损伤评分相关+组别校正P<0.05+留一法P<0.05+背景过滤+NicheNet前200”串成所有候选的必过链。直接受体使用配体-受体入口；下游调节基因使用配体活性和靶预测入口。关联和敏感性结果用于分层，不能抹去正式功能证据或原始差异。')
page('2  三队列参与分析，但不混成一个22例疾病比较')
table(['队列','接收细胞/独立样本','进入本轮的方式'],[['GSE314943','2,311个CD4；3 Saline+3 Silica','Core、Naive、ActCD4、Treg各自样本级差异；冻结AllMac用于通讯'],['GSE183682','1,413个CD4；4只独立小鼠','D7和D56各1对照+1暴露；与前队列形成10样本分层共同效应模型'],['GSE274372','494个保守TLike；12只小鼠','D0/D7/D28/D56各3；单独时间模型，不能当作后期同期盐水对照设计']],[100,165,225])
add('GSE183682中的表达、BCR和TCR文件不重复计算为生物学样本。两个时间点各只有1+1，因此本报告只给各时间点的描述性方向，未生成独立时间点的差异P值。联合模型为~study_time_stratum+condition，设计满秩，10样本、4参数、残差自由度6。它估计共同条件效应，不等于两个独立队列各自显著。')
add('GSE274372的作者CD4/Treg对象为888个细胞；沿用既往R08=0/1的保守TLike定义后保留494个。D0_R2仅8个TLike细胞，主时间模型保留；另做至少10细胞的敏感性模型（11样本）。TLike不能直接改名为纯CD4或Treg。')
add('通讯发送端：GSE314943冻结AllMac 1,796细胞（AMlike 700细胞另存覆盖证据）；GSE274372 Broad=Mac共894细胞，出现在11个样本。每个通讯分析分组再按至少10细胞决定可评估样本。GSE183682无已验证的巨噬细胞注释，本轮用于接收端差异整合，没有补造该队列的巨噬通讯结果。')
add('样本级数据覆盖总计22只小鼠，但疾病联合推断仅使用前两队列10只；第三队列的时间结果单独陈述。Core和三个亚群来自同一批小鼠，不算四个独立验证。','small')
page('3  执行阈值和统计方法可完整追溯')
table(['轴/步骤','冻结口径'],[['差异表达','原始UMI按独立样本汇总；edgeR QL，filterByExpr，TMM，robust离散度。目标效应原始P<0.05入选；|log2FC|≥0.5仅作为效应分层。'],['多重检验','各模型×状态×预设对比单独BH；保留全部测试基因。17个差异检验家族分别报告。P和FDR均保留；名义显著不写作校正后显著。'],['功能轴','Mouse GOALL中内源性凋亡及其调控、凋亡性线粒体变化和焦亡相关条目；Reactome R-MMU-109606、R-MMU-5620971；原始实验文献补充另列。'],['表达与复合物','通讯样本细胞数≥10；一个基因在至少5%细胞检出，并在同组至少50%可评估样本满足。配体及受体必需亚基均满足；不要求每个样本均检出。'],['直接受体入口','CellChatDB.mouse蛋白类LR；完整亚基覆盖。另记录IL1A/B-IL1R2的文献支持诱饵结合，不把它当激活信号。'],['下游靶入口','NicheNet v2 mouse；按完整名义差异上/下调集合分别估计活性；候选配体须有合格巨噬LR；corrected AUPR>0且前30，靶在表达背景内前250。'],['创新性','所有55基因统一检索疾病、死亡机制、CD4死亡三个场景，并补常用别名；重点方向审阅原始研究和矛盾结果。未检出不等于全球首次。']],[100,390])
add('不同队列按基因符号和明确旧名映射，重复符号求和，记录原始feature映射。联合模型仅使用两队列共同可映射基因。3+3精确置换的最小双侧P=0.10限制不适用于本次edgeR模型；未为取得更小P而替换评分检验。')
add('NicheNet活性/靶排序是网络预测，不产生独立小鼠处理效应P值；这里没有运行CellChat概率或置换检验，使用的是其先验数据库。功能知识库和通讯先验可能共享文献，不能算相互独立的机制验证。','small')
page('4  联合差异支持多个状态中的候选变化')
img('joint_heatmap.png',w=480)
add('热图只展示讨论相关基因，完整全基因结果在附件。未显示星号不等于原队列无差异；入选采用原队列或联合或独立时间分支的预设并集。Naive联合结果整体FDR较高，保留原始P的名义证据定位。')
add('IL1R2、TNFRSF9、FOXP1在Core与Treg中均通过联合模型FDR<0.05；SKIL在Treg中通过。DDIT4和GSN的重点状态属于名义显著。TNFRSF9统计强，并不自动意味着创新性优先；STAT1在GSE183682多个状态/时间方向减弱或相反，也不能由联合名义差异恢复原先完整STAT1-GSDMD机制主张。','small')
page('5  IL1R2提供直接结合入口，问题是保护性反馈')
r=stat('IL1R2','Treg');add(f'Treg：联合log2FC={r.logFC:+.3f}，P={num(r.PValue)}，FDR={num(r.FDR)}。Core联合log2FC=+3.891，P=1.78e-7，FDR=2.02e-4。GSE183682的Treg在D7与D56描述性logCPM差值分别为+2.408、+6.906，同向。')
add('原六样本Treg中，三只暴露小鼠的IL1R2阳性比例为28.1%、30.1%、27.3%；对应阳性细胞50、47、44，其中≥2 UMI细胞29、24、22。所有这些阳性Treg均检出所用T细胞标记集合中的至少一个。BG1子集仍保留升高方向。')
table(['证据段','本次证据','正确解释'],[['发送端至受体','AllMac IL1A/B和Treg IL1R2满足同组覆盖；双向分组明细均保存','这是表达兼容的诱饵结合入口'],['功能依据','2025原始研究：IL1R2结合ENO1并限制巨噬细胞糖酵解和GSDMD相关焦亡','可作为跨细胞类型的功能假说，不是已证实Treg机制'],['来源比较','全肺Broad注释中，Silica组Mac平均IL1A CP10k约0.373；IL1B在中性粒细胞更高','优先检验巨噬IL1A环境；IL1B不能称巨噬特异'],['创新问题','矽肺Treg中的IL1R2是否限制ENO1相关焦亡/损伤','“保护性反应”与“驱动损伤”必须分别实验判定']],[100,205,185])
add('IL1R2缺少TIR信号域；不能画成“IL1结合IL1R2后直接激活NF-κB/焦亡”。接收端RNA升高与保护性功能相容，也可能是活化伴随变化。第一步应在纯化Treg确认蛋白和胞内定位，随后测试IL1R2干预是否改变ENO1活性及死亡终点。')
add('原始依据：https://pmc.ncbi.nlm.nih.gov/articles/PMC12533375/','small')
page('6  DDIT4与SKIL提供两类不同的损伤调节问题')
for gene,scope in [('DDIT4','ActCD4'),('SKIL','Treg')]:
 r=stat(gene,scope);add(f'{gene} / {scope}：log2FC={r.logFC:+.3f}；P={num(r.PValue)}；FDR={num(r.FDR)}。','h2')
add('DDIT4：GSE183682两个时间点ActCD4方向均为升高（+0.331、+1.250 logCPM）。原六样本暴露组ActCD4阳性比例约68.9%-82.4%，并有较多多UMI细胞。巨噬TNF-TNFRSF1B入口在Silica组覆盖，配体活性第3、DDIT4靶排名89；其他可评估路线见附件，不将所有路线都当成真实驱动。')
add('DDIT4已有非T细胞直接干预研究连接线粒体损伤、ROS-TXNIP-NLRP3及焦亡（PMID 38771134）；另一足细胞研究支持REDD1依赖的炎症/焦亡（PMID 39920111）。本项目的新问题是：巨噬环境能否在CD4中调用这一损伤节点。它不是“从未有人研究DDIT4焦亡”。')
add('SKIL：GSE183682两个时间点Treg方向均升高（+1.971、+1.608 logCPM）。Silica组三样本Treg阳性比例约52.8%-67.9%。GO鼠注释包含DNA损伤内源性凋亡的正调控（IMP）；SnoN原始研究支持其p53应激功能，但并非矽肺Treg研究。')
add('SKIL可连接巨噬CXCL10-CXCR3（活性第7，靶192）和CCL2-CCR2（活性第14，靶235）的表达兼容路径，两组均有覆盖。全肺来源比较中Silica组Mac的CCL2均值最高，但CXCL10不特异。NicheNet没有给出配体改变SKIL的因果方向；应把CCL2/CXCL10与SKIL扰动及p53/凋亡终点结合检验。')
add('SKIL的CCL2靶排名接近250边界；改用200时这条边会丢失，但CXCL10入口仍保留，基因层面不丢失。五个优先/备选候选均有至少一条≤200或直接受体入口。')
add('文献：https://pubmed.ncbi.nlm.nih.gov/38771134/；https://pubmed.ncbi.nlm.nih.gov/39920111/；https://pmc.ncbi.nlm.nih.gov/articles/PMC3897278/','small')
page('7  FOXP1具有跨队列下降，但巨噬上游仍需定位')
r=stat('FOXP1','Core');add(f'Core联合log2FC={r.logFC:+.3f}，P={num(r.PValue)}，FDR={num(r.FDR)}；Treg联合log2FC=-1.125，P=7.58e-4，FDR=0.0274。GSE183682 Core的D7/D56均下降；Treg D7约+0.040，接近零，D56下降，不能概括为所有状态均重复。')
img('time_foxp1.png')
add('第三队列TLike的D7、D28、D56均下降；D7通过该时间对比的FDR。去掉8细胞样本后D7 log2FC=-1.934，P=2.14e-5，FDR=0.00733，方向保留。第三队列是时间证据和TLike范围，不等同于纯CD4疾病验证。')
add('原六样本Core中FOXP1阳性比例44.5%-78.0%，多UMI细胞743个；BG1子集下降方向保留。当前可用的下游路线为SEMA3E-PLXND1和CD200-CD200R1，FOXP1靶排名约138和155。两种配体在全肺其他细胞中明显更高，因此不能据此认定巨噬细胞是关键来源。')
add('已有文献证明FOXP1可抑制BH3-only BIK等凋亡基因；Treg研究中，一项实验观察到FOXP1缺失伴凋亡增加，另一项在rTreg/aTreg比较中未见凋亡差异。该矛盾保留，说明必须在矽肺具体状态检验，而不能把“FOXP1维持Treg”当作新发现或普适损伤结论。')
add('依据：PMC3741507；PMC7534899（Fig.6）；PMC6534289（S3C）。已撤稿的FOXP1缺氧/复氧论文未用作依据。','small')
page('8  GSN保留为备选，旧主线按证据重新定位')
r=stat('GSN','ActCD4');add(f'GSN / ActCD4联合log2FC={r.logFC:+.3f}，P={num(r.PValue)}，FDR={num(r.FDR)}；两个GSE183682时间点均下降。原六样本中一个对照的GSN表达明显较高，BG1后该样本仍较高，但总体下降保留。')
add('可用PROS1-AXL等下游预测入口；PROS1并非巨噬特异。GSN的胞内线粒体保护有实验依据，但分泌型gelsolin的免疫作用不同，人和小鼠的VDAC调节也存在差异。本轮基因计数不能区分这些蛋白形式，因此列为备选而非直接建立机制。')
table(['对象','本轮决定','依据'],[['TNFRSF9','恢复交集，作为已知轴参照','Core/Treg联合差异强；但2016矽肺干预、2018肺泡巨噬机制及2024时间图谱已有研究。'],['DUSP5','保留Treg交集，不再唯一候选','联合Treg P=0.0265、FDR=0.1866；D7下降而D56轻度升高；直接焦亡文献为SCLC。'],['STAT1','恢复功能文献入口，降级为参照','Core/Treg/ActCD4有联合名义差异，独立时间方向不稳；既往GSDMD来源与因果问题仍未解决。'],['IKZF2','未入当前功能交集','Treg稳态/存活研究不能自动替代内源性线粒体凋亡或焦亡直接依据。'],['GSDMD','未入本轮三项交集','功能轴成立，但本轮匹配状态的巨噬入口未达预设条件；旧低UMI/背景问题亦未消除。']],[70,135,285])
add('恢复候选不意味着恢复旧机制主张。当前55基因保留了疾病响应、功能依据和模型入口；真正的来源特异性、近端信号及死亡执行仍属于后续实验检验对象。')
page('9  分样本表达支持筛选，背景敏感性不替代细胞鉴定')
img('expression.png',w=490)
add('实心点为原始计数CP10k，空心点为背景负荷≤1%的BG1子集CP10k；每点是一只小鼠。五个候选的重点范围均保留方向，图中没有把细胞当作重复进行检验。完整阳性细胞数、≥2 UMI、T标记共检出和背景负荷逐样本列于附件。')
add('BG1是在冻结细胞集合内按预设异源RNA负荷做的描述性子集分析，不是真正的环境RNA去污染，也不能凭此确认每个低丰度RNA的来源。GSDMD和髓系标记类候选仍保留来源风险；不会把过滤后无法评估写成“通过”。')
add('竞争来源比较使用全肺Broad注释；通讯使用冻结AllMac。它们用于不同层面的证据，不拼成同一来源范围的显著因果链。该处理回应此前“AMlike配体差异与AllMac相关性跨范围拼接”的问题。')
page('10  创新性以新问题分层，而非检索不到就通过')
add('对55基因完成165项Europe PMC题名/摘要检索，分别覆盖矽肺/硅尘、内源性凋亡/焦亡、CD4/Treg死亡情境；补检4-1BB、SnoN、caspase-1、IRE1、PP2Cm等常用别名。每项保留查询式、时间、命中数及前10条记录；另对重点候选补充网页检索和原始研究。')
add('本轮排查了同名词误命中、预印本、不同细胞来源和既往已做机制。APP等短符号会命中非基因研究；检索计数只是检索记录，不是创新性得分。FOXP1既有Treg论文未必在限制了死亡关键词的题名/摘要查询中出现，因此结合原文补查。')
table(['分层','处理方式'],[['优先验证问题（4）','IL1R2、DDIT4、SKIL、FOXP1：明确指出与既有研究不同的接收细胞/疾病/来源或功能连接，并记录当前缺口。'],['备选（1）','GSN：保留方向和线粒体功能依据，但关注对照异质性、蛋白形式和跨物种差异。'],['已知通路/疾病参照','例如TNFRSF9、STAT1、CD44、S100A8、TNFRSF1A等；不因有文章一律删除，也不重复包装已知问题。'],['其他三项交集池','完整保留其P/FDR、状态、功能注释和通讯路径；本轮尚未形成足够清晰的优先新机制问题。']],[125,365])
add('创新性结论为“具有可区别于已知工作的验证问题”，不是全球首次声明。未检出文献、题名/摘要未覆盖全文、仅有模型预测等情况明确保留。表中5个优先/备选不是经统计检验得到的第四个硬阈值集合。')
page('11  历史候选逐项回溯解释了去留')
add(f'对既往11基因基准、27基因说明/展示以及新增历史主线共{len(hist)}个基因建立回溯表。{int(hist.three_axis.sum())}个恢复进入本轮三项交集：'+ '、'.join(hist.loc[hist.three_axis,'gene'])+'。全量历史高效应表及现全基因结果在附件可追溯。')
show=['IL1R2','TNFRSF9','FOXP1','GSN','DUSP5','STAT1','ZEB2','IL10','IKZF2','GSDMD','PENK','TTC39C','MDFIC','TFF1','SPP1','CCL5','NTRK3','LCN2']
rows=[]
for g in show:
 z=hist[hist.gene==g]
 if not len(z):continue
 r=z.iloc[0];rows.append([g,'是' if r.DE_axis else '否/NA','是' if r.function_axis else '待证据','是' if r.three_axis else '否',str(r.reason).split('；')[0]])
table(['基因','差异','功能','三项交集','本轮去留原因'],rows,[62,50,60,60,258])
add('“未纳入功能轴”仅指当前冻结的内源性凋亡/焦亡范围和已核查文献，不等于否认泛凋亡、炎症、迁移或代谢作用。高效应的免疫球蛋白或上皮标记不会因效应大而直接升级成CD4损伤目的基因。','small')
page('12  复核范围与下一步实验的判定目标')
add('内部复核已完成：17个差异检验家族逐一独立重算BH；12套pseudobulk计数核查整数、非负、样本顺序、映射后总量不超过原始library；155项三项交集逐一核查P、功能成员及匹配状态的通讯路径。结果见QA_statistical_integrity和交集逐项复核。')
add('当前最有价值的下一步不是再次围绕单一旧基因追加相关性，而是验证本轮可区分的机制问题：')
table(['方向','首轮验证目标','什么结果会改变优先级'],[['IL1R2/Treg','确认蛋白/胞内定位；IL1R2-ENO1及死亡执行终点','若仅为表面活化标记、无ENO1或死亡效应，降为表型标记。'],['DDIT4/ActCD4','比较巨噬条件培养/配体处理与DDIT4干预；线粒体及焦亡终点','若只有一般应激变化而无线粒体或死亡终点，收窄到应激适应。'],['SKIL/Treg','SKIL干预结合p53和内源性凋亡终点；检验CCL2/CXCL10入口','若配体不改变SKIL，保留接收端候选但移除该上游链。'],['FOXP1/CD4','纯化CD4和状态内确认下降；区分巨噬、上皮/内皮来源作用','若仅有组成/活化解释或上游来自其他细胞，修正来源命题。'],['GSN/ActCD4','区分胞内与分泌形式；复核对照异质性和线粒体保护','若效应依赖单个对照或蛋白形式不符，优先级下降。']],[92,216,182])
add('方法边界集中说明：小样本发现属于探索性；分层联合模型假设可共享疾病效应，n=1时间点不能独立检验；时间队列无同期盐水且TLike非纯CD4；捕获的活细胞转录组不直接测量死亡执行；网络不提供因果方向；背景RNA、来源竞争及功能跨细胞类型外推仍需验证。上述边界不取消已观察到的方向和差异。')
add('交付完成标准：新框架已执行，完整候选与去留可复核，关键矛盾和阴性保留。此文件是供再次独立评审的内部版本，不是经外部审核通过的客户机制报告。','small')
# Full gene pool, concise evidence at representative comparison; full multiple comparisons retained elsewhere.
ordered=summary.sort_values(['four_axis_priority','gene'],ascending=[False,True]).reset_index(drop=True)
for start in range(0,len(ordered),19):
 page(f'附录A  完整55基因交集池（{start+1}-{min(start+19,len(ordered))}）')
 add('代表对比仅用于展示：优先表使用指定重点状态；其余优先展示联合模型，再按P排序选一项。它不是新的合并P值。全部155项对比及原检验家族保留在附件。','small')
 rows=[]
 for r in ordered.iloc[start:start+19].itertuples():rows.append([r.gene,r.scope,f'{r.logFC:+.2f}',num(r.PValue),num(r.FDR),'受体' if r.any_direct_receptor else '预测靶',str(r.review_priority)])
 table(['基因','范围','log2FC','P','FDR','入口','定位'],rows,[68,55,53,65,65,62,122])
 add('标记“受体”表示该基因至少有一个状态具有直接受体/诱饵结合入口；当前代表对比的精确入口需看155项表。所有“预测靶”均为模型支持，不能称为已验证调控。时间证据只有CD44、FOXP1、HIF1A、TNFRSF1B进入三项交集。','small')
page('附录B  主要方法与原始文献')
refs=[('方法阈值与NicheNet','https://github.com/saeyslab/multinichenetr/blob/main/vignettes/pairwise_analysis_MISC.Rmd'),('NicheNet基础分析','https://github.com/saeyslab/nichenetr/blob/master/vignettes/ligand_activity_geneset.md'),('CellChatDB.mouse版本快照','https://github.com/jinworks/CellChat/blob/main/data/CellChatDB.mouse.rda'),('Reactome intrinsic pathway','https://reactome.org/content/detail/R-MMU-109606'),('Reactome pyroptosis','https://reactome.org/content/detail/R-MMU-5620971'),('IL1R2-ENO1 / macrophage pyroptosis','https://pmc.ncbi.nlm.nih.gov/articles/PMC12533375/'),('DDIT4 / ROS-TXNIP-NLRP3','https://pubmed.ncbi.nlm.nih.gov/38771134/'),('REDD1 / podocyte pyroptosis','https://pubmed.ncbi.nlm.nih.gov/39920111/'),('SnoN-p53 / stress response','https://pmc.ncbi.nlm.nih.gov/articles/PMC3897278/'),('FOXP1 / BIK','https://pmc.ncbi.nlm.nih.gov/articles/PMC3741507/'),('FOXP1 / Treg apoptosis context','https://pmc.ncbi.nlm.nih.gov/articles/PMC7534899/'),('FOXP1 / Treg apoptosis unchanged','https://pmc.ncbi.nlm.nih.gov/articles/PMC6534289/'),('GSN / neuronal mitochondrial apoptosis','https://doi.org/10.1016/j.mcn.2003.09.012'),('GSN / human-mouse difference','https://pubmed.ncbi.nlm.nih.gov/11039896/'),('DUSP5 / SCLC pyroptosis','https://pubmed.ncbi.nlm.nih.gov/40644541/'),('STAT1 / intrinsic apoptosis','https://pmc.ncbi.nlm.nih.gov/articles/PMC4374929/'),('4-1BB / CD8 survival','https://pubmed.ncbi.nlm.nih.gov/12391199/'),('4-1BB / prior silicosis intervention','https://pmc.ncbi.nlm.nih.gov/articles/PMC5039680/'),('4-1BB / prior macrophage silicosis mechanism','https://pubmed.ncbi.nlm.nih.gov/30250465/')]
for i,(label,url) in enumerate(refs,1):add(f'[{i}] {label}\n{url}','small')
add('完整查询、命中题名/摘要、原始GO证据、数据库哈希及R包版本均随复核附件保存；不能只凭此页参考文献列表判断所有候选已完成全系统综述。','small')
page('附录C  文件对应关系与复算位置')
table(['内容','附件位置/文件'],[['核心交集与优先问题','final_tables/四轴候选总表_55基因.tsv；优先验证候选_5基因.tsv'],['全对比及历史回溯','final_tables/全部155项候选对比.tsv；历史候选逐项回溯.tsv'],['完整检验家族','tables/DE_all_models.tsv；各模型独立表；model_diagnostics.tsv'],['功能与通讯','tables/function_membership_extended.tsv；candidate_complete_routes.tsv；NicheNet_activity_all.tsv；NicheNet_top250_targets.tsv'],['表达与来源','tables/candidate_sample_expression_QC.tsv；whole_lung_source_expression.tsv.gz；communication_sample_coverage.tsv'],['创新性与原始记录','final_tables/创新性检索记录_165项.tsv；创新性文献命中与摘要.tsv；literature_success/'],['冻结方案及资源版本','protocol_frozen.json；resources/；code/；logs/'],['完整性复核','tables/QA_statistical_integrity.tsv；final_tables/交集逐项复核.tsv；MANIFEST_SHA256.json']],[125,365])
add('服务器项目目录：\n/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/three_cohort_four_axis_20261017_v16','small')
add('原始矩阵、细胞级计数、edgeR拟合对象及中间输入保留于服务器inputs目录；未向下载目录复制细胞矩阵。附件提供报告复核所需统计证据与代码，原始GEO数据和项目冻结输入可用于全流程复算。','small')
add('历史S17部分服务器文件名存在既有编码异常，本轮按实际路径只读读取，保存原始路径，不改写原文件名或内容。公开资源下载曾改由本机传入；文献API失败日志保留，最终165条查询全部成功，未把失败当零命中。','small')
add('本报告与附件为v16基准；后续独立评审应针对该版本提出问题编号。','small')
pdf=D/'v16_三队列四轴重筛报告.pdf'
def footer(canvas,doc):
 canvas.saveState();canvas.setFont('CN',8);canvas.setFillColor(colors.HexColor('#536775'));canvas.drawString(48,27,'v16 | 三队列四轴重筛 | 内部复审版 | 2026-10-17');canvas.drawRightString(A4[0]-48,27,str(doc.page));canvas.restoreState()
doc=SimpleDocTemplate(str(pdf),pagesize=A4,rightMargin=48,leftMargin=48,topMargin=43,bottomMargin=43,title='v16 三队列四轴重筛报告',author='项目分析');doc.build(story,onFirstPage=footer,onLaterPages=footer)
(O/'report_output.json').write_text(json.dumps({'pdf':str(pdf),'delivery':str(D)},ensure_ascii=False,indent=2),encoding='utf-8');print(pdf)