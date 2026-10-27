from pathlib import Path
import json,shutil
O=Path(__file__).parent;B=O.parent/'lung_screening_display_20261017_v24';base=(B/'report.py').read_text(encoding='utf-8');shutil.copy2(B/'editorial_replacements.json',O/'editorial_replacements.json');ms=json.loads((B/'figure_manifest.json').read_text(encoding='utf-8'))
for m in ms:
 if m['id'] in ['F05','F06']:m['file']=str(O/'figures'/(m['id']+'.png'))
 if m['id']=='F01':m['title']='GSE314943主队列的全肺细胞图谱'
 if m['id']=='F03':m['title']='GSE314943主队列的三类核心CD4状态'
 if m['id']=='F05':m['caption']='联合差异分析使用GSE314943的6个样本及GSE183682的4个样本，共10个样本（5个对照、5个矽尘）。edgeR模型控制队列-时间层；红为FDR<0.05，黄为仅原始P<0.05，灰为其余。基因名位于绘图区外，连线端点对应实际基因坐标。全肺图谱来自主队列，补充队列图谱见前文示意图B。'
 if m['id']=='F06':m['caption']='蓝色为GSE314943的6样本疾病效应，橙色为GSE314943与GSE183682合计10样本的分层疾病效应；二者包含重叠样本。TTC39C联合DE P=0.163、FDR=0.803，其程序关联与疾病差异证据分别评价。'
(O/'figure_manifest.json').write_text(json.dumps(ms,ensure_ascii=False,indent=2),encoding='utf-8')
base=base.replace('矽肺_CD4候选研究_筛选过程修订_v24_20261017','矽肺_CD4候选研究_队列衔接修订_v25_20261017').replace('矽肺免疫全景_CD4候选与验证_报告_v24.pdf','矽肺免疫全景_CD4候选与验证_报告_v25.pdf')
# Keep numbered scientific Figures1-24 stable; add three clearly labelled explanatory displays.
helper="""
def overview_image(id,maxh=520):
 fn=O/'figures'/(id+'.png');im=PILImage.open(fn);k=min(505/im.width,maxh/im.height);story.append(Image(str(fn),width=im.width*k,height=im.height*k,hAlign='CENTER'));story.append(Spacer(1,7))
def cohort_pages():
 page();p('队列关系：6个主队列样本与4个补充样本构成联合分析','h2')
 p('“联合10样本”指GSE314943与GSE183682的样本级分层模型，共5个对照和5个矽尘样本。两个队列的细胞图谱分别展示；联合统计在相同CD4状态内汇总计数，并控制队列与时间层。')
 overview_image('S_A')
 p('示意图A。主队列GSE314943用于全肺及核心CD4图谱；GSE183682提供D7和D56的4个补充样本。三个模型层为GSE314943-D140、GSE183682-D7和GSE183682-D56。GSE274372的12个时间样本单独分析，不计入该10样本模型。','small')
 page();p('补充队列GSE183682：联合分析中4个样本的CD4图谱','h2')
 p('GSE183682纳入1,413个Core CD4细胞，包括742个Naive、398个ActCD4和273个Treg。下图按已有队列内坐标展示细胞状态，并列明参与联合模型的全部4个样本及细胞数量。')
 overview_image('S_B')
 p('示意图B。GSM5567501/502分别为D7对照/矽尘，GSM5567503/504分别为D56对照/矽尘。柱图为细胞数量描述，不以细胞数代替样本量；每时间点1+1，不作单独推断性检验。两队列UMAP的坐标独立计算，不能直接比较坐标位置或距离。','small')
def early_framework():
 page();p('差异分析前的筛选框架：三项证据交集与后续分层','h2')
 p('差异分析提供候选变化证据，随后结合目标损伤/死亡功能注释及同状态巨噬入口。三项证据交集得到55个基因，再与补充候选形成90个基因的评价集合。创新性与实验可行性参与后续分层，不作为第四个统一硬门槛。')
 overview_image('S_C')
 p('示意图C。圆形交叠仅用于说明筛选逻辑，面积不代表集合大小；55为按细胞状态和对比匹配后去重的三项交集基因数，不完整列示各互斥区域数量。功能维度来自GO/Reactome注释及直接功能文献，不以本项目GO富集P值替代单基因功能依据。逐步阈值、集合合并和候选去向详见后文筛选章节。','small')
"""
base=base.replace("p('矽肺免疫微环境与\\nCD4⁺ T细胞候选研究'",helper+"\np('矽肺免疫微环境与\\nCD4⁺ T细胞候选研究'",1)
base=base.replace("page();p('分析方法','h')","cohort_pages()\npage();p('分析方法','h')",1)
base=base.replace("for m in manifest:\n if m['id']=='F21':screening_methods()","for m in manifest:\n if m['id']=='F05':early_framework()\n if m['id']=='F21':screening_methods()")
base=base.replace("p('图示阅读顺序','h2');p('图1-2", "p('图示阅读顺序','h2');p('示意图A/B说明队列关系与补充队列图谱；示意图C在差异分析前介绍交集框架。图1-2")
base=base.replace("'figures':len(manifest)","'figures':len(manifest),'explanatory_displays':3,'total_visuals':len(manifest)+3")
(O/'report.py').write_text(base,encoding='utf-8');extra=[{'id':'S_A','title':'图谱范围与联合统计范围','role':'cohort relationship'},{'id':'S_B','title':'GSE183682补充队列CD4图谱','role':'existing cohort-specific UMAP and sample counts'},{'id':'S_C','title':'差异分析前的三项交集框架','role':'conceptual overlap; not area-proportional quantitative Venn'}];(O/'extra_manifest.json').write_text(json.dumps(extra,ensure_ascii=False,indent=2),encoding='utf-8');print('REPORT_PREPARED')