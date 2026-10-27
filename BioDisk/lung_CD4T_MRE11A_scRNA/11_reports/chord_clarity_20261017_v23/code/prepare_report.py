from pathlib import Path
import json,shutil
O=Path(__file__).parent;B=O.parent/'lung_client_layout_20261017_v22';base=(B/'report.py').read_text(encoding='utf-8');shutil.copy2(B/'editorial_replacements.json',O/'editorial_replacements.json');ms=json.loads((B/'figure_manifest.json').read_text(encoding='utf-8'))
for m in ms:
 if m['id'] in ['F11','F12','F17','F20']:m['file']=str(O/'figures'/(m['id']+'.png'))
 if m['id']=='F11':m['caption']+=' 弦色对应接收状态：青色Naive、橙色ActCD4、紫色Treg；宽箭头尖端为接收端。每条弦汇总该来源至相应状态的多条候选相互作用。'
 if m['id']=='F12':m['caption']='展示6条Treg配体-受体候选关系，每条等权。不同配体使用不同颜色；L为发送端配体，R为接收端受体，宽箭头尖端指向R。图下逐条列出配体→受体对应关系。节点按对应关系排列以减少交叉。IL1R2为诱饵受体，结合关系不等同于胞内信号传导；本图未表示实测激活或通讯强度。'
 if m['id']=='F17':
  m['title']='邻域主干辅助显示Naive与ActCD4的状态排序';m['caption']='1,734个细胞的既有K15邻域图为单一连通分量。按DPT分为12个分位组，汇总组间邻接权重并按组大小归一化，取最大生成树形成11条主干边；节点位于最接近组内UMAP中位数的细胞。箭头指向较高DPT。主干为描述性图抽象，不是独立谱系推断；DPT排序、统计及不包含Treg的范围均保持一致。'
 if m['id']=='F20':
  m['title']='证据矩阵区分配体-候选连接及受体关系';m['caption']='行表示配体及括号内受体组分，列为候选基因及目标状态。蓝色为NicheNet靶先验，并标注原靶排名；橙色为直接诱饵受体关系；灰色为本展示范围内无合格连接，不代表所有背景下均不存在。排名不等于效应量或P值。保留原8条代表路线；TTC39C未满足对应联合名义差异筛选条件，列内保持空缺。'
(O/'figure_manifest.json').write_text(json.dumps(ms,ensure_ascii=False,indent=2),encoding='utf-8')
start=base.index("p('矽肺免疫微环境与\\nCD4⁺ T细胞候选研究'");end=base.index("page();p('研究设计与数据范围'",start)
cover="""p('矽肺免疫微环境与\\nCD4⁺ T细胞候选研究','h');p('单细胞图谱、功能程序、细胞通讯及实验候选整合报告','h2');p('2026年10月17日｜研究报告')
p('研究目的','h2');p('本研究围绕“矽肺巨噬细胞环境与CD4⁺ T细胞状态及损伤反应有何联系”展开，利用公开单细胞数据筛选可进入基础实验的候选基因，并确定优先验证的细胞亚群和功能终点。')
p('使用什么数据，完成哪些分析','h2');p('以GSE314943的3个对照和3个矽尘样本为主队列，构建37,155个全肺细胞的图谱，重点分析其中2,311个CD4细胞。另纳入GSE183682的4个样本开展分层联合分析，并以GSE274372的12个时间点样本提供时间背景；联合模型包含主队列，不作为另一次独立验证。')
p('分析依次包括免疫谱系与CD4亚群组成、样本级差异表达、GO/KEGG及GSEA、巨噬细胞至CD4的配体-受体候选、转录因子网络、焦亡/凋亡及线粒体相关程序评分、状态拟时序和NicheNet配体-靶基因预测。结合功能文献，对90个基因进行综合评价，并重点分析9个基因的10个细胞状态组合。')
p('主要发现与推荐方向','h2');p('疾病组呈Naive比例下降、ActCD4和Treg比例上升的方向；对应组成检验未达到P<0.05。综合筛选得到52个符合计算条件的基因，但各候选的证据侧重不同。结合疾病差异、程序关联和既往功能研究，建议优先验证以下3个方向。')
tab([['候选与定位','直接依据','首轮实验终点'],['IL1R2 / Treg\\n优先验证','联合疾病上调：logFC=3.290，P=1.66×10^-6，FDR=0.000662；有Treg功能干预研究','Treg功能及损伤耐受，按结果再评价代谢机制'],['TTC39C / ActCD4\\n损伤方向探索备选','ApoExec名义关联：rho=0.943，P=0.00347，FDR=0.0556；联合疾病差异未达名义显著','存活、凋亡与线粒体状态；先确认干预作用方向'],['SKIL / Treg\\n状态调控备选','TAct名义关联：rho=0.965，P=0.00208，FDR=0.0744；有Treg分化研究','Treg状态与功能，同时检测SnoN蛋白及SKI补偿']],[105,225,175])
p('现有结果支持候选选择及预实验设计，尚未证明“巨噬信号→候选基因→CD4线粒体损伤或焦亡”的完整因果链。下一步用巨噬刺激与候选干预的组合实验，检验候选是否实际改变T细胞表型。','small')
"""
base=base[:start]+cover+base[end:]
base=base.replace('矽肺_CD4候选研究_纵向修订_v22_20261017','矽肺_CD4候选研究_摘要图形修订_v23_20261017').replace('矽肺免疫全景_CD4候选与验证_纵向报告_v22.pdf','矽肺免疫全景_CD4候选与验证_报告_v23.pdf')
insert="""
results['F17']='四套DPT排序的相关为0.904-0.999。叠加的邻域主干连接Naive/ActCD4的主要状态区域，帮助阅读DPT方向；线条来自原始邻域连接的描述性抽象，不作为真实分化过程或损伤进展的证明。'
results['F20']='矩阵直接显示IL1R2的诱饵受体关系与SKIL的配体-靶先验连接。TTC39C在本展示对应的联合差异及网络条件下没有合格连接，其程序关联证据仍单独评价。'
"""
base=base.replace('for m in manifest:\n template=',insert+'\nfor m in manifest:\n template=')
r=json.loads((O/'editorial_replacements.json').read_text(encoding='utf-8'));r.append(['图中4个分位区间只用于展示。','图中4个分位区间只用于展示。UMAP主干将K15邻域图中的细胞按DPT分为12个分位组，汇总组间连接强度并除以两组细胞数乘积的平方根，以最大生成树显示主要邻域连接；箭头沿既有DPT增大方向。'])
(O/'editorial_replacements.json').write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf-8');(O/'report.py').write_text(base,encoding='utf-8');print('PREPARED_REPORT')