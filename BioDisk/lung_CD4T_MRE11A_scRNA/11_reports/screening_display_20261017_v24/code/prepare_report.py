from pathlib import Path
import json,shutil
O=Path(__file__).parent;B=O.parent/'lung_chord_clarity_20261017_v23';base=(B/'report.py').read_text(encoding='utf-8');shutil.copy2(B/'editorial_replacements.json',O/'editorial_replacements.json');ms=json.loads((B/'figure_manifest.json').read_text(encoding='utf-8'))
for m in ms:
 if m['id']=='F05':m['file']=str(O/'figures/F05.png');m['caption']+=' 标记基因采用绘图区外的独立标签列，连线端点对应实际基因坐标。'
old21=ms[20].copy();old22=ms[21].copy();old21['id']='F23';old22['id']='F24'
new=[{'id':'F21','module':'候选筛选','title':'三项证据交集与补充候选合并构建评价集合','file':str(O/'figures/F21.png'),'caption':'依次取目标效应名义差异、目标损伤/死亡功能与同状态巨噬入口的交集。数字同时列去重基因数与基因-状态-对比记录数。55个交集基因与42个研究问题相关补充候选合并，重叠7个，形成90个综合评价基因；该步骤是并集扩展。补充候选不是额外独立队列，交集不产生新的联合P或集合FDR。'}, {'id':'F22','module':'候选筛选','title':'资格分层与实验推荐分别依据证据完整度和研究价值','file':str(O/'figures/F22.png'),'caption':'90个基因按同状态证据、细胞来源及作用角色综合评价，形成52个计算支持、37个待补和1个未达到名义差异条件的基因。待补项在首次不满足相应条件时计数，不重复累加。9个重点基因来自完整90基因集合，包括5个计算支持及4个待补；FOXP1涉及两个状态，共10个组合。创新性、功能文献与实验可行性用于形成预实验建议，不作为统一二元交集。'}]
ms=ms[:20]+new+[old21,old22];(O/'figure_manifest.json').write_text(json.dumps(ms,ensure_ascii=False,indent=2),encoding='utf-8')
base=base.replace('assert len(manifest)==22','assert len(manifest)==24').replace('矽肺_CD4候选研究_摘要图形修订_v23_20261017','矽肺_CD4候选研究_筛选过程修订_v24_20261017').replace('矽肺免疫全景_CD4候选与验证_报告_v23.pdf','矽肺免疫全景_CD4候选与验证_报告_v24.pdf')
base=base.replace('图19-20 NicheNet；图21-22候选与验证。','图19-20 NicheNet；图21-22候选筛选流程；图23-24候选与验证。')
# Existing candidate panels become Figures23/24; add two screening result panels.
add="""
results['F23']=results['F21']
results['F24']=results['F22']
results['F21']='主要差异模型共涉及9,091个去重基因；2,765个在至少一个目标对比中达到原始P<0.05。与目标功能集合相交后保留111个，再匹配相同接收状态的巨噬入口得到55个。结合补充候选形成90基因评价集合。'
results['F22']='综合评价最终支持52个计算候选，37个保留为证据待补，1个未达到当前名义差异条件。重点9基因不是52基因的继续硬筛子集；TTC39C虽然处于功能证据待补层，仍可凭凋亡程序关联及功能实验先例进入探索性预实验。'
def screening_methods():
 page();p('候选筛选：证据维度、阈值与集合关系','h')
 p('筛选首先整合差异、目标损伤/死亡功能和巨噬来源入口三项证据，再进行细胞身份与作用角色评价。创新性及实验可行性用于确定研究方向，不要求所有候选在四项指标上同时通过一个统一阈值。')
 tab([['维度','实际判定条件','参考证据与解释'],['差异表达','有效样本级模型的目标效应原始P<0.05；按基因-状态-对比保留记录','FDR、|log2FC|≥0.5、跨队列方向及敏感性用于评价强度；不要求所有队列同时显著'],['目标损伤/死亡功能','GO/Reactome内源性凋亡、线粒体损伤相关过程或焦亡集合；明确的直接功能文献可扩展','两类机制取并集；不要求基因同时参与凋亡和焦亡。一般凋亡/状态证据不足以证明目标细胞自身损伤'],['巨噬来源入口','配体和完整受体组分在同一条件及接收状态内满足覆盖；候选为受体组分或具有对应配体-靶先验','每样本至少10个目标细胞，阳性比例≥5%，且组内至少一半合格样本满足。网络先验不等于功能激活'],['创新性与可行性','逐项比较已有机制、细胞背景、研究空间与可观察实验表型','作为分层依据。已有发表不自动排除，未检索到文献不自动证明创新；不产生二元创新性P值']],[80,215,210])
 p('交集按同一细胞状态匹配；时间证据单列','h2')
 p('主要DE范围包括主队列4个细胞范围、联合队列4个范围及时间队列3个时间对比和1个总体时间检验，共12个家族。名义差异在这些目标对比间取并集，时间效应与疾病效应分别解释；功能和网络匹配完成后再按基因去重，避免跨状态拼接。','small')
 p('初筛与综合评价的网络范围','h2')
 p('初筛网络分支采用校正AUPR>0的合格配体，按活动排序取前30个，每配体取表达背景内前250个靶基因；直接受体分支单独保留。综合评价时扩展检查有正向活动的合格配体，网络排名作为参考，不因超过30/250而统一排除。','small')
 p('补充集合、资格分层与预实验推荐','h2')
 p('42个补充候选围绕T细胞状态、细胞死亡及相关调控问题纳入，与55个交集基因合并后统一评价；来源及成员逐项列入附表。细胞来源不清或主要作用于其他细胞的输出效应基因进入待补层。重点分析选择IL1R2、DDIT4、SKIL、FOXP1、GSN、AHNAK、TTC39C、MDFIC、ST6GALNAC3；预实验顺序进一步依据功能文献、程序关联和可操作性确定。','small')
"""
base=base.replace('for m in manifest:\n template=',add+'\nfor m in manifest:\n if m[\'id\']==\'F21\':screening_methods()\n template=')
(O/'report.py').write_text(base,encoding='utf-8');print('MANIFEST',len(ms))