from pathlib import Path
import json,shutil,hashlib
O=Path(__file__).parent;V=O.parent/'lung_client_editorial_20261017_v21';base=(V/'report.py').read_text(encoding='utf-8')
shutil.copy2(V/'editorial_replacements.json',O/'editorial_replacements.json')
r=json.loads((O/'editorial_replacements.json').read_text(encoding='utf-8'));r.append(['图11展示合格相互作用数量；图12展示样本级发送配体表达。','图11通过弦图展示巨噬细胞来源与CD4状态之间满足表达覆盖的相互作用数量，图12展示Treg相关配体-受体候选连接。'])
(O/'editorial_replacements.json').write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf-8')
ms=json.loads((V/'figure_manifest.json').read_text(encoding='utf-8'))
for m in ms:
 m['file']=str(O/'figures'/(m['id']+'.png'))
 if m['id'] in ['F07','F08']:m['caption']+=' 点旁文字标注原始P值；完整家族FDR保留于统计表。'
 if m['id'] in ['F09','F10']:
  sc='Treg' if m['id']=='F09' else 'ActCD4';m['title']=sc+'的KEGG Top 5多通路富集曲线';m['caption']='采用enrichplot::gseaplot2绘制运行富集曲线、基因命中位置及排名指标三层图。在原始P<0.05的KEGG结果中按|NES|降序展示前5项，同值按通路ID排序；图下列出各通路原始P与完整家族FDR。图例位于独立区域，颜色对应同名通路；排名仍为sign(logFC)×sqrt(QL F)，全部统计值来自原分析。P估计下限为1e-8。'
 if m['id']=='F11':
  m['title']='巨噬细胞至CD4状态的表达相容通讯候选';m['caption']='按发送细胞范围和疾病分组分别绘制弦图，方向为巨噬细胞至CD4状态。弦宽表示面板内合格相互作用数量的相对分布，面板下列出各状态的绝对数量；各面板圆周分别归一化，不宜跨面板直接比较绝对弦宽。AllMac与AMlike为重叠来源，分别展示；Core与其三个状态不重复累加。PROS1路线未纳入。数量不等于实测通讯强度。'
 if m['id']=='F12':
  m['title']='Treg相关配体-受体候选连接';m['caption']='展示已评价路线中的6条Treg配体-受体候选关系，每条关系赋相同权重。橙色为发送端配体，蓝色为Treg受体；箭头表示配体至受体的关系方向，不表示实测激活。IL1R2为诱饵受体，其结合关系不等同于胞内信号传导。受体必要组分及来源范围见完整路线表。'
(O/'figure_manifest.json').write_text(json.dumps(ms,ensure_ascii=False,indent=2),encoding='utf-8')
base=base.replace('矽肺_CD4候选研究_正式报告_v21_20261017','矽肺_CD4候选研究_纵向修订_v22_20261017').replace('矽肺免疫全景_CD4候选与验证_正式报告_v21.pdf','矽肺免疫全景_CD4候选与验证_纵向报告_v22.pdf')
base=base.replace("template='portrait' if m['id'] in ['F15','F16'] else 'landscape'","template='portrait'").replace("maxw=505 if template=='portrait' else 752;maxh=465 if template=='portrait' else 315","maxw=505;maxh=540")
update="""
results['F09']='Treg的KEGG富集曲线同时呈现疾病上调端与下调端的通路分布。核糖体相关集合偏向下调端；其他通路需结合其具体成员及细胞状态解释。'
results['F10']='ActCD4的多通路曲线显示不同功能集合在完整基因排序中的富集位置。正、负NES分别对应疾病上调端和下调端，图中并列报告每条通路的P值及FDR。'
results['F11']='巨噬细胞与Naive、ActCD4和Treg之间均存在满足所设表达覆盖的配体-受体候选。分组弦图展示这些候选的状态分布，并保留各发送范围的独立定义。'
results['F12']='Treg候选关系包括IL1A/IL1B-IL1R2、CCL2-CCR2、CXCL10-CXCR3、CD200-CD200R1及SEMA3E-PLXND1，可用于安排后续来源、受体与近端响应验证。'
"""
base=base.replace('for m in manifest:\n template=',update+'\nfor m in manifest:\n template=')
(O/'report.py').write_text(base,encoding='utf-8')
ref=Path(r'D:\Project_0106_scripts_20260430\Project_0106_scripts_20260430\r.09_GSEA.r');shutil.copy2(ref,O/'reference_r.09_GSEA.r')
plan={'change_scope':['all portrait pages','legend non-overlap','GO/KEGG raw P annotations','circlize communication chords','gseaplot2 top5 overlays using supplied reference'],'analysis_inputs':'v20 saved results; no differential/enrichment rerun','GSEA_selection':'raw P<0.05, descending abs(NES), ID tie-break, 5 KEGG pathways per state','GSEA_reference':str(ref),'reference_sha256':hashlib.sha256(ref.read_bytes()).hexdigest(),'GSEA_reference_adaptation':'mouse saved statistics and ranks; retains three panels and relative heights1.5/0.4/1; external legend and statistic table; no GSE16561/human annotation applied'}
(O/'change_scope.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2),encoding='utf-8')
print('REPORT_PREPARED')