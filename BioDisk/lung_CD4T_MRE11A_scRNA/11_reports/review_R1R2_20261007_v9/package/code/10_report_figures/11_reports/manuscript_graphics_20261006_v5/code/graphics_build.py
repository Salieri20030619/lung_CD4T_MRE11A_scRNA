from pathlib import Path
import json,re
B=Path(__file__).parent
src=(B/'narrative_document.py').read_text(encoding='utf8').replace("D=B/'narrative_delivery'","D=B/'graphics_delivery'")
exec(compile(src,'graphics_document_foundation','exec'))
G=json.loads((D/'provenance/GSEA_figure_manifest.json').read_text(encoding='utf8'))
h('1.10　UMAP与GSEA可视化')
textblock('''UMAP直接使用冻结Seurat对象已有的二维坐标。全肺对象为GSE288887整合集合，含其复用的GSE135893与GSE227136来源；沿用原广谱标记注释lineage_prelim。全肺与免疫子集使用相同坐标，CD4亚群与来源面板使用原CD4对象的另一套坐标；本轮未重跑降维、整合或聚类。CD4细胞标识、供者、疾病组、来源和亚群均与疾病分析清单逐项核对。点云用于展示细胞组成和来源分布，疾病比较仍以供者模型为依据。
GSEA曲线直接读取原保存的gseaResult对象，使用原来源校正伪合并排序sign(logFC)×√QLF、基因集及检验结果。分别在Core和Primary的GO-BP、KEGG完整已检验条目中，按|NES|降序选取5条；并列时按通路ID排序，不按P值或BH筛选，也不因名称剔除核糖体等条目。enrichplot::gseaplot2绘制加权running ES、各通路命中条带、完整排序统计量三层，权重指数保持原值1；rel_heights=(1.5,0.4,1)、subplots=1:3、pvalue_table=FALSE、ES_geom='line'、base_size=11。各曲线的极值ES与原enrichmentScore核对，命中数与原setSize核对。独立图为10×8英寸PDF及300 dpi PNG，图下保留原NES、P和BH；BH家族仍为相应分析范围×数据库的全部原检验条目，不对展示的5条重新校正。''')

base_p=p;base_h=h;base_fig=fig
old_to_new={1:'2',2:'3',3:'4',4:'5',5:'6',6:'S2',7:'9',8:'10'}
def refs(s):
 s=re.sub(r'图([1-8])(?!\d)',lambda m:'图'+old_to_new[int(m.group(1))],s)
 return s.replace('{NEW_GSEA}','图7、图8').replace('{OLD_SUMMARY}','图S2').replace('{PRIMARY_GSEA}','图S3、图S4').replace('{NEW_UMAP}','图1')
def p(t,style=None):return base_p(refs(t),style)
def newfig(file,title,caption,supp=False):
 global p
 previous=p;p=base_p
 try:return base_fig(file,title,caption,supp)
 finally:p=previous

umap_caption='A，全肺冻结UMAP含20供者、25份样本、192575个细胞，颜色为原广谱标记初步注释：Epi，上皮；Endo，内皮；Str，基质；My，髓系；T，T细胞；NK，自然杀伤细胞；B，B细胞；Mast，肥大细胞。B，在A的同一坐标与坐标范围内展示57706个免疫细胞；My为广谱髓系标签，未将其全部改称巨噬细胞。C，独立冻结CD4嵌入含19供者、24份样本、6523个细胞，CCR7/GZMA/IL7R分别对应CD4-CCR7/CD4-GZMA/CD4-IL7R，另含Treg和CD4state。D，与C相同的坐标按来源着色：G135=GSE135893，G227=GSE227136，G288=GSE288887新样本。A/B来源为原全肺scored对象，C/D来源为CD4_final对象；人类全肺与CD4使用各自嵌入，未重新整合。点的空间分离不作为疾病显著性检验。'
def h(t,level=2):
 if t.startswith('2.2'):
  newfig('Human_UMAP_atlas_CD4.png','人类肺组织全景、主要免疫类群与CD4亚群的真实UMAP',umap_caption)
 return base_h(t,level)

translation={
 'KEGG_RIBOSOME':'核糖体','KEGG_OXIDATIVE_PHOSPHORYLATION':'氧化磷酸化（OXPHOS）',
 'KEGG_PARKINSONS_DISEASE':'帕金森病基因集','KEGG_CARDIAC_MUSCLE_CONTRACTION':'心肌收缩','KEGG_HUNTINGTONS_DISEASE':'亨廷顿病基因集',
 'GOBP_NEGATIVE_REGULATION_OF_SYNCYTIUM_FORMATION_BY_PLASMA_MEMBRANE_FUSION':'质膜融合形成合胞体的负调控',
 'GOBP_REGULATION_OF_SYNCYTIUM_FORMATION_BY_PLASMA_MEMBRANE_FUSION':'质膜融合形成合胞体的调控',
 'GOBP_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT':'ATP合成耦联电子传递','GOBP_OXIDATIVE_PHOSPHORYLATION':'氧化磷酸化（OXPHOS）',
 'GOBP_PROTON_MOTIVE_FORCE_DRIVEN_ATP_SYNTHESIS':'质子驱动力驱动的ATP合成','GOBP_CYTOPLASMIC_TRANSLATION':'细胞质翻译'
}
def gsea_caption(nm):
 z=G[nm];scope='核心CD4（排除CD4state）' if nm.startswith('Core') else '全体CD4（保留CD4state）';db='KEGG' if nm.endswith('KEGG') else 'GO-BP'
 names='；'.join(r['display_code']+'，'+translation[r['ID']] for r in z['selected'])+'。'
 return f'{scope}的{db}原来源校正分析，{z["ranked_genes"]}个排序基因、{z["tested"]}个原检验条目。展示完整结果中|NES|最大的5条，不等同P值最小的5条；未设置显著性展示门槛，本图所选5条均通过原BH<0.05。上层为真实加权running ES，中层为逐通路命中条带，下层为完整排序统计量；正排序表示CWP方向。颜色与通路代码一致，图下列原NES、P及BH，BH仍覆盖全部{z["tested"]}条。'+names+('KEGG疾病名称为基因集名称，不用于判断患者伴随疾病。' if db=='KEGG' else 'reg.表示调控，-reg.表示负调控；通路名称描述基因集，本图不单独判定该生物过程已发生。')

summary_caption=None
def fig(file,title,caption,supp=False):
 global summary_caption
 if file=='fig2_pathways.png':
  summary_caption=refs(caption)
  newfig(G['Core_KEGG']['stem']+'.png','人类核心CD4的KEGG Top5多通路富集曲线',gsea_caption('Core_KEGG'))
  newfig(G['Core_GO']['stem']+'.png','人类核心CD4的GO-BP Top5多通路富集曲线',gsea_caption('Core_GO'))
 else:newfig(file,title,refs(caption),supp)

src=(B/'narrative_results.py').read_text(encoding='utf8')
src=src.replace('各范围的供者和细胞覆盖见图1及表1。','真实UMAP见{NEW_UMAP}，各范围的供者和细胞覆盖见图1及表1。')
src=src.replace('图6将GSEA、供者GSVA及既有剔除敏感性并列展示。','GSEA曲线见{NEW_GSEA}，供者GSVA及既有剔除敏感性结果见{OLD_SUMMARY}，全体CD4曲线见{PRIMARY_GSEA}。')
exec(compile(src,'graphics_results','exec'))
supp=(B/'manuscript_supplement.py').read_text(encoding='utf8').replace('方法1.6','方法1.9').replace('不仅本表7程序','不仅所展示的7个程序')
supp=supp.replace('doc.add_page_break();','')
exec(compile(supp,'graphics_supplement','exec'))
base_h('补充图：通路汇总与全体CD4富集曲线',1)
newfig('fig2_pathways.png','通路富集汇总、供者GSVA与呼吸链敏感性',summary_caption,supp=True)
for nm in ['Primary_KEGG','Primary_GO']:
 doc.add_page_break()
 newfig(G[nm]['stem']+'.png','人类全体CD4的'+('KEGG' if nm.endswith('KEGG')else'GO-BP')+' Top5多通路富集曲线',gsea_caption(nm),supp=True)
p=base_p;h=base_h;fig=base_fig
exec(compile((B/'manuscript_finish.py').read_text(encoding='utf8'),'graphics_finish','exec'))
assert len(figures)==14 and len(tables)==15
