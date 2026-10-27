from pathlib import Path
import ast,re
# Preserve the validated document foundation while changing the scientific order.
src=Path('manuscript_document.py').read_text(encoding='utf8')
prefix,body=src.split("h('摘要',1)",1)
prefix=prefix.replace("D=B/'manuscript_delivery'","D=B/'narrative_delivery'")
exec(compile(prefix,'manuscript_document_foundation','exec'))
h('摘要',1)
textblock('''目的：以尘肺患者的免疫细胞图谱为起点，分析CD4⁺ T细胞亚群及其表型变化，识别与接收者疾病程序相符的巨噬细胞信号，并分别评价IKZF2、STAT1及其相关功能程序。
方法：复用已核验的单细胞对象和既有结果。人类CD4 RNA队列含19名供者（10名对照、9名煤工尘肺患者），以供者为统计单位分析组成及程序差异；巨噬细胞通讯结合配体–受体表达、NicheNet配体活性和同供者关联。候选表达与SCENIC调控子活性分别评价，使用三个随机种子、去重叠基因及供者Bootstrap进行敏感性分析。GSEA、GSVA、UCell和表达匹配评分用于刻画功能程序。小鼠主队列每组3个文库，用于亚群、程序和信号方向的对照。
结果：人类CD4细胞分为CCR7、GZMA、IL7R、Treg及CD4state五个亚群。GZMA调节性相关程序和CCR7抑制相关程序呈名义疾病关联；来源与深度校正后P分别为0.008304和0.015698，对应12程序BH为0.099645和0.188374。巨噬信号分析中，IL15在人类核心CD4及IL7R的名义上调靶集均排名1，校正AUPR为0.012640和0.012361；GZMA的FDR上调靶集则以CLSTN1、ADAM15等候选居前。小鼠核心CD4的FDR上调靶集中Il15排名3。两个转录因子呈不同证据特征：人类GZMA亚群IKZF2 RNA的来源校正效应为+0.044530（P=0.002109，FDR=0.025304），加入深度后FDR=0.051376；STAT1在GZMA亚群与炎症小体程序的三个种子均为名义正相关（P=0.008587–0.030860，FDR=0.256882–0.344046）。小鼠核心CD4的两候选调控子活性均呈上调。功能程序方面，人类核心CD4的KEGG氧化磷酸化GSEA为NES=2.774、FDR=2.06×10⁻¹³，去核糖体及线粒体编码基因后仍保留富集。IKZF2与焦亡执行相关评分在三个种子均负相关（n=17，ρ=−0.830、−0.838、−0.664；P=0.000446、0.000349、0.013365；原BH=0.02567、0.02883、0.11077），并伴随凋亡执行评分的名义负相关；其背景敏感性和各评分方法的差异在讨论中统一评估。
结论：人类患者数据呈现亚群特异的CD4转录状态、巨噬信号优先级和呼吸链程序富集，IKZF2与STAT1分别提供不同的表达及调控关联。小鼠证据补充了相应亚群和程序的方向比较。这些发现为后续选择细胞范围、候选和功能终点提供依据，研究结论限定于观察性转录及通讯证据。''')
p('关键词：煤工尘肺；CD4 T细胞；巨噬细胞；IKZF2；STAT1；GSEA；SCENIC；NicheNet','Caption')
method_source=src[src.index("h('1　研究背景与方法',1)"):]
parts=re.split(r"(?=^h\('1\.\d+)",method_source,flags=re.M)
exec(compile(parts[0],'methods_heading','exec'))
blocks={re.search(r"h\('1\.(\d+)",x).group(1):x for x in parts[1:]}
order=['1','2','3','4','9','7','8','5','6']
for new,old in enumerate(order,1):
 b=blocks[old]
 b=re.sub(r"h\('1\.\d+",f"h('1.{new}",b,count=1)
 b=b.replace('分析沿“疾病表型—细胞内程序—调控证据—发送端与接收端吻合”的顺序组织。','分析依次呈现免疫全景与CD4表型、巨噬细胞通讯、两个独立候选及功能程序。')
 b=b.replace('基因集内容保存在分析对象中，未获得可独立核实的 MSigDB 发布号，因此不据软件号推断数据库发布版本。','基因集内容及成员保存在分析对象中。')
 b=b.replace('历史脚本调用 UCell 默认接口，未显式冻结 maxRank 参数，不能把小鼠的显式参数移植给人类。','人类评分沿用原始对象与既有UCell调用。')
 exec(compile(b,f'methods_{new}','exec'))
