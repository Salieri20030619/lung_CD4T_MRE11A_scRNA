from pathlib import Path
import json,shutil,pandas as pd
H=Path(r'C:\Users\31937');O=H/'lung_selection_joint_revision_20261017_v30';B=H/'lung_citation_revision_20261017_v29';s=(B/'report.py').read_text(encoding='utf-8').replace('矽肺_CD4候选研究_引文修订_v29_20261017','矽肺_CD4候选研究_筛选与联合分析_v30_20261017').replace('矽肺免疫全景_CD4候选与验证_报告_v29.pdf','矽肺免疫全景_CD4候选与验证_报告_v30.pdf')
for n in ['editorial_replacements.json','extra_manifest.json']:shutil.copy2(B/n,O/n)
shutil.copy2(B/'figures/S_C.png',O/'figures/S_C.png')
m=json.loads((B/'figure_manifest.json').read_text(encoding='utf-8'))
for a in m:
 if a['id'] in ['F04','F15','F16']:a['file']=str(O/'figures'/(a['id']+'.png'))
 if a['id']=='F04':a['title']='联合10样本支持Naive减少与Treg增加';a['caption']='上图展示全部10个独立样本的Core CD4内组成；C/S为对照/矽尘，314/183分别为GSE314943/GSE183682。下图比较主队列3+3与分层联合5+5的比例效应；联合模型控制队列-时间层，80种保持层内组别数的分配，3个状态构成BH家族。组间效应以比例单位表示，0.10对应10个百分点。'
 if a['id']=='F15':a['title']='联合分析未支持焦亡相关程序的疾病增强';a['caption']='各状态并列主队列3+3及分层联合5+5的疾病效应、原始P及FDR。原始模型20种组别分配，联合模型80种分层分配；两模型分别在完整52项程序-状态家族内BH校正。效应为评分单位，不是实测死亡率。'
 if a['id']=='F16':a['title']='联合Treg凋亡程序上升，线粒体应激缺少组间支持';a['caption']='模型、独立样本及检验家族同图15。MitoISR为线粒体应激转录程序，ApoExec为凋亡相关程序。Treg联合P=0.050未满足预设P<0.05，不按显示精度改称显著。分层共同效应与原六样本均值差具有不同估计范围。'
(O/'figure_manifest.json').write_text(json.dumps(m,ensure_ascii=False,indent=2),encoding='utf-8')
# Summary and navigation reflect all valid samples, rather than only six-sample descriptive signals.
s=s.replace('初始样CD4（Naive）比例下降、活化样CD4（ActCD4）及调节性T细胞（Treg）比例上升，组成检验P分别为0.10、0.20、0.10。','分层联合5+5样本显示初始样CD4（Naive）比例下降38.52个百分点、调节性T细胞（Treg）上升25.61个百分点，两者P=0.025、FDR=0.0375；活化样CD4（ActCD4）上升12.91个百分点，P=0.150。')
s=s.replace('主队列三类CD4状态的线粒体应激程序均向上（P=0.80、0.80、0.40），ActCD4、Treg焦亡相关评分P为0.30、0.20，均未检出显著组间差异。','联合5+5中ActCD4、Treg焦亡程序P分别为1.000、0.675，线粒体应激P为0.725、0.375；Treg凋亡程序P=0.050、FDR=0.144，未通过预设P<0.05。')
s=s.replace('Naive比例下降，ActCD4/Treg上升；组成检验未达名义显著。Treg候选存在疾病差异及功能程序变化。','联合5+5支持Naive下降、Treg上升（均P=0.025、FDR=0.0375）；ActCD4比例P=0.150。Treg候选另有分子差异。')
s=s.replace('主队列线粒体应激、凋亡与焦亡程序未检出显著组间差异；候选关联提供另一层探索线索。','分层5+5未支持焦亡或线粒体应激程序普遍增强；Treg凋亡P=0.050。候选关联单独提供探索线索。')
s=s.replace('三候选为何推荐｜来源、分层与研究价值','候选选择全程｜补充来源、定向评估与最终推荐')
a=s.index("page();p('候选选择全程｜");b=s.index("page();p('研究设计与数据范围'",a)
new="""page();p('候选选择全程｜补充来源、定向评估与最终推荐','h')
p('候选选择包含集合构建、研究对象提名和深入评价三个步骤。差异/功能/巨噬信号交集提供55个基因；补充清单用来核对尚未被该交集覆盖的研究对象。随后选定9个对象追加表型与损伤关联，最后结合结果和实验可测性形成3个方向。')
tab([['步骤','具体构成','该步解决的问题'],['补充42个','11个候选比较对象 + 27个高效应核查信号 + 4个定向补核对象，三组不重叠','补回待评估线索，检查高效应是否来自目标T细胞；不是42个均已具备损伤功能的基因'],['合并90个','55个三项交集 + 42个补充 − 7个重叠；分为52个计算支持、37个待补、1个未达条件','完整记录资格与证据缺口，避免用后续排名删除其他基因'],['定向深入评估9个','功能机制候选5个：IL1R2、DDIT4、SKIL、FOXP1、GSN；补充探索4个：AHNAK、TTC39C、MDFIC、ST6GALNAC3','这是研究对象提名，在追加关联分析前确定；不是90个经过一个统计门槛后剩9个'],['推荐3个方向','IL1R2：Treg与IL-1环境；TTC39C：存活/损伤；SKIL：Treg状态','依据九候选分析中的支持、主要阴性/敏感性结果与功能文献，确定首轮可检验问题']],[85,240,180])
p('为什么保留补充集合','h2');p('11个对象来自候选比较清单；27个信号来自高效应基因的逐项来源与文献核查表；另定向加入DUSP5、GSDMD、CXCR6、MRE11。补充清单包含免疫球蛋白及组织相关信号，是为了完成来源核查，不能仅凭大效应认定为T细胞机制靶点。完整42基因及来源分组见附录A。')
p('如何判断最后三个','h2');p('IL1R2与SKIL具备联合疾病差异和明确Treg研究切入点；TTC39C由凋亡关联及功能先例支持损伤探索。九候选逐项结果与其他六者的保留理由见后文对照表。推荐表示当前预实验价值，不是所有指标均居前三。')
"""
s=s[:a]+new+s[b:]
# Add sample-scope page before proportions.
pos=s.index('def early_framework():')
fun="""def sample_scope_page():
 page();p('统计范围｜哪些结果使用3+3，哪些使用5+5','h')
 p('现有疾病数据包含主队列6个和补充队列4个独立样本。能够对齐的样本级差异、程序评分和关联采用分层联合10样本；需要已验证细胞来源或特定图谱结构的分析保留其实际队列范围。')
 tab([['分析','实际使用范围','检验与解释'],['疾病差异表达','主队列6及联合10分别报告','edgeR准似然模型；P不受20种标签分配的最小值限制'],['CD4组成','主队列3+3，并补充分层5+5','比例为样本内Core CD4分母；联合共同效应控制队列/时间层，3状态独立校正'],['焦亡/凋亡/线粒体评分','主队列3+3与分层5+5并列','共同程序评分可对齐；联合模型已有完整结果。图15—16同时显示'],['候选程序关联','联合10为M1；M2进一步调整疾病','M1层内2,880种秩排列；M2只有36种有效条件排列，不能只按总n推断分辨率'],['图谱、发送端通讯、TF、轨迹','主要为主队列6样本；补充队列图谱另列','补充队列缺少对应完整发送端/调控与轨迹分析，不能仅合并细胞就称相同验证']],[105,170,230])
 p('3+3的某些P值为什么不可能小于0.05','h2');p('双侧组标签完全枚举共有C(6,3)=20种；观察分组与互换分组效应绝对值相同，最小P=2/20=0.10。分层5+5按三个队列-时间层保留组别数，共20×2×2=80种，最小P=2/80=0.025。该限制适用于这里的组间比例/评分检验，不适用于全部分析。')
 p('因此，3+3的P=0.10不能单独证明没有变化；联合数据仍需看实际效应与方向。联合模型假设可比较的共同疾病效应，跨队列不一致会削弱估计，不能为获得较小P而忽略分层。','small')

"""
s=s[:pos]+fun+s[pos:];s=s.replace("if m['id']=='F03':state_guide()","if m['id']=='F03':state_guide()\n if m['id']=='F04':sample_scope_page()")
pos=s.index('for m in manifest:')
upd={
'F04':'联合5+5的分层比例效应显示Naive下降38.52个百分点、Treg上升25.61个百分点，两者P=0.025、FDR=0.0375；ActCD4增加12.91个百分点，P=0.150、FDR=0.150。该结果支持优先深入Treg状态；3+3的变化方向保留，但其双侧精确P最小为0.10。',
'F06':'Treg中IL1R2联合logFC=3.290、P=1.66×10^-6、FDR=0.000662；SKIL logFC=1.438、P=2.15×10^-5、FDR=0.00314，构成两者的疾病差异依据。TTC39C联合DE P=0.163、FDR=0.803，因此其推荐必须由后续损伤关联与功能先例解释。',
'F09':'Treg核糖体集合负向富集（NES=-3.228），细胞因子-受体集合正向富集（NES=1.973，P=0.000108、FDR=0.00685），其核心成员包括IL1R2。富集把Treg差异与免疫环境联系起来；它属于通路背景，不能单独证明某个候选作用。',
'F12':'IL1A/IL1B与IL1R2的诱饵结合为IL1R2提供明确的环境切入点；CCL2-CCR2及CXCL10-CXCR3为后续SKIL响应先验提供接收端线索。此页解决信号相容性，具体候选的预测和程序关联由后续结果分别评价。',
'F15':'联合5+5中ActCD4焦亡程序效应为+0.00152、P=1.000、FDR=1.000，Treg为+0.01465、P=0.675、FDR=0.836；均未支持疾病增强。相比3+3，联合效应明显减弱，不能把当前结果概括为“增加样本后焦亡会显著”。',
'F16':'联合Treg的ApoExec效应为+0.08148、P=0.050、FDR=0.144，未通过P<0.05；ActCD4为+0.01495、P=0.400、FDR=0.594。联合MitoISR在ActCD4/Treg为+0.00774/+0.02016，P=0.725/0.375，FDR=0.857/0.574。损伤程序没有普遍的疾病组间支持，后续候选关联应作为另一层问题解释。',
'F19':'联合Treg上调目标集以ICAM1和TNF的校正AUPR居前，ActCD4以WNT4和IL10居前。这里排序的是解释响应集合的配体，不是候选基因优劣；SKIL对应的具体配体-靶先验见下一图。',
'F20':'IL1R2的IL-1诱饵结合与SKIL的CXCL10、CCL2、SPP1靶先验提供不同上游切入点。TTC39C未满足该联合差异网络输入条件，因此以程序关联和实验可测性作损伤探索；H2-DM记录仅作为抗原呈递信息。'}
s=s[:pos]+'results.update('+repr(upd)+')\n'+s[pos:]
# Two readable nine-gene comparison pages precede final integrated evidence.
a=s.index("page();p('三候选证据连接｜")
comp="""page();p('九候选深入评价Ⅰ｜五个功能机制候选的结果分化','h')
p('这五个对象从三项交集及功能文献中定向提名，可分别提出免疫调节、应激、凋亡易感性或细胞保护问题。追加分析在提名之后进行，用来检验这些问题，而不是倒推提名门槛。M1/M2含义沿用前文。')
tab([['对象','进入深入评价的研究问题','已有结果如何影响推荐'],['IL1R2 / Treg','Treg免疫环境与候选功能；直接诱饵受体关系','联合DE FDR=0.000662；Inhib M1 rho=0.885、P=0.00903、FDR=0.0744，M2 P=0.278；PyroExec P=0.988。结合Treg研究先例，优先表型/IL-1环境验证，不定为已证实焦亡靶点。'],['SKIL / Treg','状态调控及应激程序','联合DE FDR=0.00314；TAct M1 rho=0.965、P=0.00208、FDR=0.0744，M2 P=0.0278、FDR=0.8681；ApoExec M2 P=0.917。保留状态调控路线，损伤连接仍待检验。'],['DDIT4 / ActCD4','应激信号是否连接目标损伤','联合DE P=0.00337、FDR=0.288；PyroExec M1 rho=0.126、P=0.751、FDR=0.939。疾病差异没有转化为当前焦亡程序支持，保留应激备选。'],['FOXP1 / Core、Treg','状态重塑与凋亡易感性','Core-PyroExec M1 rho=-0.822、P=0.0219、FDR=0.146；等权状态组成敏感性P=0.494。组成影响限制其作为直接损伤候选，保留状态参照。'],['GSN / ActCD4','胞内保护与来源相关信号','ApoExec M1 rho=-0.459、P=0.264、FDR=0.459；M2 P=1.000。当前损伤关联不足，保留备选而不优先投入。']],[90,145,270])
p('五者属于计算支持层，但具有不同证据强弱；资格通过不意味着所有损伤检验阳性。FOXP1两个状态构成两个组合，因此九基因共有十个基因-状态组合。','small')
page();p('九候选深入评价Ⅱ｜四个补充对象与三条推荐路线','h')
p('四个补充对象用于核查T细胞状态或潜在损伤线索，均处于功能证据待补层；它们不因未进入三项交集而被直接丢弃。下表并列后续关联及限制。')
tab([['对象','追加评价的切入问题','结果与当前取舍'],['AHNAK / Naive','T细胞表型与氧化代谢联系','OXPHOS M1 rho=0.811、P=0.0417、FDR=0.244；M2 P=0.0833。提供代谢/状态线索，尚不足以优先定位线粒体损伤。'],['TTC39C / ActCD4','存活、凋亡及线粒体状态','ApoExec M1 rho=0.943、P=0.00347、FDR=0.0556；M2 rho=0.922、P=0.0278、FDR=0.556。另有TLike-MitoISR时间调整线索及肺癌细胞干预先例，推荐损伤探索；联合DE与上游仍不足。'],['MDFIC / ActCD4','补充凋亡相关候选','ApoExec M1 rho=0.947、P=0.00278、FDR=0.0556，强度并不低于TTC39C；M2 P=0.0833、FDR=0.794，外部方向和功能背景限制迁移，保留损伤备选。'],['ST6GALNAC3 / ActCD4','应激与氧化代谢程序','MitoROS M1 rho=0.872、P=0.00972、FDR=0.0972；M2 P=0.222、FDR=0.889。关联对疾病调整敏感，目标损伤直接功能依据不足，保留探索。']],[95,135,275])
p('由结果形成的推荐','h2');p('IL1R2以联合疾病差异、直接IL-1结合及Treg功能背景支持表型路线；SKIL以联合差异和较稳定的活化程序关联支持状态路线；TTC39C在补充对象中提供凋亡、时间应激及可参考功能实验三类线索，承担损伤探索。TTC39C不是因原始相关最强而胜出，MDFIC仍作为备选保留。')
p('三条路线用于验证同一巨噬—T细胞主线的不同环节。选择不等同全部90基因的统一总分前三，也不把名义关联写成已验证机制。完整九候选效应、P/FDR、模型和检验家族见配套表。','small')
"""
s=s[:a]+comp+s[a:]
# Replace redundant remaining-candidate paragraph with synthesis tied to actual results.
a=s.index("p('其余候选继续按证据分层保留','h2')");b=s.index("page();p('功能文献",a)
s=s[:a]+"p('综合结论：先验证表型，再连接损伤','h2');p('联合样本支持Naive/Treg组成变化，IL1R2与SKIL有疾病表达差异；候选程序关联将IL1R2、SKIL、TTC39C分别连接到免疫调节、活化和凋亡问题。联合组间评分尚不支持普遍焦亡或线粒体应激增强，因此推荐以目标T细胞表型和独立损伤读数验证联系。其余六个深入评价对象按上述结果保留，不由一个阴性检验判定无生物学作用。')\n"+s[b:]
# Sampling boundary discussion rather than treating unobtainable significance as a biological negative.
s=s.replace('精确枚举保证既定排列空间计算完整，不证明可交换性或模型假设成立。MDFIC等强效应而P有限的结果不能全部归因于分辨率。','3+3双侧组间精确检验最小P=0.10，不能用未过0.05否定明确效应。分层5+5最小P=0.025，但共同疾病效应依赖跨队列可比性；程序效应减弱反映结果差异，不能仅用样本量解释。精确枚举不证明可交换性假设成立。')
# Detailed provenance with the complete42 list in a customer-readable appendix.
pos=s.index("screening_methods()\nfor m in manifest:")
src=pd.read_csv(O/'tables/supplement42_provenance.tsv',sep='\t');groups={g:'、'.join(z.gene.tolist()) for g,z in src.groupby('origin_group')}
add="\npage();p('附录A｜42个补充对象的来源与名单','h')\np('补充集合由候选比较清单、高效应信号核查表及定向补核名单构成。各组是研究对象来源，不是新的独立统计筛选；统一进入90基因评价后再记录资格和去向。')\n"
for g in ['候选比较对象11','高效应核查信号27','定向补核对象4']:add+="p("+repr(g)+",'h2');p("+repr(groups[g])+")\n"
add+="p('11个与27个清单无重叠；4个定向对象均为额外加入，合计42。与55个三项交集重叠7个，故总评价集合为90。原始成员、来源表对应及后续资格逐基因保存在supplement42_provenance.tsv。','small')\n"
s=s[:pos]+s[pos:].replace('screening_methods()\nfor m in manifest:','screening_methods()'+add+'for m in manifest:',1)
# Explicit model description for new proportion sensitivity and existing score display.
needle="p('GO/KEGG与GSEA','h2')"
s=s.replace(needle,"p('联合细胞比例补充分析使用两队列相同Core CD4三状态定义，以每样本各状态细胞数/Core细胞数为结局，线性模型控制队列-时间层并估计共同疾病效应；保持每层组别数枚举80种分配。Naive、ActCD4、Treg共3项为一个BH家族，主队列3项另行校正。程序评分继续使用既定完整52项家族。','small')\n"+needle)
(O/'report.py').write_text(s,encoding='utf-8');print('REPORT_PREPARED')