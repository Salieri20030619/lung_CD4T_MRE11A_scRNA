from pathlib import Path
import json,shutil,ast,hashlib
H=Path(r'C:\Users\31937');B=H/'lung_readability_revision_20261017_v28';O=H/'lung_citation_revision_20261017_v29';O.mkdir(exist_ok=True)
for folder in ['figures','evidence']:(O/folder).mkdir(exist_ok=True)
for n in ['figure_manifest.json','editorial_replacements.json','extra_manifest.json']:shutil.copy2(B/n,O/n)
shutil.copy2(B/'figures/S_C.png',O/'figures/S_C.png')
s=(B/'report.py').read_text(encoding='utf-8').replace('矽肺_CD4候选研究_可读性修订_v28_20261017','矽肺_CD4候选研究_引文修订_v29_20261017').replace('矽肺免疫全景_CD4候选与验证_报告_v28.pdf','矽肺免疫全景_CD4候选与验证_报告_v29.pdf')
a=s.index("page();p('功能文献为实验提供依据与方向约束'");b=s.index("page();p('实验验证建议｜",a)
new="""page();p('功能文献为实验提供依据与方向约束','h')
p('IL1R2：以Treg及IL-1环境响应为切入点','h2');p('Treg特异Il1r2缺失的肿瘤研究报告克隆扩增和效应特征改变，相关单细胞/TCR分析包含MC38肿瘤及抗PD-1处理背景。[1] Tfr研究支持其调节生发中心的IL-1环境，不能直接等同全部成熟Treg。[2]')
p('文献[3]报告了巨噬细胞/脓毒症背景下的IL1R2—ENO1、糖酵解与焦亡联系。矽尘暴露小鼠中性粒细胞中已有Il1r2转录表达报道，属于表达证据。[4] 本项目研究问题仍是矽肺相关巨噬环境对CD4/Treg功能的影响。')
p('文献[3]注：原文Results与Discussion关于ENO1遗传降低效应的方向表述不一致。本报告保留其研究背景，不据此预设ENO1干预方向，也不替原作者判定哪处表述有误。','small')
p('TTC39C：按实测判断存活/凋亡作用方向','h2');p('人A549与NCI-H1299细胞经shRNA降低TTC39C后，Annexin V-APC单染阳性增加。[5] 这提供功能先例，不能将本项目TTC39C与ApoExec的相关性直接解释为CD4中的促凋亡作用。')
p('SKIL：区分蛋白、分化背景与状态效应','h2');p('SKIL编码SnoN，SKI为另一基因。Arkadia缺失的诱导Treg体系中，Ski/SnoN双靶向能救援分化，单靶向不足，提示该背景下的冗余。[6] PM2.5暴露及肺上皮细胞模型中已有ROS/SnoN与肺纤维化联系，不能直接外推为矽尘Treg机制。[7]')
p('检测资源与适用范围','h2');p('Tfr原研究使用4E2进行IL1R2流式检测。[2] 厂家另有小鼠CD121b流式及TTC39C、SnoN的WB资料。厂家、货号、物种、验证材料及官方链接见附录E；检测抗体的适用性与本项目原代T细胞功能干预效力分别核对。')
"""
s=s[:a]+new+s[b:]
s=s.replace('H2-DM相关数据库记录另列抗原呈递网络信息，不作为普通表面配体-CD4作用轴推荐。','H2-DM相关记录另列抗原呈递网络信息；其中H2-DMb1的溶酶体膜定位与抗原肽装载背景见Reactome条目[10]。H2-DMb2保留为项目数据库先验，不以该单一H2-DMb1条目证明其直接作用。')
m=json.loads((O/'figure_manifest.json').read_text(encoding='utf-8'))
for x in m:
 if x['id']=='F20':x['caption']+=' H2-DMb1细胞定位的数据库背景见[10]，该条目不作为H2-DMb2直接功能证据。'
(O/'figure_manifest.json').write_text(json.dumps(m,ensure_ascii=False,indent=2),encoding='utf-8')
s=s.replace('LDH、IL-1β可来自巨噬细胞或残留条件培养基，不能据此定位T细胞死亡。[11]','矽尘刺激巨噬细胞的研究检测到LDH、IL-1β变化。[11] 因此本方案需排除巨噬来源及培养基残留，混合上清不能单独归属T细胞死亡。')
s=s.replace('原始功能研究涉及肿瘤Treg—IL-1—成纤维细胞环境及Tfr免疫调节；ENO1/焦亡机制来自巨噬细胞。[1–3] 本项目优先检验这些作用能否迁移到矽尘相关巨噬—Treg环境。','肿瘤Treg与成纤维细胞环境研究[1]、Tfr的IL-1调节研究[2]提供不同功能背景。文献[3]报告巨噬细胞中的ENO1/焦亡联系，其ENO1干预方向文字不一致见前述文献注记。本项目检验这些背景能否迁移到矽尘相关巨噬—Treg环境。')
s=s.replace('肺癌细胞研究中，降低TTC39C后凋亡增加，为存活/凋亡方向提供先例。[5]','肺癌细胞研究中，shRNA降低TTC39C后Annexin V-APC阳性增加，为存活/凋亡方向提供先例。[5]')
s=s.replace('模型阶段核对目标T细胞相关gasdermin蛋白、裂解性死亡读数和检测条件。按表达与结果选择GSDMD/GSDME等执行轴，不默认采用巨噬NLRP3—caspase-1—GSDMD流程。[12–13] 焦亡机制结合目标T细胞执行层、膜损伤/死亡以及后续执行轴干预；Annexin V、LDH、IL-1β或caspase-3单项均不足以独立判定。[8]','Gasdermin成孔[12]与特定条件下caspase-3切割GSDME介导焦亡[13]已有原始研究。本方案在模型阶段核对目标T细胞相关蛋白、裂解性死亡读数和检测条件，按实际表达与结果选择执行轴。依据机制区分死亡类型的原则[8]，本方案组合执行层、膜损伤/死亡及后续执行轴干预；不以Annexin V、LDH、IL-1β或caspase-3单项完成机制判定。')
# Unified bibliography: seven journal items on first page, remaining sources on second.
refs=[
('[1]','Chen L, et al. IL1R2 increases regulatory T cell population in the tumor microenvironment by enhancing MHC-II expression on cancer-associated fibroblasts. Journal for ImmunoTherapy of Cancer. 2022;10(9):e004585.','https://doi.org/10.1136/jitc-2022-004585'),
('[2]','Pyrillou K, et al. Loss of T follicular regulatory cell-derived IL-1R2 augments germinal center reactions via increased IL-1. JCI Insight. 2024;9(5):e174005.','https://doi.org/10.1172/jci.insight.174005'),
('[3]','Tan C, et al. Critical Role of IL1R2-ENO1 Interaction in Inhibiting Glycolysis-Mediated Pyroptosis for Protection Against Lethal Sepsis. Advanced Science. 2025;12(39):e02297.','https://doi.org/10.1002/advs.202502297'),
('[4]','Lam M, et al. Epithelial NLRP3 drives silica-induced lung injury and fibrosis through IL-18 and pro-fibrotic neutrophil recruitment. Particle and Fibre Toxicology. 2026;23:30.','https://doi.org/10.1186/s12989-026-00682-9'),
('[5]','Rong H, et al. Ttc39c is a potential target for the treatment of lung cancer. BMC Pulmonary Medicine. 2022;22:391.','https://doi.org/10.1186/s12890-022-02173-x'),
('[6]','Xu H, et al. Arkadia-SKI/SnoN signaling differentially regulates TGF-β-induced iTreg and Th17 cell differentiation. Journal of Experimental Medicine. 2021;218(11):e20210777.','https://doi.org/10.1084/jem.20210777'),
('[7]','Li X, et al. PM2.5 regulates TGF-β1/Smads-mediated pulmonary fibrosis via ROS/SnoN in vitro and in vivo. International Immunopharmacology. 2025;153:114477.','https://doi.org/10.1016/j.intimp.2025.114477'),
('[8]','Galluzzi L, et al. Molecular mechanisms of cell death: recommendations of the Nomenclature Committee on Cell Death 2018. Cell Death & Differentiation. 2018;25(3):486-541.','https://doi.org/10.1038/s41418-017-0012-4'),
('[9]','Thermo Fisher Scientific. Tetramethylrhodamine, Methyl Ester, Perchlorate (TMRM) - FAQs. Catalog number T668. 网页；访问日期：2026-10-17.','https://www.thermofisher.com/order/catalog/product/T668/faqs'),
('[10]','Reactome. H2-DMb1 [lysosomal membrane]. Stable identifier: R-MMU-2130673. Mus musculus. 数据库条目；访问日期：2026-10-17.','https://www.reactome.org/content/detail/R-MMU-2130673'),
('[11]','Yin H, et al. Acute Silica Exposure Triggers Pulmonary Inflammation Through Macrophage Pyroptosis: An Experimental Simulation. Frontiers in Immunology. 2022;13:874459.','https://doi.org/10.3389/fimmu.2022.874459'),
('[12]','Ding J, et al. Pore-forming activity and structural autoinhibition of the gasdermin family. Nature. 2016;535(7610):111-116.','https://doi.org/10.1038/nature18590'),
('[13]','Wang Y, et al. Chemotherapy drugs induce pyroptosis through caspase-3 cleavage of a gasdermin. Nature. 2017;547(7661):99-103.','https://doi.org/10.1038/nature22393')]
a=s.index('refs=[');b=s.index("p('补充材料','h2')",a)
s=s[:a]+'refs='+repr(refs)+"\nfor i,(n,title,url) in enumerate(refs):\n if i==7:page();p('参考文献（续）','h')\n p(n+' '+title+'\\n'+url,'small')\n"+s[b:]
# Make reagent evidence readable without a separate attachment.
pos=s.index('# Layout uses')
append="""page();p('附录E｜候选蛋白检测资源与适用范围','h')
p('以下为检测资源核对，不是本项目原代T细胞中的已验证操作方案。应用和物种依据厂家公开资料，验证材料单列；正式选用时确认具体制剂、批次和用途。核对日期：2026年10月17日。')
tab([['靶标','厂家与货号','应用/物种','厂家验证材料','本项目解释'],['IL1R2','BD；4E2；554450（PE）','流式；小鼠（QC）','厂家常规测试流式；原Tfr研究另有4E2检测先例[2]','用于检测资源参考。厂家描述克隆可阻断IL-1结合，但该PE制剂不能自动视为本项目活细胞功能阻断方案。[R1]'],['TTC39C','Santa Cruz；G-6；sc-390575','WB、IP、IF、ELISA；小鼠/大鼠/人','小鼠TTC39C转染293T与未转染对照；另列细胞裂解物WB','尚非原代CD4验证；需确认目标条带、靶向效率和阴性对照。[R2]'],['SnoN / SKIL','Santa Cruz；B-3；sc-136958','WB等；小鼠/大鼠/人','小鼠SnoN转染293T与未转染对照','厂家声明识别SnoN和SnoN2；区分条带及亚型，原代Treg应用需验证。[R3]']],[65,110,100,115,115])
p('用途与机制证据分别管理','h2');p('流式检测、受体阻断、细胞清除和基因干预不是同一种用途。文献[6]提供SnoN/SKI抗体及蛋白研究背景，当前替代货号的适用性须按自身厂家资料核对。蛋白可检出不等于候选具有预期功能。')
p('官方来源','h2')
p('[R1] BD Biosciences. PE Rat Anti-Mouse CD121b, 554450. 网页；访问日期：2026-10-17.\\nhttps://www.bdbiosciences.com/en-gb/products/reagents/flow-cytometry-reagents/research-reagents/single-color-antibodies-ruo/pe-rat-anti-mouse-cd121b.554450','small')
p('[R2] Santa Cruz Biotechnology. TTC39C (G-6): sc-390575. 产品说明书；访问日期：2026-10-17.\\nhttps://datasheets.scbt.com/sc-390575.pdf','small')
p('[R3] Santa Cruz Biotechnology. SnoN Antibody (B-3): sc-136958. 网页；访问日期：2026-10-17.\\nhttps://www.scbt.com/p/snon-antibody-b-3','small')
"""
s=s[:pos]+append+s[pos:]
(O/'report.py').write_text(s,encoding='utf-8');(O/'references.json').write_text(json.dumps(refs,ensure_ascii=False,indent=2),encoding='utf-8')
resources=[{'target':'IL1R2','manufacturer':'BD Biosciences','clone':'4E2','catalog':'554450','application':'Flow cytometry (Routinely Tested)','reactivity':'Mouse (QC Testing)','material':'Manufacturer flow-testing data; separate original Tfr clone precedent','scope':'PE detection reagent; clone-level blocking claim not equivalent to validated live-cell functional protocol','url':'https://www.bdbiosciences.com/en-gb/products/reagents/flow-cytometry-reagents/research-reagents/single-color-antibodies-ruo/pe-rat-anti-mouse-cd121b.554450'}, {'target':'TTC39C','manufacturer':'Santa Cruz Biotechnology','clone':'G-6','catalog':'sc-390575','application':'WB; IP; IF; ELISA','reactivity':'Mouse; rat; human','material':'Mouse TTC39C-transfected versus non-transfected 293T; cell lysates','scope':'Not primary CD4 validation','url':'https://datasheets.scbt.com/sc-390575.pdf'}, {'target':'SnoN','manufacturer':'Santa Cruz Biotechnology','clone':'B-3','catalog':'sc-136958','application':'WB; IP; IF; IHC(P); ELISA','reactivity':'Mouse; rat; human','material':'Mouse SnoN-transfected versus non-transfected 293T','scope':'Recognizes SnoN and SnoN2; not primary Treg validation','url':'https://www.scbt.com/p/snon-antibody-b-3'}]
(O/'evidence/reagent_resources.json').write_text(json.dumps(resources,ensure_ascii=False,indent=2),encoding='utf-8')
(O/'change_scope.json').write_text(json.dumps({'base':'v28','ref3_action':'annotation only; no attempted resolution of source discrepancy','statistics_changed':False,'candidate_order_changed':False,'base_sha256':hashlib.sha256(Path(json.loads((B/'report_delivery.json').read_text(encoding='utf-8'))['pdf']).read_bytes()).hexdigest()},indent=2),encoding='utf-8');print('PREPARED')