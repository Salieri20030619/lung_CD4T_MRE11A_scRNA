h('7  Bootstrap：供者抽样稳定性及调整的代价',1)
h('7.1  固定问题、输入与重抽样流程')
p('这部分完整恢复20261004已完成的Bootstrap。它回答“现有供者名单变化时，相关强度、符号及调整差异有多不确定”，没有新增候选、阈值或终点。旧探索限定两组配对：MacAM IL18与IL7R IKZF2；MacAM TNF与Core STAT1。每组分别检查配体–TF、TF–StateEffector、配体–StateEffector三个相关。三段符号吻合不是中介检验，也不要求独立候选必须沿此链成立。')
p('StateEffector是前半段冻结的效应状态分数：去除TF调控子重叠成员后保留PRF1、IFNG、NKG7；先在冻结Core细胞集合内按基因标准化log1p CP10k，再求细胞及供者均值。它不是从疾病DE选择的靶集。前半段30个疾病比较没有名义支持，因此下述关联本身不赋予该分数“疾病表型”身份。')
table('Bootstrap实际规格',['项目','冻结设置与含义'],[
['名单 / 读数','RNA19；与AUC匹配的RNA17；AUC seed17/29/43。RNA为供者平均CP10k，AUC为既有DropOverlap供者平均活性。RNA19与匹配17是不同问题。'],
['独立单位','整行供者连同发送端、接收端、来源/疾病、深度及各seed共同抽样；不是抽细胞，不增加独立样本量。'],
['主抽样','按来源×疾病分层有放回，每层维持原n。RNA19层大小2/8/2/7，匹配17为1/7/2/7；17人中唯一G135对照在每次分层抽样中固定出现。'],
['5000次与配对','每种名单/方案5000次。RNA17、三AUC seed、两条探索及不同调整共享抽样索引，以计算成对Δρ；RNA19使用另一组索引。'],
['随机性','NumPy PCG64：主分层19/17种子2026100401/02；未分层敏感性03/04。未重新训练SCENIC或重算AUCell。'],
['调整与计算','Unadjusted为Spearman；B控制G、S及两端深度的秩；B_Treg仅TNF加入Treg比例的秩。每次抽样重新求秩并标准化连续协变量后残差相关。'],
['区间 / 失败','有效抽样的百分位95%CI（线性分位插值）；有效数<1000不出区间。秩亏、df<3、常量或残差近零为失败，不补抽。报告5000次尝试中有效数。'],
['性质','条件于既有来源与疾病构成的抽样不确定性。借鉴多层数据保留独立单位的原则[12]；本项目实际是供者分层配对Bootstrap，并非另作细胞层级抽样。']],[3.6,13.8])
h('7.2  已有结果：同时保留不跨零与跨零区间')
def boot_table(axis,adj):
 q=df('bootstrap_effect_summary');q=q[(q.scheme=='stratified')&(q.axis==axis)&(q.adjustment==adj)];rows=[]
 for modality in ['RNA19','RNA17','AUC17','AUC29','AUC43']:
  v=q[q.modality==modality];row=[modality]
  for edge in ['ligand_TF','TF_phenotype','ligand_phenotype']:
   r=v[v.edge==edge].iloc[0];row.append(f'{signed(r.original_effect,3)} [{signed(r.CI95_low,3)}, {signed(r.CI95_high,3)}]')
  row.append('/'.join(str(int(x)) for x in v.valid.unique())+'/5000');rows.append(row)
 table(axis+'：分层Bootstrap '+adj,['读数','配体–TF：ρ [CI]','TF–Effector：ρ [CI]','配体–Effector：ρ [CI]','有效/尝试'],rows,[2,4.2,4.2,4.2,2.8],note='CI为既有5000次抽样的百分位区间；ρ为原名单效应。RNA19 n19，RNA17及AUC n17；AUC17/29/43指seed，不是新增名单。原方法不产生P/FDR，故不列这两栏。[A5]')
boot_table('IL18','B')
fig('R03_bootstrap_IL18.png','IL18–IKZF2的区间与读数比较','主分层/B模型；点=原名单ρ，线=95%Bootstrap CI，黑短线=抽样中位数，虚线=0；颜色区分RNA19、RNA17和AUC三种子。L–TF=MacAM IL18与IL7R IKZF2；TF–E=IKZF2与StateEffector；L–E=IL18与StateEffector。各行5000/5000有效。所有L–TF区间跨零，TF–E调整后区间为负。[A5]')
p('IL18的原配体–IKZF2效应为正，但B模型中所有读数的95%CI均跨零。TF–Effector为负且上述区间不跨零；IL18–Effector仍跨零。因此支持的是特定TF–程序关联的抽样方向，不能写成IL18已通过该TF抑制效应程序。AUC seed29−17的配体–TF差为−0.219，CI[−0.617,+0.030]；seed43−17为−0.191，CI[−0.558,+0.074]。种子间点估计不同，成对区间未确认差异。')
boot_table('TNF','B');boot_table('TNF','B_Treg')
fig('R03_bootstrap_TNF.png','TNF–STAT1的B模型与Treg调整','点=原ρ、线=95%Bootstrap CI、黑短线=抽样中位数、虚线=0。B控制疾病、来源及两端深度；B_Treg再控制Treg比例。L–TF=MacAM TNF与Core STAT1；TF–E=STAT1与StateEffector；L–E=TNF与StateEffector。RNA19 n19；匹配RNA与AUC n17。B每行5000有效；B_Treg RNA19为5000，匹配17为4998。颜色区分读数。[A5]')
p('TNF的B模型中配体–Effector区间保持正向，配体–TF及TF–Effector区间跨零。加入Treg比例后，配体–TF各读数区间均为正；RNA17和AUC seed17的TF–Effector区间也为正，但其余读数仍跨零。这是同一名单上条件相关的变化，不能直接称为“Treg介导”或“排除Treg后机制确证”。所有15个B_Treg−B成对Δρ区间均跨零，说明模型点估计增强尚不能等同于调整效应已明确。')
h('7.3  稳定性诊断及应如何读“符号兼容率”')
q=df('bootstrap_all_seeds_chain_summary');q=q[(q.scheme=='stratified')&(((q.axis=='IL18')&(q.adjustment=='B'))|((q.axis=='TNF')&q.adjustment.isin(['B','B_Treg'])))]
table('三种子同次抽样的符号兼容统计',['探索 / 调整','兼容 / 有效','比例','有效 / 尝试'],[[r.axis+' / '+r.adjustment,f'{int(r.all_seeds_compatible_count)}/{int(r.valid)}',f'{r.all_seeds_compatible_fraction:.2%}',f'{int(r.valid)}/{int(r.attempted)}'] for r in q.itertuples()],[4.2,4.4,3.4,5.4],note='“兼容”是既定三段符号的代数条件；不是该机制为真的概率、不是P/FDR、不是三次独立重复。三种子来自同一数据。[A5]')
p('分层抽样的独特供者数中位数仅RNA19为13、匹配17为12；重复行使抽样有效信息少于表面行数。匹配17的B_Treg有16.44%尝试在按独特供者诊断时支持自由度不足3，此诊断没有被用作删除区间的二次门槛。')
p('未分层敏感性会抽掉稀少来源/疾病组合：RNA19的B有效3793/5000（75.86%）、B_Treg3792/5000；匹配17分别2807/5000（56.14%）和2805/5000。因此其CI仅条件于可估计抽样，不能用高失败率结果替代主结果或宣称更稳定。全部未调整、B/B_Treg、两抽样方案及成对差表均保留在服务器。[A5]')
p('Bootstrap对主线的贡献是量化来源受限的小样本不确定性，保留真实区间与失败原因。它没有修复来源与疾病混杂，也没有把未通过疾病比较的终点变为疾病程序。')

h('8  NicheNet：从巨噬表达背景到患者疾病靶集',1)
h('8.1  四类证据必须分开')
table('通讯证据层级',['层级','输入与输出','可以 / 不可以说明什么'],[
['LR表达兼容','发送细胞配体与接收细胞受体/必要复合物成员达到覆盖门槛。','支持具备表达条件；低检出记覆盖不足，不判生物学不存在。'],
['配体活性排序','NicheNet先验的ligand–target权重与疾病靶集、表达背景比较。','预测哪个配体较能解释目标集合；排名/AUPR没有等价P/FDR。'],
['靶程序吻合','高权重预测靶确实落入独立单位DE定义的程序，查看实际基因。','支持具体程序的一致性；先验权重不是该患者已验证的受体后信号。'],
['同供者关联','两端足量细胞的供者，配体RNA与接收者程序进行偏秩相关。','关联的P和BH单独报告；同数据选靶并非独立验证，也不能证明因果。']],[3.2,7.1,7.1])
p('NicheNet方法通过整合的配体–靶先验将细胞间信号与接收端基因程序联系起来[10]。本项目先定义疾病DE与表达背景，再对所有符合条件的巨噬来源/接收亚群计算。TNF、IL18用于与既有探索对照，未限制全局候选只有二者。T/NK仅为发送来源竞争敏感性，不能替代巨噬发送端。')
h('8.2  资源、基因映射及有限规格')
table('人类NicheNet实际参数',['项目','已冻结设置'],[
['软件 / 人prior','nichenetr 2.2.1.1、edgeR 4.8.2；经典human v1，Zenodo 3260758。使用正式已安装资源；人类基因符号精确匹配，缺失/重复映射记录。'],
['覆盖范围','7个接收范围×All/BG1×Base/Source×FDR05/Nominal05×Up/Down=112个规格状态。发送端MacAM、MonoMac、MacWide；MacWide是union，非第三个独立复制。'],
['细胞 / 表达门槛','接收者每供者≥20细胞；配体和受体在相应端≥10%细胞、且至少2名CWP供者检出。5%仅作为既有LR覆盖敏感性，不新搜阈值或重排全网。'],
['疾病DE','供者pseudobulk edgeR robust QL。Base counts~G；Source counts~G+S；原全库总量×TMM为offset，不再加平均深度变量。'],
['靶集','FDR<0.05或P<0.05，且|log2FC|≥0.5；Up、Down独立。FDR是该DE全受检基因BH，名义靶集不改称FDR靶集。'],
['背景','受检基因∩在任一疾病组至少2位合格供者中≥10%表达∩prior可映射基因。不是全基因组，也不是与靶集一样的疾病筛选。'],
['可计算性','背景≥500、靶≥10、背景负类≥10、潜在配体≥2。否则保留TARGETSET_LT10/DE不足等状态，不赋零分。'],
['排序 / 来源','合格发送来源的潜在配体union用于同一接收规格排名，再分别标注来源覆盖；不是每个来源另造独立排名。688是资源中的配体总数，不是每个接收端可计算数。'],
['prior审计','LT SHA256：4bcbd2ab9f357905…；LR：5e46ad6b932b7b03…；完整hash、文件版本和成员表在服务器provenance及原资源清单。']],[3.7,13.7])
p('下调靶集也用同一个不带调控方向的先验。某配体对Down集排名高不等于该配体在生物学上抑制这些基因。corrected AUPR校正背景基准；若全列表均为负，排名第一也不能描述为强正活性。')
q=df('stage2_targetset_status');q=q[(q['filter']=='All')&(q.model=='Source')];rows=[]
for receiver in ['Primary','Core','CD4-CCR7','CD4-GZMA','CD4-IL7R','Treg','CD4state']:
 v=q[q.receiver==receiver]
 if len(v)==0:continue
 row=[receiver];first=v.iloc[0];row.append(str(int(first.background_n)) if pd.notna(first.background_n) else 'DE不足')
 for test in ['FDR05','Nominal05']:
  parts=[]
  for direction in ['Up','Down']:
   r=v[(v.test==test)&(v.direction==direction)].iloc[0];parts.append(str(int(r.geneset_n)) if pd.notna(r.geneset_n) else '—')
  row.append('/'.join(parts))
 row.append(str(int(first.candidate_n)) if pd.notna(first.candidate_n) else '—');rows.append(row)
table('人类Source/All：实际背景与靶数',['接收范围','背景基因','FDR Up/Down','名义 Up/Down','候选配体'],rows,[3.5,2.8,3.8,3.8,3.5],note='每个目标集<10即不计算该NicheNet排名。CD4state仅1位Cells20对照，DE不足；其细胞存在不等于该比较可估计。Primary与Core有嵌套，不能当独立验证。[A6]')
h('8.3  人类实际配体活性：有排名，也有不支持')
q=df('human_NN_main');q=q[(q.model=='Source')&(q.direction=='Up')];rows=[]
for receiver,test,ligands in [('Core','Nominal05',['IL15','SPP1','TNF','IL18']),('CD4-IL7R','Nominal05',['IL15','IL18','TNF']),('CD4-GZMA','FDR05',['CLSTN1','ADAM15','ICAM2','IL15','TNF','IL18']),('CD4-GZMA','Nominal05',['IL15']),('Treg','Nominal05',['AIMP1','IL15'])]:
 for ligand in ligands:
  v=q[(q.receiver==receiver)&(q.test==test)&(q.test_ligand==ligand)]
  if len(v):
   r=v.iloc[0];rows.append([receiver,test,ligand,int(r['rank']),num(r.aupr_corrected,6),int(r.geneset_n)])
table('人类All/Source/Up关键结果（同表保留负分）',['接收范围','靶集','配体','排名','corrected AUPR','靶数'],rows,[3.1,2.7,2.6,1.5,4,3.5],note='排名为该规格潜在配体内部排序；不是显著性。FDR05与Nominal05不能合并选择最佳结果。本表从完整候选表按既有重点及排名前列摘录，未重算排名。[A6]')
fig('R04_human_NicheNet.png','人类疾病靶集与候选配体活性','All、Source、Up；色值为corrected AUPR，正/负颜色按图例；文字标示原排名，NA表示该配体在此规格没有活动分数，不代表生物学不存在。各列均为名义上调靶集，GZMA严格FDR靶集在上表另列。有效接收供者Primary/Core n17、CCR7 n9、GZMA n12、IL7R n16、Treg n10，发送支持按MacAM、MonoMac及其union另列。该图没有P/FDR含义。[A6]')
p('Core和IL7R名义上调靶集均以IL15排名1，但它们对应的FDR上调靶数分别3和0，不能称FDR靶集支持。IL7R的IL18虽排名5，corrected AUPR为−0.001073；TNF为−0.006103。GZMA FDR上调10基因可计算，前列为CLSTN1/ADAM15等，IL15/TNF/IL18没有在这一严格靶集居前。GZMA名义集IL15虽排名1但分数−0.004430，Treg名义集第一名AIMP1也为负；这类结果保留为相对排序，不包装成正预测支持。')
q=df('human_LT_examples');rows=[]
for receiver,ligand in [('Core','IL15'),('CD4-IL7R','IL15'),('CD4-GZMA','IL15'),('CD4-GZMA','ICAM2')]:
 v=q[(q.receiver==receiver)&(q.ligand==ligand)].sort_values('weight',ascending=False);rows.append([receiver+'/'+ligand,'；'.join(f'{r.target} ({r.weight:.4f})' for r in v.itertuples())])
table('实际疾病靶集内的预测靶例子',['接收范围 / 配体','先验权重靠前的实际靶成员'],rows,[4.5,12.9],note='均为All/Source/Nominal05/Up，数字是人v1 ligand–target权重；本表证明基因在疾病靶集与先验交集中，不证明直接结合、受体激活或因果。完整靶、背景及ligand–target成员留在服务器。[A6]')
h('8.4  同供者关联：9项BH结果及其脆弱性')
p('配对关联是另一个统计层。接收DE标记Base/Source决定所用疾病程序，但所有下表关联均控制G、S和两端深度秩；两端各≥20细胞，每组至少2位供者、残差df≥3。BH在接收范围/filter/DE模型/靶集/方向内跨配体和发送来源。完整表保留负结果，并非只输出下列9项。')
q=df('human_NN_BH');rows=[]
for r in q.itertuples():rows.append([r.receiver.replace('CD4-','')+'/'+getattr(r,'filter'),r.model+'/'+r.test.replace('Nominal05','Nom').replace('FDR05','FDR')+'/'+r.direction,r.sender+'/'+r.ligand,f'{r.n}({r.n_case}/{r.n_control});{r.df}',signed(r.rho,3),num(r.P,6),num(r.FDR,6)])
table('既有BH<0.05配对关联的全部9行',['接收/过滤','DE/靶/方向','发送/配体','n(病/对);df','ρ','P','BH FDR'],rows,[2.3,3.3,3.2,2.4,1.6,2.2,2.4],note='All/BG1重叠病例；MacWide为union。关联负号为所用程序分数与配体负相关，不自动指向保护作用。全部为同数据选择后的探索性关联。[A6,A7]')
q=df('NN_BH_influence');rows=[]
for r in q.itertuples():rows.append([r.receiver.replace('CD4-','')+'/'+getattr(r,'filter')+'/'+r.ligand,f'{r.computed}/{r.deletions}',int(r.nonidentifiable),f'[{r.rho_min:.3f}, {r.rho_max:.3f}]',int(r.sign_flips),num(r.max_P,5)])
table('上述关联逐供者删除诊断',['接收/过滤/配体','可计算/删除','不可识别','可计算ρ范围','变号数','最大P'],rows,[5,2.7,2.4,3,1.8,2.5],note='仅展示既有leave-one-donor-out诊断；无重新筛选阈值。不可识别不能当作通过稳定性检验。未为删除模型另造FDR。[A7]')
p('CCR7的DSC2、CXCL2及GZMA-BG1的VTN仅9次删除中1次可计算，原ρ接近±1不能克服df=3的设计脆弱性。CXCL2对应Source/FDR/Up靶集仅7个，不满足NicheNet靶数门槛；其配对关联显著不等于NicheNet活性通过。ICAM2、VEGFB、GSTP1各自的删除结果均可计算且无变号；RTN4有1次不可识别。它们仍是小样本同数据关联，不因删除稳定而升级为机制。')
h('8.5  既有小鼠NicheNet与人类对照')
p('小鼠使用正式mouse v2先验（Zenodo 7074291），按小鼠符号计算；跨物种对照通过已核验同源符号映射，不能仅把首字母大写当作已完成映射。人v1与鼠v2先验不同，因此只比较候选与程序方向，不直接比较排名的数值差或AUPR大小。')
q=df('M3_geneset_status');rows=[]
for (receiver,arm),v in q.groupby(['receiver','filter'],sort=False):
 row=[receiver+'/'+arm,int(v.background.iloc[0])]
 for test in ['FDR05','Nominal05']:row.append('/'.join(str(int(v[(v.test==test)&(v.direction==direction)].geneset.iloc[0])) for direction in ['Up','Down']))
 row.append(int(v.ligands.iloc[0]));rows.append(row)
table('鼠M3：Core与ActCD4的靶/背景规模',['接收/过滤','背景','FDR Up/Down','名义 Up/Down','潜在配体'],rows,[4,3,4,4,2.4],note='GSE314943，3+3文库。FDR、名义两套靶集分别计算；ActCD4的FDR靶过少，不能由名义计算完成来补足其证据。旧M5进一步保留Pct10/Pct5和Full/Drop7的有限敏感性。[A6]')
q=df('mouse_NN');q=q[(q['filter']=='All')&(q.gate=='Pct10')&(q.direction=='Up')&(q.omit=='Full')];rows=[]
for receiver,test,ligands in [('Core','FDR05',['Cd86','Cd274','Il15','Tnf','Il18']),('ActCD4','Nominal05',['Il15','Tnf','Il18'])]:
 for ligand in ligands:
  v=q[(q.receiver==receiver)&(q.test==test)&(q.test_ligand.str.upper()==ligand.upper())]
  if len(v):
   r=v.iloc[0];rows.append([receiver+'/'+test,r.test_ligand,int(r['rank']),num(r.aupr_corrected,6),int(r.geneset_n)])
table('鼠All/Pct10/Full/Up实际活动结果',['接收/靶集','配体','排名','corrected AUPR','靶数'],rows,[5.2,3.2,2.2,4,2.8],note='同一配体的Full/Drop7、BG1和Pct5结果在附录及服务器保留。排序没有推断性P。Core混合信号与ActCD4名义信号不得合并为同一机制。[A6]')
p('人Core/IL7R的IL15名义排名与鼠Core FDR上调集IL15排名3构成可讨论的跨物种候选对照；鼠Core排名1/2为Cd86/Cd274，并非IL15独占。人GZMA严格靶集优先级与上述范围不同。鼠ActCD4缺乏可计算的FDR靶集，名义结果只能作为方向探索。IL18/TNF已有旧假设和配对相关，不能因其先被研究而覆盖全候选结果。')
p('本模块把患者疾病程序、受体表达背景、所有合格巨噬发送端及完整负结果放在一起，提出可复核候选。它不把小鼠网络完成程度当作生物学优先级，也未把人鼠预测串成IKZF2→STAT1或反向通路。')

h('9  独立候选的补充敏感性与外部证据',1)
h('9.1  删除成员、原分数与背景拆分')
p('IKZF2–IL7R关联的固定高权重成员删除已完成：10%/20%删除×3seed共6规格，ρ范围−0.6216至−0.8021，P=0.001698–0.030936；相对冻结原分数的Δρ范围−0.08809至+0.01325，成对jackknife CI均包含0。关联未显示依赖这一固定少数成员；删除得分成员不是基因敲除，不支持成员非必要或实验干预效果。')
p('GZMA TReg在控制IKZF2读数后的24个既有敏感性规格及原/背景程序分解均已保留。此类调整用于检查统计关联独立性，不能把共表达调成“直接功能”。IKZF2与STAT1分别保留自己的适用亚群、RNA/AUC证据及负结果；不以能否连成串联信号决定保留哪位候选。')
h('9.2  其他发送候选的实际边界')
table('既有补充关联及拆分结果',['候选 / 接收程序','原Source关联','解释边界'],[
['MacAM HBEGF / IL7R Mitophagy','n11（5/6），ρ=+0.8477，P=0.007860，FDR=0.406863。','原目标分数Source P=0.151865；双端深度P=0.062711。不能将背景扣除后的关联直接称疾病程序。'],
['MacAM ALCAM / IL7R Mitophagy','n11（5/6），ρ=−0.7951，P=0.018329，FDR=0.597377。','原目标Source P=0.920254；双端深度P=0.093690，背景成分贡献明显。'],
['MacAM HMGB1 / GZMA TReg','n9（3/6），ρ=+0.7922，P=0.033684，FDR=0.765531。','原目标关联P=0.881532，原TReg疾病P=0.465611；不足以称疾病相关机制。'],
['MonoMac SERPINA1 / IL7R程序','混合MonoMac来源中有全族校正后的既有线索。','MacAM对应ρ=−0.0566、P=0.894；不能改称纯肺泡巨噬信号。'],
['TXLNA / 小n模型','n6近完美相关出现极小t近似P。','既有精确P=0.002778、条件置换P=0.02083；全家族FDR未通过，避免边界相关误读。']],[4.5,6,6.9],note='均为以前完成的补充结果，本轮未新增拟合或按显著性筛候选。线粒体/自噬等属于探索，不是人类主线必须通过的条件。[A7]')
h('9.3  外部人数据能补什么')
p('GSE174725的2暴露者/3矽肺患者BALF可为肺免疫状态提供独立方向对照，但不是19供者肺组织CWP队列的同质复制；小样本、采样部位和疾病标签差异均保留。GSE296682的两位患者无对照，可检查发送端定位、受体表达和细胞类型覆盖，不能检验病例–对照差异。外部数据的完整既有定位和覆盖表在服务器，不以存在可视化代替患者疾病效应。')

h('10  与相关研究比较及可支持的结论',1)
table('邻近研究与本项目的关系',['研究与真实材料','已有机制 / 发现','本项目可对应与不可外推'],[
['Cassel 2008 [1]；鼠巨噬/遗传缺陷及矽肺模型','NLRP3相关炎症小体参与矽尘反应和肺病。','为Inflam与IL18提供背景；不证明人GZMA STAT1或IL7R IKZF2由IL18驱动。'],
['Peng 2022 [2]；人BALF2+3、鼠IFNγR缺失及巨噬实验','IFNγ信号不足可加剧实验性病变。','提醒STAT1方向有细胞/阶段依赖；不能把髓系受体缺失写成CD4 STAT1缺失。'],
['You 2024 [3]；鼠驻留记忆CD4、转移与共生实验','CD69⁺CD103⁻与CD103⁺亚群功能不同；IL7参与维持。','支持细分CD4；不能直接映射到本项目GZMA/IL7R标签或替代为IL15。'],
['Zhang 2019 [4]；鼠7/14天DHI药物干预','全肺pSTAT1/pSTAT3及Th反应随多靶点药物改变。','邻近药理证据；不是细胞特异或TF特异的因果验证。'],
['Ding 2023 [5]；患者血清及巨噬外泌体/T细胞实验','外泌体miR125a-5p、TRAF6与T细胞平衡变化。','巨噬–T已有机制研究；血清外泌体不能一律归为巨噬来源。仅按已核验摘要描述，未补造样本量。'],
['Akimova 2011 [6]；人/鼠T细胞','HELIOS可与活化、增殖相关，并非唯一Treg来源标记。','IKZF2不能单独定义Treg；也不能据此否认其调控功能。'],
['Szabo 2026 [7]；人多组织，关键CRISPR为2名2–3月龄男婴脾CD3细胞','HELIOS限制部分效应程序，支持婴儿记忆T细胞状态。','与IKZF2–效应程序负关联方向相容；不是两名成人肺CWP患者的敲除验证。'],
['Afkarian 2002 [8]；鼠初始CD4、Stat1/Stat4缺失及T-bet表达','IFNγ–STAT1–T-bet与IL12R表达/分化。','支持STAT1的CD4生物学；不证明任意巨噬配体或IKZF2串联。']],[5.4,5.5,6.5])
p('目前较明确的是：存在特定来源模型下的人GZMA IKZF2 RNA信号、IL7R IKZF2调控子–程序关联，以及人GZMA STAT1调控子–Inflam关联；小鼠的RNA、AUC及部分配体预测提供相应程序的对照。它们的名义/FDR层级不同，且在人类Core疾病AUC和STAT1 RNA等位置保留真实阴性。')
p('主要缺口是疾病与来源高度耦合、亚群稀疏、多数程序和配体关联仍属探索、部分BH配对依赖低自由度设计，以及受体激活/蛋白活性/功能扰动未验证。负结果只限定所用名单、读数和模型；未计算或覆盖不足不写成“无作用”。')
p('留待讨论的目的基因阶段：分别围绕IKZF2和STAT1选择与患者证据相符的亚群、状态和验证读数，再讨论巨噬发送端是否应进入实验。可以独立检验两个TF；不预先要求IKZF2→STAT1或STAT1→IKZF2。实验体系、干预与优先级由客户决定。本轮没有启动新增调控子、Bootstrap、机制整合或大网格分析。')

h('11  可复现性、统计家族与补充结果索引',1)
h('11.1  历史模型差异与原始六候选')
p('原始六候选（Il1r2、Penk、Ttc39c、Foxp1、Tnfrsf9、Lag3）作为项目起点保留；后来关注IKZF2/STAT1没有抹去这些原发现。下表是最初Core/Full edgeR原结果的复核，不与后续保留全库offset的mechanism_step3结果拼表。原阶段在表达过滤后重新计算有效库大小（keep.lib.sizes=FALSE）；后续模型使用完整原始库offset，因此TF末位系数/P/FDR略有差别。')
q=df('original_Core_Full_focal_refit_vs_original')
table('原Core/Full候选结果（按原统计口径）',['基因','log2FC','P','原全基因FDR'],[[r.symbol_original,signed(r.logFC_original,6),num(r.PValue_original,8),num(r.FDR_original,6)] for r in q.itertuples()],[4.5,4.2,4.5,4.2],note='GSE314943，3+3；Source不存在于此鼠单队列设计。Full/NoBG是基因集合处理，All/BG1是细胞背景过滤。Penk极小P的输出精度按原存储保留，不补造有效位数。[A2]')
h('11.2  Original、DropOverlap与BG1：不选择最好种子')
q=df('old_TF_group');rows=[]
for sp,scope,model in [('Mouse','Core','Group'),('Mouse','ActCD4','Group'),('Human','Core','SourceDep')]:
 for arm in ['All','BG1']:
  for tf in ['IKZF2','STAT1']:
   v=q[(q.species==sp)&(q.scope==scope)&(q.model==model)&(q['filter']==arm)&(q.tf.str.upper()==tf)];row=[sp+'/'+scope+'/'+tf,arm]
   for variant in ['Original','DropOverlap']:
    z=v[v.variant==variant].sort_values('seed');row.append(' / '.join(f'{r.coef:+.4f} ({num(r.p,5)})' for r in z.itertuples()) or '未收录')
   rows.append(row)
table('关键疾病AUC的全部种子、评分与过滤对照',['范围/TF','过滤','Original：β (P)，17/29/43','DropOverlap：β (P)，17/29/43'],rows,[4.3,1.5,5.8,5.8],note='鼠Group、人SourceDep沿用旧SCENIC2。鼠n6；人Core All n17、BG1 n15；原未计算组间FDR，故不列空栏。三种子逐项保留，不能用其中最小P作为统一结论。[A4]')
h('11.3  小鼠NicheNet有限敏感性')
q=df('M5_focused_ligand_ranks');q=q[(q['filter']=='All')&(q.gate=='Pct10')&(q.direction=='Up')&q.test_ligand.str.upper().isin(['IL15','TNF'])];rows=[]
for (receiver,test,ligand),v in q.groupby(['receiver','test','test_ligand'],sort=False):
 row=[receiver+'/'+test,ligand]
 for omit in ['Full','Drop7']:
  z=v[v.omit==omit];row.append(f'rank {int(z.iloc[0]["rank"])}; {z.iloc[0].aupr_corrected:+.6f}; targets {int(z.iloc[0].geneset_n)}' if len(z) else '靶集不足/未计算')
 rows.append(row)
table('固定IL15/TNF：Full与Drop7对照',['接收/靶集','配体','Full：排名 / corrected AUPR / 靶数','Drop7：排名 / corrected AUPR / 靶数'],rows,[4.3,2,5.6,5.5],note='其余BG1、Pct5及Down规格完整保存。Drop7在靶和背景同时去7个既有候选，检验目标集合依赖，不是新的疾病DE或新配体筛选。[A6]')
h('11.4  家族与完整成员保存范围')
table('P/FDR及全量结果的审计规则',['分析族','家族 / 状态','本报告呈现'],[
['原RNA edgeR','每个原范围、过滤及设计中的全部受检基因；Full/NoBG和后续All/BG1各自保存。','保留原FDR；不借用新OLS FDR。'],
['M1/M2表型','比例BH在分母/模型内各亚群；程序旧修复版在scope/filter/model内12程序；扩展表另有跨scope家族。','明确使用哪一族；同P不同家族不称P变化。'],
['20261005 TF疾病OLS','readout/variant/seed/filter/model内跨scope×2TF。','保留实际FDR及未通过者；Source与SourceDep分开。'],
['旧SCENIC2疾病AUC','原规则每组≥10才计算BH；当前鼠/人组间未计算。','删除整列空FDR并表注；不写为FDR阴性。'],
['旧SCENIC2 TF–程序','species/filter/seed/variant/scope/model内所有TF×feature，可计算n≥10。','保留实际FDR与名义P。'],
['独立IKZF2 excess/删除','原固定探索未计算多重FDR；删除CI为固定模型敏感性。','不造空列；与SCENIC程序BH分开。'],
['NicheNet DE/活动/配对','DE：全受检基因BH；活动：无推断P；配对：receiver/filter/DEmodel/test/direction内跨配体×发送端。','FDR靶集、负活动分数及BH关联分开。'],
['Bootstrap','5000次/规格的有效抽样区间与成对差；没有新增P/FDR。','报告原ρ、CI、有效数及失败诊断。']],[4,7.2,6.2])
p('服务器保存：M1/M2全部504模型、前半段125组模型及30疾病比较；人类112规格状态、全配体表和不支持/不足状态；靶集、背景、ligand–target及regulon完整成员；Bootstrap 150条效应汇总、210条成对差、诊断和原5000次结果；SCENIC1–4网络/参数/覆盖、删除及方向结果。Word/PDF选取解释核心问题的聚合表图；供者逐行详情和复现归档留在服务器，不额外上传Library。')
h('11.5  来源路径与修订责任')
table('证据索引（均相对项目根目录）',['编号','路径 / 内容'],[
['A1','11_reports/human_first_20261004_v1/results；technical_delivery_20261002_ikzf2_stat1_v1/revision_v2/final_release及stepwise/front_half既有表型模型。'],
['A2','既有RNA原表、mechanism_step3及独立修订中的原counts精确复核；11_reports/independent_candidates_revision_20261005_v1/results/source_evidence。'],
['A3','11_reports/independent_candidates_revision_20261005_v1/results：human_TF疾病OLS、IL7R固定excess、成员删除、GZMA TReg独立性。'],
['A4','02_analysis/scenic1_20260929_103513_xRoEnr；scenic2_20261001_092655_ZKNWD6；scenic3_20261001_123452_Kvp9h5；scenic4_20261002_091545_iazXf1。'],
['A5','11_reports/bootstrap_TNF_IL18_20261004_v1：原脚本、抽样索引、完整效应/成对差/失败诊断与旧图。'],
['A6','11_reports/human_first_20261004_v1/results；technical_delivery最终归档及mechanism_step5_20260928_153825_Q4sB1x/tables：人/鼠NicheNet、geneset_members与prior清单。'],
['A7','missing_links、mouse_sender及independent_candidates_revision既有源证据/敏感性；本修订provenance/report_source_manifest.json逐文件SHA256。']],[1.6,15.8])
p('项目根：/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA。修订目录：11_reports/phase_delivery_IKZF2_STAT1_20261006_v2。本轮仅核查、抽取聚合证据、重写报告及重绘聚合图；科学模型、拟合、网络及Bootstrap新增运行数均为0。报告生成代码、源文件hash清单、聚合摘录与版面QA保存在修订目录。')
p('本版修复旧9页交付中的结构性遗漏：恢复背景与文献比较；完整说明SCENIC、Bootstrap、NicheNet；按分析族写明公式与深度变换；删除未计算的空FDR栏；恢复名义积极结果及对应FDR边界；独立框定IKZF2和STAT1。')
