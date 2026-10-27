# Executed inside report.py after the old-network comparison section.
jqa=json.loads((A/'joint_network/QA_complete.json').read_text(encoding='utf-8'))
jcov=rd('joint_network_TF_coverage.tsv')
jassoc=rd('joint_network_TF_associations.tsv')
jsum=rd('joint_network_candidate_summary.tsv')
fig('JN01','完整5+5重新建网完成两个范围、三个种子',
    f'All与BG1分别完成种子17/29/43的GRNBoost2及motif剪枝。至少2/3种子重复、去重叠后可评分regulon分别为{jqa["All_scored"]}和{jqa["BG1_scored"]}个。',
    '每次训练保留全部10个独立样本，每样本最多250个细胞；样本较少时不重复抽样。全体保留细胞用于AUCell评分。2/3稳定性反映细胞下采样/随机种子复现，不是新增独立样本验证。')
fig('JN02','当前联合网络显示TF活性变化',
    '新网络使用共同测得基因和统一motif资源重新推断，再在独立样本层面检验疾病效应。',
    'All/BG1各自重新建网；稳定边至少2/3种子，候选/程序成员及TF自身从评分签名移除。*表示原始P＜0.05；完整regulon家族FDR保留。')
fig('JN03','稳定调控先验与候选样本关联分别呈现',
    '图中先按种子支持和平均权重为每个主候选展示至多5条稳定motif边，再叠加同亚群RNA-活性关联结果。',
    '连线粗细区分2/3与3/3种子支持；红线为同亚群、队列/时间调整后原始P＜0.05，灰线为未支持或不可评估。边的存在不等于直接结合或因果调控。')
new('联合新网络的主候选TF证据与敏感性')
rows=[]
for r in jsum[jsum.branch=='primary'].itertuples():
 rows.append([r.gene+' / '+r.scope,r.variant,int(r.tested),int(r.nominal),int(r.FDR)])
table(['候选/亚群','范围','可评估连接','P＜0.05','FDR＜0.05'],rows,[135,70,100,95,99])
para(f'联合新建网络独立复核通过{jqa["checks"]}项；其中{jqa["independent_primary_TF_tests"]}项主候选RNA/TF效应与P值完成独立复算。原52项核查继续有效。')
z=jassoc[(jassoc.variant=='Raw')&(jassoc.model=='JointStratum')&jassoc.kind.eq('jointnet_TF_candidate_primary')&(jassoc.P<.05)].sort_values('P')
if len(z):
 para('主分析名义连接如下；全部结果及条件调整、BG1敏感性保留在统计表。','small')
 table(['候选','TF','ρ','P / FDR'],[[r.gene,r.tf,fmt(r.rho),fmt(r.P)+' / '+fmt(r.FDR)] for r in z.head(6).itertuples()],[95,100,85,219])
else:
 para('当前联合新网络的主分析仍没有达到P＜0.05的主候选RNA-TF连接；这一结果加强了“调控连接尚未闭合”的定位，不能用结构先验取代样本统计。')

new('IL2RA-NFKB2名义关联在疾病调整后未保留')
rows=[]
for var in ['Raw','BG1']:
 for model in ['JointStratum','JointStratumCondition']:
  r=jassoc[(jassoc.variant==var)&(jassoc.model==model)&(jassoc.gene=='IL2RA')&(jassoc.scope=='Treg')&(jassoc.tf=='NFKB2')].iloc[0]
  rows.append([var,'队列/时间＋疾病' if 'Condition' in model else '队列/时间',fmt(r.rho),fmt(r.P),fmt(r.FDR)])
table(['范围','调整变量','ρ','P','FDR'],rows,[65,140,95,95,104])
para('主分析IL2RA-NFKB2为名义关联，未过FDR；BG1范围通过其完整家族FDR。进一步调整疾病组别后，两个范围均未保留P＜0.05。它支持伴随疾病状态变化的调控线索，尚不支持独立于疾病组别的稳定连接。')
para('NFKB1和NFKB2的Treg疾病活性变化均为正向、原始P=0.025，主分析FDR=0.350。程序成员与候选已从regulon评分移除，关联仍属同批样本中的探索性评价。')
new('NFKB2更偏向Treg活化与应激关联，焦亡程序未获支持')
pp=jassoc[(jassoc.scope=='Treg')&(jassoc.variant=='Raw')&(jassoc.model=='JointStratum')&(jassoc.tf=='NFKB2')&jassoc.kind.eq('jointnet_focalTF_program')]
table(['程序','ρ','P','FDR'],[[mod,fmt(pp[pp.program==mod].iloc[0].rho),fmt(pp[pp.program==mod].iloc[0].P),fmt(pp[pp.program==mod].iloc[0].FDR)] for mod in ['TAct','MitoISR','ApoExec','PyroExec']],[150,110,115,124])
para('同来源上游比较：AllMac发送端、Treg接收端、同一批10样本，均为队列/时间调整主模型。','small')
bridge=rd('IL2RA_NFKB2_matched_source_bridge.tsv')
table(['配体','配体-IL2RA ρ','P / FDR','配体-NFKB2 ρ','P / FDR'],
 [[r.ligand,fmt(r.ligand_IL2RA_rho),fmt(r.ligand_IL2RA_P)+' / '+fmt(r.ligand_IL2RA_FDR),fmt(r.ligand_NFKB2_rho),fmt(r.ligand_NFKB2_P)+' / '+fmt(r.ligand_NFKB2_FDR)] for r in bridge.itertuples()],
 [65,100,115,100,119])
para('IL18的两个负向关联与SPP1/IL6/IGF1的部分正向关联分别保留；不能跨配体拼接显著段，也不能把同时变化解释为配体刺激导致候选或死亡改变。','small')
