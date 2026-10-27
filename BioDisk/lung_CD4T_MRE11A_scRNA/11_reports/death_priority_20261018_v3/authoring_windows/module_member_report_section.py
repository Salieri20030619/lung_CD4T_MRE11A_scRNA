new('程序成员移除区分了稳定代谢关联与成员敏感的损伤线索')
mm=rd('module_member_sensitivity_summary.tsv')
zz=mm[(mm.variant=='Raw')&(mm.model=='JointStratum')]
table(['候选/程序','ρ范围','仍为P＜0.05','最减弱成员'],
      [[r.gene+' / '+r.program,fmt(r.rho_min)+' - '+fmt(r.rho_max),str(int(r.nominal_omissions))+'/'+str(int(r.members)),r.most_attenuating_member] for r in zz.itertuples()],
      [165,125,105,104])
para('五项关系在逐一移除评分成员后均保留正向相关，但名义显著性稳定程度不同。IL2RA-糖酵解为117/117；HIF1A-PyroExec为2/5，ARL6IP1-MitoROS为7/12。后两者保留探索性定位，并明确成员敏感性。')
para('IL2RA-ApoExec主分析本身未达P＜0.05；不能因某些成员移除后P变小而改写主结论。OXPHOS为19/24名义支持，仍需结合疾病调整后减弱的结果。')
para('本项为基于前述发现追加的成员敏感性，不是重新筛选。Raw/BG1和两个模型共660项检验，20个完整家族分别校正；原主分析P/FDR不变。移除次数不构成新的生物学重复。','small')
