from pathlib import Path
import json,pandas as pd,numpy as np
D=Path('manuscript_delivery');E=json.load(open('revision_delivery/report_evidence_aggregate.json',encoding='utf8'));F=json.load(open(D/'added_evidence_aggregate.json',encoding='utf8'))
checks=[]
def check(name,actual,expected,tol=1e-6):
 ok=bool(np.isclose(actual,expected,atol=tol,rtol=tol));checks.append({'name':name,'actual':float(actual),'expected':float(expected),'pass':ok})
d=pd.DataFrame(F['all_TF_program_associations'])
for ar,feature,expected in [('All','PyroExec',[[-.8299822646,.0004460527,.0256710373],[-.8378323756,.0003494645,.0288308242],[-.6637767864,.0133649723,.1107740823]]),('All','ApoExec',[[-.6866216386,.0095375802,.0825990437],[-.7266373105,.0049019064,.0833631335],[-.6075152115,.0276431687,.1546143334]])]:
 for seed,vals in zip([17,29,43],expected):
  r=d[(d.species=='Human')&(d.scope=='Core')&(d.model=='SourceDep')&(d.tf=='IKZF2')&(d['filter']==ar)&(d.variant=='DropOverlap')&(d.seed==seed)&(d.feature=='score_'+feature)].iloc[0]
  for key,ex in zip(['rho','p','FDR'],vals):check(f'IKZF2 {ar} {feature} {seed} {key}',r[key],ex)
  check(f'{feature} n seed{seed}',r['n'],17)
u=pd.DataFrame(F['human_UCell_module_group_tests']);check('human UCell 36 comparisons',len(u),36);check('human UCell nominal count',(u.p_value<.05).sum(),0)
s=pd.DataFrame(F['human_GSVA_summary']);check('human GSVA BH discoveries',s.BH05.sum(),0)
g=pd.DataFrame(F['human_mitochondrial_pathway_stability']);check('stable mito comparisons',g.all_fdr_lt_005.sum(),14)
m=pd.DataFrame(F['mouse_UCell_score_group_tests']);check('mouse UCell 32 comparisons',len(m),32);check('mouse UCell nominal count',(m.PValue<.05).sum(),0)
nn=pd.DataFrame(E['human_NN_main']);r=nn[(nn.receiver=='Core')&(nn.model=='Source')&(nn.test=='Nominal05')&(nn.direction=='Up')&(nn.test_ligand=='IL15')].iloc[0];check('human core IL15 AUPR',r.aupr_corrected,.012640,1e-6)
check('lineage total',sum([85473,40770,29503,14853,8626,7954,4854,542]),192575)
r=pd.DataFrame(E['human_TF_disease']).query("readout=='RNA' and scope=='CD4-GZMA' and tf=='IKZF2' and model=='Source'").iloc[0];check('human GZMA RNA P',r.P,.002109,1e-6);check('human GZMA RNA BH',r.FDR,.025304,1e-6)
(D/'QA/numeric_audit.json').write_text(json.dumps({'passed':sum(x['pass'] for x in checks),'total':len(checks),'checks':checks},indent=2),encoding='utf8');print(json.dumps({'passed':sum(x['pass'] for x in checks),'total':len(checks),'failed':[x for x in checks if not x['pass']]}))
