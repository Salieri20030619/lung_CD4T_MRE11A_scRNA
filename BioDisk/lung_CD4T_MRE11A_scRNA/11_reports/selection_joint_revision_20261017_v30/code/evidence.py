from pathlib import Path
import pandas as pd,numpy as np,json,itertools,zipfile,io,hashlib
from scipy.stats import false_discovery_control
H=Path(r'C:\Users\31937');O=H/'lung_selection_joint_revision_20261017_v30';O.mkdir(exist_ok=True)
for n in ['evidence','tables','plot_data','figures']:(O/n).mkdir(exist_ok=True)
B=H/'lung_three_cohort_20261017_v16';z=zipfile.ZipFile(B/'data/resources/legacy_v13_evidence.zip');aa=pd.read_csv(io.BytesIO(z.read('S17_table_02.tsv')),sep='\t');bb=pd.read_csv(io.BytesIO(z.read('S17_table_06.tsv')),sep='\t');xx=set(aa['候选基因']);yy=set(bb['基因']);extra={'DUSP5','GSDMD','CXCR6','MRE11'}
assert len(xx)==11 and len(yy)==27 and not xx&yy and len(extra-(xx|yy))==4
sup=pd.read_csv(B/'final_tables/历史候选逐项回溯.tsv',sep='\t');assert xx|yy|extra==set(sup.gene)
rows=[{'gene':g,'origin_group':'候选比较对象11' if g in xx else '高效应核查信号27' if g in yy else '定向补核对象4','source_table':'S17_table_02' if g in xx else 'S17_table_06' if g in yy else 'finalize_tables.py explicit set','source_role':'candidate comparison' if g in xx else 'review high effect and cell origin' if g in yy else 'targeted follow-up'} for g in sorted(set(sup.gene))]
pd.DataFrame(rows).merge(sup,on='gene',validate='one_to_one').to_csv(O/'tables/supplement42_provenance.tsv',sep='\t',index=False)
for fn in ['S17_table_02.tsv','S17_table_05.tsv','S17_table_06.tsv']:(O/'evidence'/fn).write_bytes(z.read(fn))
# Reproduce existing group scores with explicitly stratified label enumeration.
def test(z,y,joint):
 z=z.reset_index(drop=True);v=z[y].to_numpy(float);cond=z.group.eq('Silica').to_numpy(float)
 labs=(z.dataset.astype(str)+'_'+z.day.astype(str)) if joint else pd.Series(['all']*len(z))
 C=np.column_stack([np.ones(len(z)),pd.get_dummies(labs,drop_first=True,dtype=float).to_numpy()]);M=np.eye(len(z))-C@np.linalg.pinv(C);xc=M@cond;rv=M@v;effect=float(xc@rv/(xc@xc))
 groups=[np.flatnonzero(labs.eq(g)) for g in sorted(labs.unique())];choices=[list(itertools.combinations(g,int(cond[g].sum()))) for g in groups];null=[]
 for parts in itertools.product(*choices):
  q=np.zeros(len(z))
  for ids in parts:q[list(ids)]=1
  q=M@q;null.append(float(q@rv/(q@q)))
 return effect,float(np.mean(np.abs(null)>=abs(effect)-1e-12)),len(null)
a=pd.read_csv(H/'lung_append_analysis_20261017_v18/tables/sample_gene_and_program_scores.tsv',sep='\t');e=pd.read_csv(H/'lung_append_analysis_20261017_v18/tables/programme_group_effects.tsv',sep='\t');checks=[]
for co in ['314','Joint']:
 z0=a[a.variant.eq('Raw')&a.dataset.isin(['GSE314943'] if co=='314' else ['GSE314943','GSE183682'])]
 fam=e[e.variant.eq('Raw')&e.cohort.eq(co)]
 for r in fam.itertuples():
  z=z0[z0.scope.eq(r.scope)].sort_values(['dataset','sample']);ef,pv,n=test(z,'M_'+r.program,co=='Joint');assert abs(ef-r.effect)<1e-10 and abs(pv-r.P)<1e-12
  checks.append({'cohort':co,'scope':r.scope,'program':r.program,'effect':ef,'P':pv,'original_FDR':r.FDR,'allocations':n,'family':r.family})
c=pd.DataFrame(checks)
for co,z in c.groupby('cohort'):assert np.allclose(false_discovery_control(z.P),z.original_FDR)
c.to_csv(O/'evidence/group_score_recalculation.tsv',sep='\t',index=False)
e[e.variant.eq('Raw')&e.cohort.isin(['314','Joint'])].to_csv(O/'tables/group_scores_both_models.tsv',sep='\t',index=False)
# New supplementary proportion comparison, same sample-level fraction, 3-test family per model.
i=pd.read_csv(B/'data/tables/sample_inventory.tsv',sep='\t');i=i[i.dataset.isin(['GSE314943','GSE183682'])&i.scope.isin(['Naive','ActCD4','Treg','Core'])]
w=i.pivot(index=['dataset','sample','group','day'],columns='scope',values='cells').reset_index();assert len(w)==10;assert (w[['Naive','ActCD4','Treg']].sum(axis=1)==w.Core).all()
for sc in ['Naive','ActCD4','Treg']:w[sc+'_fraction']=w[sc]/w.Core
w.to_csv(O/'plot_data/F04.tsv',sep='\t',index=False);pr=[]
for co in ['314','Joint']:
 z=w if co=='Joint' else w[w.dataset.eq('GSE314943')]
 for sc in ['Naive','ActCD4','Treg']:
  ef,pv,n=test(z,sc+'_fraction',co=='Joint');pr.append({'cohort':co,'scope':sc,'effect':ef,'P':pv,'allocations':n,'n':len(z),'family':co+'|CD4_proportion_3states'})
pr=pd.DataFrame(pr)
for co,z in pr.groupby('cohort'):pr.loc[z.index,'FDR']=false_discovery_control(z.P)
pr.to_csv(O/'tables/proportion_joint_supplement.tsv',sep='\t',index=False);print('PROPORTIONS',pr.to_string(index=False));print('SOURCE42',11,27,4,'GROUP_SCORE_QA',len(c))
# Named nine-gene panel and relevant complete-model statistics, not new filtering.
assoc=pd.read_csv(H/'lung_append_analysis_20261017_v18/tables/all_associations.tsv',sep='\t');spec=[('IL1R2','Treg','Inhib'),('DDIT4','ActCD4','PyroExec'),('SKIL','Treg','TAct'),('FOXP1','Core','PyroExec'),('GSN','ActCD4','ApoExec'),('AHNAK','Naive','OXPHOS'),('TTC39C','ActCD4','ApoExec'),('MDFIC','ActCD4','ApoExec'),('ST6GALNAC3','ActCD4','MitoROS')];parts=[]
for g,sc,mod in spec:
 q=assoc[assoc.variant.eq('Raw')&assoc.cohort.eq('Joint')&assoc.gene.eq(g)&assoc.scope.eq(sc)&assoc.program.eq(mod)&assoc.model.isin(['JointStratum','JointStratumCondition'])];assert len(q)==2;parts.append(q)
pd.concat(parts).to_csv(O/'tables/nine_gene_comparison_statistics.tsv',sep='\t',index=False)
meta={'source42_counts':[11,27,4],'source42_exact_match':True,'group_score_checks':len(c),'group_score_values_unchanged':True,'main_3plus3_assignments':20,'joint_5plus5_stratified_assignments':80,'two_sided_min_P':[.1,.025],'new_analysis':'3 joint sample-level CD4 proportion tests; 3-state BH family; original n6 reproduced','candidate_panel':'5 mechanism leads plus4 exploratory nominees fixed before additional association calculations; not automatic top9 among90'};(O/'evidence/analysis_scope.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding='utf-8')