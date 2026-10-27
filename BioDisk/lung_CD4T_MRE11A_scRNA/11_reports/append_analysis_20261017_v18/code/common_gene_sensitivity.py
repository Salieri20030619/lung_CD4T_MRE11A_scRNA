from pathlib import Path
import ast,json,pandas as pd,numpy as np,itertools,math
from scipy import stats
O=Path(r'C:\Users\31937\lung_append_analysis_20261017_v18');score=pd.read_csv(O/'tables/sample_gene_and_program_scores.tsv',sep='\t');expr=pd.read_csv(O/'tables/sample_selected_gene_expression.tsv.gz',sep='\t');cov=pd.read_csv(O/'tables/module_coverage.tsv',sep='\t');r=pd.read_csv(O/'tables/all_associations.tsv',sep='\t');plan=json.loads((O/'protocol_frozen.json').read_text());affected=['MitoISR','Mitophagy','Glycolysis'];members={};audit=[]
for mod in affected:
 z=cov[cov.module==mod];common=set.intersection(*(set(s.split(';')) for s in z.mapped));members[mod]=sorted(common)
 for ds in z.dataset:audit.append(dict(module=mod,dataset=ds,common_genes=len(common),removed=';'.join(sorted(set(z.loc[z.dataset==ds,'mapped'].iloc[0].split(';'))-common)),members=';'.join(sorted(common))))
pd.DataFrame(audit).to_csv(O/'tables/common_gene_member_audit.tsv',sep='\t',index=False)
new=score.copy()
for i,row in new.iterrows():
 if row.variant not in ['Raw','BG1']:continue
 z=expr[(expr.dataset==row.dataset)&(expr.variant==row.variant)&(expr.scope==row.scope)&(expr['sample']==row['sample'])].set_index('gene')
 for mod,genes in members.items():
  vals=z.reindex(genes).CP10k.to_numpy();assert np.isfinite(vals).all();new.at[i,'M_'+mod]=float(np.log1p(vals).mean())
for i,row in new[new.variant.isin(['FixedComp','BG1FixedComp'])].iterrows():
 base='Raw' if row.variant=='FixedComp' else 'BG1';z=new[(new.dataset==row.dataset)&(new.variant==base)&new.scope.isin(['Naive','ActCD4','Treg'])&(new['sample']==row['sample'])];assert len(z)==3
 for mod in affected:new.at[i,'M_'+mod]=float(z['M_'+mod].mean())
new.to_csv(O/'tables/common_gene_sample_scores.tsv',sep='\t',index=False)
# Reuse the same frozen statistic functions; only these three panels change.
module=ast.parse((O/'analyze.py').read_text(encoding='utf-8'));funcs=[n for n in module.body if isinstance(n,ast.FunctionDef) and n.name in ['design','blocks','permutations','coefficient','assoc','finalize']];exec(compile(ast.Module(body=funcs,type_ignores=[]),'frozen_statistics','exec'))
pcache={};results=[];testids=set();aa=pd.read_csv(O/'tables/clean_TF_sample_AUC.tsv',sep='\t');d0=r[r.y.isin(['M_'+v for v in affected])]
for row in d0.itertuples():
 ds=['GSE314943'] if row.cohort=='314' else ['GSE314943','GSE183682'] if row.cohort=='Joint' else ['GSE274372'];d=new[new.dataset.isin(ds)&new.variant.eq(row.variant)&new.scope.eq(row.scope)].sort_values(['dataset','sample']).reset_index(drop=True)
 if row.cohort=='274_min10':d=d[d.cells>=10].reset_index(drop=True)
 if row.x.startswith('TF_'):
  ap=aa[(aa.variant==row.variant)&(aa.scope==row.scope)].pivot(index='sample',columns='tf',values='AUC').add_prefix('TF_');d=d.merge(ap,on='sample',how='left')
 assoc(d,row.x,row.y,row.cohort,'CommonGene_'+row.variant,row.kind,row.scope,row.model,row.sender,{'gene':row.gene,'program':row.program,'original_variant':row.variant,'sensitivity_reason':'mapping_coverage_discrepancy'})
x=finalize(results,'common_gene_sensitivity.tsv')
record={'reason':'Three panels have unequal mapped membership across cohorts; all affected comparisons rescored on common genes to evaluate comparability. Primary results unchanged.','affected_panels':affected,'planned_comparisons':len(x),'status':'post-preflight mapping sensitivity, not replacement of primary analysis','family_rule':'Separate same-kind/cohort/variant/model families for this sensitivity; not used to relabel original FDR.'};(O/'coverage_sensitivity_record.json').write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8')
print('COMMON_GENES',[(k,len(v)) for k,v in members.items()],'TESTS',len(x));print(x[(x.cohort=='274')&(x.model=='TimeDayAdjusted')&x.gene.eq('TTC39C')][['program','rho','P','FDR']].to_string(index=False))