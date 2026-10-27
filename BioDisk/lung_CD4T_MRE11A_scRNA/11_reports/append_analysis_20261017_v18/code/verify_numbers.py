from pathlib import Path
import pandas as pd,numpy as np,itertools,json,hashlib,zipfile,math
from scipy import stats
O=Path(r'C:\Users\31937\lung_append_analysis_20261017_v18');r=pd.read_csv(O/'tables/all_associations.tsv',sep='\t');s=pd.read_csv(O/'tables/sample_gene_and_program_scores.tsv',sep='\t');checks=[]
for name in ['all_associations.tsv','programme_group_effects.tsv','focal_TF_disease_effects.tsv','common_gene_sensitivity.tsv']:
 d=pd.read_csv(O/'tables'/name,sep='\t')
 for fam,z in d.groupby('family'):
  ok=z.P.notna();p=z.loc[ok,'P'].to_numpy();expected=stats.false_discovery_control(p) if len(p) else np.array([]);err=float(np.max(np.abs(expected-z.loc[ok,'FDR']))) if len(p) else 0;assert err<1e-12;checks.append(dict(file=name,check='BH_recalculation',family=fam,n=len(p),max_error=err,status='PASS'))
# Independently recompute key correlations/permutation p using explicit least-squares residuals per permutation.
for gene,sc,yy,cohort,model in [('TTC39C','ActCD4','M_ApoExec','Joint','JointStratum'),('TTC39C','ActCD4','M_ApoExec','Joint','JointStratumCondition'),('IL1R2','Treg','G_ENO1','Joint','JointStratum'),('FOXP1','Core','M_PyroExec','314','Pooled')]:
 d=s[s.variant.eq('Raw')&s.scope.eq(sc)&s.dataset.isin(['GSE314943','GSE183682'] if cohort=='Joint' else ['GSE314943'])].sort_values(['dataset','sample']).reset_index(drop=True);cx=np.ones((len(d),1));lab=d.dataset.astype(str)+'_'+d.day.astype(str) if cohort=='Joint' else pd.Series(['all']*len(d))
 if cohort=='Joint':cx=np.column_stack([cx,pd.get_dummies(lab,drop_first=True).to_numpy(float)])
 if 'Condition' in model:cx=np.column_stack([cx,d.group.eq('Silica').to_numpy(float)]);lab=lab+'_'+d.group
 xr=stats.rankdata(d['G_'+gene]);yr=stats.rankdata(d[yy]);rx=xr-cx@np.linalg.lstsq(cx,xr,rcond=None)[0];ry=yr-cx@np.linalg.lstsq(cx,yr,rcond=None)[0];rho=float(np.corrcoef(rx,ry)[0,1]);groups=[np.flatnonzero(lab.to_numpy()==b) for b in sorted(lab.unique())];ex=0;nn=0
 for parts in itertools.product(*[list(itertools.permutations(g)) for g in groups]):
  yp=yr.copy()
  for g,p in zip(groups,parts):yp[g]=yr[list(p)]
  residual=yp-cx@np.linalg.lstsq(cx,yp,rcond=None)[0];rr=np.corrcoef(rx,residual)[0,1];ex+=abs(rr)>=abs(rho)-1e-12;nn+=1
 pval=ex/nn;z=r[(r.cohort==cohort)&r.variant.eq('Raw')&r.scope.eq(sc)&r.model.eq(model)&r.x.eq('G_'+gene)&r.y.eq(yy)].iloc[0];assert abs(rho-z.rho)<1e-12 and abs(pval-z.P)<1e-12;checks.append(dict(file='all_associations.tsv',check='independent_explicit_permutation',family=f'{gene}/{sc}/{yy}/{model}',n=nn,max_error=abs(pval-z.P),status='PASS'))
# No candidate self-inclusion or shared target/program genes in the revised AUC signatures.
plan=json.loads((O/'protocol_frozen.json').read_text());cand=set(g for g,sc in plan['focus_pairs']);m=pd.read_csv(O/'tables/module_member_audit.tsv',sep='\t');reg=pd.read_csv(O/'tables/clean_regulon_members.tsv',sep='\t');assert not set(m.loc[~m.removed_candidate,'gene'])&cand;assert not set(reg.loc[reg.retained,'canonical'])&(cand|set(m.gene));assert not s.duplicated(['dataset','variant','scope','sample']).any();assert len(r)==2544 and r.P.between(0,1).all()
old=Path(r'C:\Users\31937\Downloads\线粒体焦亡_双基准复评_v17_20261017\v17_双基准二次评估报告.pdf');oldsha=hashlib.sha256(old.read_bytes()).hexdigest();j=json.loads(Path(r'C:\Users\31937\lung_dual_benchmark_20261017_v17\completion.json').read_text(encoding='utf-8'));assert j['outputs'][old.name]==oldsha
pd.DataFrame(checks).to_csv(O/'tables/numerical_QA.tsv',sep='\t',index=False);out={'status':'PASS','BH_family_checks':sum(c['check']=='BH_recalculation' for c in checks),'independent_P_checks':4,'all_candidates_removed_from_programmes':True,'all_programme_genes_and_candidates_removed_from_TF_signatures':True,'baseline_v17_pdf_sha256':oldsha,'baseline_v17_unchanged':True,'primary_associations':2544,'primary_effects':221,'TF_effects':len(pd.read_csv(O/'tables/focal_TF_disease_effects.tsv',sep='\t')),'coverage_sensitivity_associations':500};(O/'QA.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out,indent=2))