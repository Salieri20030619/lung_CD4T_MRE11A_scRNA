from pathlib import Path
import pandas as pd,numpy as np,itertools,json,re,hashlib
from scipy.stats import rankdata,false_discovery_control
from scipy.io import mmread
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/four_axis_screen_20261017_v15';M=P/'02_analysis/mechanism_step1_20260927_152218_s4V5eF';S=['Saline_1','Saline_2','Saline_3','Silica_1','Silica_2','Silica_3'];rd=lambda p:pd.read_csv(p,sep='\t')
assert hashlib.sha256((O/'PLAN_FROZEN.json').read_bytes()).hexdigest()=='1852715bf9c6a5a38aae4e0dae903833a40b7562df26902cef06989894073ac5'
checks=[]
for f in [O/'03_damage/all_gene_program_tests.tsv.gz']+list((O/'04_macrophage').glob('*_all_links.tsv.gz')):
 if f.name=='consolidated_all_links.tsv.gz':continue
 a=rd(f)
 for fam,b in a.groupby('family'):
  ok=b.P.notna();q=false_discovery_control(b.loc[ok,'P'].to_numpy()) if ok.any() else np.array([]);err=float(np.max(abs(q-b.loc[ok,'BH'].to_numpy()))) if len(q) else 0;checks.append(dict(file=f.name,family=fam,total=len(b),valid=int(ok.sum()),max_BH_error=err,pass_BH=err<1e-10))
assert all(c['pass_BH'] for c in checks)
pd.DataFrame(checks).to_csv(O/'provenance/BH_independent_audit.tsv',sep='\t',index=False)
# Independent selected exact tests from original raw count matrices.
md=rd(M/'comm_input/CD4_cells.tsv');gn=rd(M/'comm_input/CD4_genes.tsv');g=pd.Index(gn.canonical);x=mmread(M/'comm_input/CD4_counts.mtx.gz').T.tocsr();cp=x.multiply((1e4/md.library.to_numpy())[:,None]).tocsr();lg=cp.copy();lg.data=np.log1p(lg.data);mod=rd(P/'11_reports/methods_complete_20261016_v13/supplement/S20_all_fixed_module_genes.tsv');ctrl=rd(P/'11_reports/methods_complete_20261016_v13/supplement/S20_all_fixed_module_controls.tsv.gz');ctrl=ctrl[ctrl.dataset=='GSE314943'];tab=rd(O/'03_damage/all_gene_program_tests.tsv.gz');results=[]
def exact(xx,yy,adj):
 rx=rankdata(xx);ry=rankdata(yy)
 for inds in ([range(3),range(3,6)] if adj else [range(6)]):
  inds=list(inds);rx[inds]-=rx[inds].mean();ry[inds]-=ry[inds].mean()
 den=np.linalg.norm(rx)*np.linalg.norm(ry)
 if den<1e-12:return np.nan,np.nan
 r=float(np.dot(rx,ry)/den);perm=(list(a)+list(b) for a in itertools.permutations(range(3)) for b in itertools.permutations(range(3,6))) if adj else itertools.permutations(range(6));values=[np.dot(rx[list(v)],ry)/den for v in perm];p=sum(abs(v)>=abs(r)-1e-12 for v in values)/len(values);return r,p
for sc in ['Naive','ActCD4','Treg']:
 for gene,module in [('GSDMD','PyroExec'),('ATF4','MitoISR'),('GPX4','MitoROS'),('STAT1','MitoROS')]:
  rows=tab[(tab.scope==sc)&(tab.gene==gene)&(tab.module==module)&(tab.branch=='All')]
  if rows.empty:continue
  members=set(mod.loc[mod.module==module,'gene'])-{gene};controls=set(ctrl.loc[ctrl.module==module,'gene'])-{gene};mi=g.get_indexer(sorted(members));mi=mi[mi>=0];ci=g.get_indexer(sorted(controls));ci=ci[ci>=0];xx=[];yy=[]
  for sample in S:
   m=md['sample'].eq(sample).to_numpy()&md.subtype.eq(sc).to_numpy();xx.append(float(cp[m,g.get_loc(gene)].mean()));yy.append(float(lg[m][:,mi].mean()-lg[m][:,ci].mean()))
  for adj in [False,True]:
   r,p=exact(xx,yy,adj);ref=rows[rows.model==('Adjusted' if adj else 'Pooled')].iloc[0];ok=(np.isnan(p) and pd.isna(ref.P)) or (abs(p-ref.P)<1e-10 and abs(r-ref.rho)<1e-10);assert ok,(sc,gene,module,r,p,ref.P);results.append(dict(scope=sc,gene=gene,module=module,model=ref.model,rho=r,P=p,status='PASS'))
pd.DataFrame(results).to_csv(O/'provenance/raw_count_exact_test_audit.tsv',sep='\t',index=False)
# Count rosters against frozen historical BG1, before any new quality filters.
old=rd(M/'tables/sample_scores.tsv');old=old[(old.dataset=='GSE314943')&old['filter'].isin(['All','BG1'])];old['scope']=old.scope.replace({'All':'Core'});new=rd(O/'01_quality/sample_scope_roster.tsv');j=old.merge(new,left_on=['scope','filter','sample'],right_on=['scope','branch','sample'],suffixes=('_old','_new'));assert len(j)==48 and (j.cells_old==j.cells_new).all()
status={'BH_families':len(checks),'BH_all_pass':True,'independent_raw_count_exact_tests':len(results),'original_BG1_All_sample_counts_reproduced':len(j),'frozen_plan_hash_unchanged':True,'status':'PASS'};(O/'provenance/numeric_validation.json').write_text(json.dumps(status,indent=2),encoding='utf-8');print(json.dumps(status))