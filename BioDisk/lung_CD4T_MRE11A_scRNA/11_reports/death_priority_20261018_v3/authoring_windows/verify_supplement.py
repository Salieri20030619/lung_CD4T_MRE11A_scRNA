from pathlib import Path
import json,hashlib,itertools,math
import pandas as pd,numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
P=Path(__file__).parent;T=P/'tables';checks=[];nums=[]
def ok(name,test,detail=''):
 checks.append(dict(check=name,passed=bool(test),detail=detail))
 if not test:raise AssertionError(name+': '+detail)
def bh(p):return multipletests(np.asarray(p,float),method='fdr_bh')[1]
for f in sorted(T.glob('*_DE.tsv')):
 if f.name=='human_G227_candidate_DE.tsv':continue
 x=pd.read_csv(f,sep='\t')
 if 'PValue' not in x or 'FDR' not in x:continue
 ok(f.name+' unique_gene',x.gene.is_unique)
 good=x.PValue.notna();ok(f.name+' complete_family_BH',np.allclose(bh(x.loc[good,'PValue']),x.loc[good,'FDR'],atol=1e-10,rtol=1e-7))
ef=pd.read_csv(T/'external_module_effects.tsv',sep='\t');s=pd.read_csv(T/'external_module_scores.tsv',sep='\t')
for r in ef.itertuples():
 z=s[(s.contrast==r.contrast)&(s.program==r.program)]
 group=z.group.eq('KO').to_numpy(float)
 X=np.column_stack([np.ones(len(z)),pd.get_dummies(z.donor,drop_first=True,dtype=float).to_numpy(),group]) if r.contrast.startswith('GSE324038') else np.column_stack([np.ones(len(z)),group])
 y=z.score.to_numpy();coef=np.linalg.lstsq(X,y,rcond=None)[0];res=y-X@coef;df=len(z)-np.linalg.matrix_rank(X);se=np.sqrt((res@res)/df*np.linalg.inv(X.T@X)[-1,-1]);p=2*stats.t.sf(abs(coef[-1]/se),df)
 nums.append(dict(contrast=r.contrast,program=r.program,effect_independent=coef[-1],P_independent=p,effect_error=coef[-1]-r.effect,P_error=p-r.P))
ok('external 42 program effects independent linear refit',len(nums)==42 and max(abs(x['P_error']) for x in nums)<1e-9 and max(abs(x['effect_error']) for x in nums)<1e-9)
for family,z in ef.groupby('family'):ok('external module BH '+family,np.allclose(bh(z.P),z.FDR,atol=1e-10))
pd.DataFrame(nums).to_csv(T/'independent_external_module_QA.tsv',sep='\t',index=False)
for name in ['external_module_GSEA.tsv','external_joint_regulon_GSEA.tsv','mouse_HIF1A_perturbation_GSEA.tsv','GSE121883_thymus_KEGG_GSEA.tsv','GSE121883_spleen_KEGG_GSEA.tsv']:
 d=pd.read_csv(T/name,sep='\t')
 for fam,z in d.groupby('family'):
  ok('GSEA '+fam+' finite P',z.pval.notna().all() and z.pval.between(0,1).all())
  ok('GSEA '+fam+' BH',np.allclose(bh(z.pval),z.padj,atol=1e-10,rtol=1e-7))
for name in ['human_candidate_associations.tsv','human_targeted_member_associations.tsv']:
 d=pd.read_csv(T/name,sep='\t')
 for family,z in d.groupby('family'):
  z=z[z.P.notna()];ok(name+' BH '+family,np.allclose(bh(z.P),z.FDR,atol=1e-10))
h=pd.read_csv(T/'human_sample_scores.tsv',sep='\t');a=pd.read_csv(T/'human_candidate_associations.tsv',sep='\t')
for r in a[(a.cohort=='G227')&(a.min_cells==10)&a.model.eq('Source')].itertuples():
 z=h[(h.series=='G227')&(h.subtype==r.scope)&(h.cells>=10)]
 rho=stats.spearmanr(z['G_'+r.gene],z[r.program]).statistic
 ok('human independent rho '+r.gene+' '+r.program,abs(rho-r.rho)<1e-10)
bm=pd.read_csv(T/'BAL174725_descriptive_effects.tsv',sep='\t');ok('BAL no invented P',bm.P.isna().all())
m=pd.read_csv(P/'analysis/tables/module_members.tsv',sep='\t');c=pd.read_csv(P/'analysis/tables/frozen_candidate_pairs.tsv',sep='\t')
ok('all 7 candidates excluded from modules',not (m.retained.astype(str).str.lower().eq('true')&m.gene.isin(c.gene)).any())
q=pd.read_csv(T/'external_fit_QC.tsv',sep='\t').set_index('contrast')
ok('human four independent donors not eight',q.loc['GSE324038_HumanPanT_Hypoxia','n_biological']==4)
ok('IL2RA independent 8 and 10',q.loc['GSE121883_spleen','n_biological']==8 and q.loc['GSE121883_thymus','n_biological']==10)
# Original report preserved byte-for-byte.
v2=Path(r'C:\Users\31937\Downloads\矽肺_凋亡焦亡优先_5加5证据补全_v2_20261018\矽肺_CD4凋亡焦亡优先_5加5生信证据补全_v2.pdf')
ok('v2 report preserved',hashlib.sha256(v2.read_bytes()).hexdigest()=='a9d113c9bc7d8c163d970528d43fb3425ffa13954a9f30fcd4d749c97e99c17a')
network=(T/'network_LOO_candidate_associations.tsv').exists()
if network:
 d=pd.read_csv(T/'network_LOO_candidate_associations.tsv',sep='\t')
 for fam,z in d.groupby('family'):
  z=z[z.P.notna()]
  if len(z):ok('LOO BH '+fam,np.allclose(bh(z.P),z.FDR,atol=1e-10))
 audit=pd.read_csv(T/'network_LOO_fit_audit.tsv',sep='\t');ok('twenty completed network folds',len(audit)==20 and audit.groupby('arm').size().eq(10).all())
 ok('LOO nine independent mice',d.loc[d.status=='OK','n'].eq(9).all())
if network:
 exec((P/'verify_network_independent.py').read_text(encoding='utf-8'))
 exec((P/'verify_reference_independent.py').read_text(encoding='utf-8'))
pd.DataFrame(checks).to_csv(T/'supplement_QA_checks.tsv',sep='\t',index=False)
result={'status':'success' if network else 'external_pass_network_pending','checks':len(checks),'independent_module_refits':len(nums),'network_checked':network,'failed':sum(not x['passed'] for x in checks)}
(P/'supplement_QA.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result))
