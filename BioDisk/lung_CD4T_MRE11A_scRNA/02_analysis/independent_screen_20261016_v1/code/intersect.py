from pathlib import Path
import pandas as pd,numpy as np,json
O=Path('C:/codex_test/output/silicosis_independent_v1/independent_screen_20261016_v1')
assert (O/'COMPLETED').exists() and (O/'exit_status').read_text().strip()=='0'
t=pd.read_csv(O/'trajectory/all_gene_and_program_tests.tsv.gz',sep='\t')
r=pd.read_csv(O/'time/Omnibus.tsv',sep='\t')
alias={'MRE11A':'MRE11','ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','PARK2':'PRKN'}
r['feature']=r.symbol.str.upper().replace(alias)
assert r.feature.notna().all() and r.feature.is_unique
r=r.rename(columns={'PValue':'time_P','FDR':'time_FDR'})
for c in ['D7_D0','D28_D0','D56_D0','D28_D7','D56_D28']:
 z=pd.read_csv(O/f'time/{c}.tsv',sep='\t')[['gene','logFC','PValue','FDR']].rename(columns={k:c+'_'+k for k in ['logFC','PValue','FDR']});r=r.merge(z,on='gene',validate='one_to_one')
g=t[t.type=='gene'];p=g[g.run=='K15_R1'].copy();stable=g.assign(passP=g.P.lt(.05)).groupby(['scope','feature']).agg(configurations=('run','nunique'),nominal_configurations=('passP','sum')).reset_index();p=p.merge(stable,on=['scope','feature'],validate='one_to_one');allr=p.merge(r,on='feature',validate='many_to_one');hit=allr[(allr.P<.05)&(allr.time_P<.05)].copy()
hit['D28_D7_endpoint_agreement']=np.sign(hit.endpoint_difference)==np.sign(hit.D28_D7_logFC)
hit['time_reversal']=np.sign(hit.D28_D7_logFC)!=np.sign(hit.D56_D28_logFC)
hit=hit.sort_values(['scope','nominal_configurations','curve_range'],ascending=[True,False,False])
(O/'intersection').mkdir(exist_ok=True)
allr.to_csv(O/'intersection/all_jointly_tested_genes.tsv',sep='\t',index=False,encoding='utf-8');hit.to_csv(O/'intersection/nominal_intersection.tsv',sep='\t',index=False,encoding='utf-8')
rows=[]
for scope,b in p.groupby('scope'):
 h=hit[hit.scope==scope];j=allr[allr.scope==scope]
 rows.append(dict(scope=scope,trajectory_tested=len(b),trajectory_nominal=int(b.P.lt(.05).sum()),jointly_tested=len(j),intersection=len(h),stable4=int(h.nominal_configurations.eq(4).sum()),both_FDR05=int(((h.BH<.05)&(h.time_FDR<.05)).sum())))
summary=pd.DataFrame(rows);summary.to_csv(O/'intersection/screen_counts.tsv',sep='\t',index=False)
print(summary.to_string(index=False));print('TIME',len(r),int(r.time_P.lt(.05).sum()));print('TOP\n',hit.groupby('scope').head(12)[['scope','feature','curve_range','P','BH','time_P','time_FDR','nominal_configurations']].to_string(index=False))
print('SCORES\n',pd.read_csv(O/'scores/subtype_rankings.tsv',sep='\t').query('case_score_rank==1')[['module','scope','case_mean','difference']].to_string(index=False));print('PROGRAMS\n',t[(t.run=='K15_R1')&(t.type=='program')][['scope','feature','P','BH','endpoint_difference']].to_string(index=False))
