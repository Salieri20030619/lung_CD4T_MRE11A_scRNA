from pathlib import Path
import pandas as pd,numpy as np,itertools,json
from scipy import stats
ROOT=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/stepwise_diagnostic_20261002_v1')
d=pd.read_csv(ROOT/'results/sender_all_paired_samples.tsv.gz',sep='\t');d=d[d.variant=='Matched']
rows=[]
for (arm,scope,seed,tf),z in d.groupby(['filter','scope','seed','tf']):
 for before,after in [(1,5),(5,10),(10,20)]:
  b=z[(z.receiver_cells>=before)&(z.sender_cells>=before)]
  end=b[(b.receiver_cells>=after)&(b.sender_cells>=after)]
  removed=sorted(set(b['sample'])-set(end['sample']))
  for k in range(len(removed)+1):
   for subset in itertools.combinations(removed,k):
    q=b[~b['sample'].isin(subset)].sort_values('sample')
    rho,p=stats.spearmanr(q.sender_mean_CP10k,q.mean_auc)
    rows.append(dict(filter=arm,scope=scope,seed=seed,tf=tf,threshold_before=before,threshold_after=after,all_excluded_at_transition=';'.join(removed),removed_subset=';'.join(subset) or 'NONE',n=len(q),n_control=int(q.group.eq('Control').sum()),n_case=int(q.group.eq('CWP').sum()),sample_ids=';'.join(q['sample']),rho=rho,p=p,df=len(q)-2,p_method='Spearman_t_approximation_exploratory',selection='All subsets of subjects removed by each predeclared threshold transition; no outcome-based omission'))
r=pd.DataFrame(rows)
keys=['filter','scope','seed','tf','sample_ids'];u=r.drop_duplicates(keys).copy();u['FDR_exploratory']=np.nan
v=u.p.dropna().index;pv=u.loc[v].p.to_numpy();order=np.argsort(pv);q=np.empty(len(pv));q[order]=np.minimum.accumulate((pv[order]*len(pv)/np.arange(1,len(pv)+1))[::-1])[::-1];u.loc[v,'FDR_exploratory']=np.minimum(q,1)
r=r.merge(u[keys+['FDR_exploratory']],on=keys,validate='many_to_one')
r.to_csv(ROOT/'results/filter_transition_subset_removals.tsv',sep='\t',index=False)
print('TRANSITION_SUBSETS_PASS',len(r))
print(r[(r['filter']=='All')&(r.scope=='CD4-IL7R')&(r.tf=='IKZF2')&(r.threshold_before==10)][['seed','removed_subset','n','rho','p']].to_string(index=False))
