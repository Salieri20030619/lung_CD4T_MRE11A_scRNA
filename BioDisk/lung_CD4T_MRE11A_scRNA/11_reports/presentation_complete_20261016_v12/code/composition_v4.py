from pathlib import Path
import sys,json,itertools
import numpy as np,pandas as pd
from scipy.stats import rankdata,false_discovery_control
O=Path(sys.argv[1]);samples=['Saline_1','Saline_2','Saline_3','Silica_1','Silica_2','Silica_3'];states=['Naive','ActCD4','Treg'];mods=['PyroExec','Inflam','ApoExec','OXPHOS','MitoISR','MitoROS','Mitophagy']
a=pd.read_csv(O/'03_receiver/TF_all_regulon_sample_activity.tsv.gz',sep='\t');p=pd.read_csv(O/'01_phenotype/sample_program_scores.tsv',sep='\t');p=p[p['filter']=='All'];r=pd.read_csv(O/'inputs/sample_scope_roster.tsv',sep='\t');r=r[r.lineage=='CD4'];ct=r[r.scope!='Core'].pivot(index='sample',columns='scope',values='cells').loc[samples,states];weights=ct.div(ct.sum(axis=1),axis=0).mean(axis=0)
pd.DataFrame({'state':states,'fixed_weight':weights.values}).to_csv(O/'06_support/Core_fixed_composition_weights.tsv',sep='\t',index=False)
X=None;Y=np.zeros((6,len(mods)))
for state in states:
 b=a[(a.scope==state)&(a.seed==17)].pivot(index='sample',columns='regulon',values='mean_auc').loc[samples];X=b*weights[state] if X is None else X+b*weights[state];Y+=p[p.scope==state].set_index('sample').loc[samples,mods].to_numpy()*weights[state]
X.to_csv(O/'06_support/Core_fixed_composition_TF_inputs.tsv',sep='\t',index_label='sample');pd.DataFrame(Y,index=samples,columns=mods).to_csv(O/'06_support/Core_fixed_composition_program_inputs.tsv',sep='\t',index_label='sample')
res=[]
for adjust in [False,True]:
 xx=rankdata(X.to_numpy(),axis=0);yy=rankdata(Y,axis=0)
 if adjust:
  for sl in [slice(0,3),slice(3,6)]:xx[sl]-=xx[sl].mean(axis=0);yy[sl]-=yy[sl].mean(axis=0)
  perm=np.array([list(i)+list(j) for i in itertools.permutations(range(3)) for j in itertools.permutations(range(3,6))])
 else:xx-=xx.mean(axis=0);yy-=yy.mean(axis=0);perm=np.array(list(itertools.permutations(range(6))))
 den=np.sqrt((xx*xx).sum(axis=0))[:,None]*np.sqrt((yy*yy).sum(axis=0))[None,:];den[den<1e-12]=np.nan;rho=xx.T@yy/den;null=np.einsum('psa,sb->pab',xx[perm],yy,optimize=True)/den;pv=(np.abs(null)>=np.abs(rho)[None]-1e-12).mean(axis=0);pv[~np.isfinite(rho)]=np.nan;q=np.full(pv.size,np.nan);ok=np.isfinite(pv.ravel());q[ok]=false_discovery_control(pv.ravel()[ok]);frame=pd.DataFrame({'regulon':np.repeat(X.columns,len(mods)),'module':np.tile(mods,len(X.columns)),'rho':rho.ravel(),'P':pv.ravel(),'BH':q,'model':'ConditionAdjusted' if adjust else 'Pooled','family':'Core_fixedcomposition_allTF7_'+str(adjust)});res.append(frame)
res=pd.concat(res);res.to_csv(O/'06_support/Core_fixed_composition_TF_program.tsv',sep='\t',index=False,na_rep='NA');f=res[res.regulon.str.split('(').str[0].str.upper().isin(['IKZF2','STAT1'])];f.to_csv(O/'06_support/Core_fixed_composition_focal.tsv',sep='\t',index=False,na_rep='NA');print(f[f.module.isin(['ApoExec','MitoROS','PyroExec'])].to_string(index=False))
# Prior edges retain signs and confidence, never a claim of measured signaling.
net=pd.read_csv(O/'inputs/TF_target_prior.tsv',sep='\t');member=pd.read_csv(O.parents[1]/'02_analysis/mechanism_step1_20260927_152218_s4V5eF/tables/module_definition.tsv',sep='\t');b=net[net.tf.isin(['IKZF2','STAT1'])].merge(member,on=None,left_on='target',right_on='gene',how='left');b.to_csv(O/'03_receiver/focal_TF_prior_targets_and_modules.tsv',sep='\t',index=False,na_rep='NA')