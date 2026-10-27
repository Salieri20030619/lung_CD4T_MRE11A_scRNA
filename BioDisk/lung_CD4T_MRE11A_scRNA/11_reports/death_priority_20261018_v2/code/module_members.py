import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
from pathlib import Path
import numpy as np,pandas as pd,json
from scipy import sparse
import statslib as S
O=S.O;B=O.parent/'three_cohort_four_axis_20261017_v16'
plan=[('HIF1A','Naive','PyroExec'),('ARL6IP1','Naive','MitoROS'),('IL2RA','Treg','OXPHOS'),('IL2RA','Treg','Glycolysis'),('IL2RA','Treg','ApoExec')]
(O/'inputs/module_member_sensitivity_protocol.json').write_text(json.dumps({'status':'frozen_before_member_results','trigger':'posthoc assessment of whether five observed primary-candidate program relationships depend on one member; not a new discovery gate','pairs':plan,'models':['JointStratum','JointStratumCondition'],'variants':['Raw','BG1'],'method':'leave-one-member mean logCP10k; candidates already excluded; one family per gene-scope-module-variant-model; primary P/FDR unchanged','statistics':'same independent10 samples and exact permutation rules; member omission is not independent validation'},indent=2),encoding='utf-8')
mods=pd.read_csv(O/'tables/module_members.tsv',sep='\t');sets={k:list(z.loc[z.retained,'gene']) for k,z in mods.groupby('module')}
data={}
for ds in ['GSE314943','GSE183682']:
 stem=ds+'__Core';x=sparse.load_npz(B/'inputs'/(stem+'.cellcounts.npz'));m=pd.read_csv(B/'inputs'/(stem+'.cells.tsv.gz'),sep='\t');g=pd.read_csv(B/'inputs'/(stem+'.genes.txt'),header=None)[0].tolist();data[ds]=(x,m,{v:i for i,v in enumerate(g)})
for gene,scope,mod in plan:
 members=sorted(set(sets[mod]));rows=[]
 for var in ['Raw','BG1']:
  rows=[]
  for ds,(x,m,ix) in data.items():
   take=m.subtype.eq(scope).to_numpy()
   if var=='BG1':take&=m.background_fraction.to_numpy()<=.01
   for sam,mm in m.loc[take].groupby('sample'):
    ids=mm.index.to_numpy();counts=np.asarray(x[ids].sum(0)).ravel();log=np.log1p(counts/mm.library.sum()*1e4);vals=log[[ix[g] for g in members]]
    r=dict(dataset=ds,sample=sam,group=mm.group.iloc[0],day=int(mm.day.iloc[0]),cells=len(mm),G=log[ix[gene]],full_score=float(vals.mean()))
    for j,g in enumerate(members):r['omit_'+g]=float((vals.sum()-vals[j])/(len(vals)-1))
    rows.append(r)
  d=pd.DataFrame(rows).sort_values(['dataset','sample']).reset_index(drop=True);assert len(d)==10
  for model in ['JointStratum','JointStratumCondition']:
   for member in members:S.assoc(d,'G','omit_'+member,'Joint',var,'member_sensitivity_'+gene+'_'+scope+'_'+mod,scope,model,extra={'gene':gene,'program':mod,'omitted_member':member,'module_size':len(members)})
a=S.finalize(S.results,'module_member_sensitivity.tsv');out=[]
for (gene,scope,mod,var,model),z in a.groupby(['gene','scope','program','variant','model']):
 out.append(dict(gene=gene,scope=scope,program=mod,variant=var,model=model,members=len(z),rho_min=z.rho.min(),rho_max=z.rho.max(),P_min=z.P.min(),P_max=z.P.max(),nominal_omissions=int((z.P<.05).sum()),positive_omissions=int((z.rho>0).sum()),most_attenuating_member=z.loc[z.rho.idxmin(),'omitted_member']))
s=pd.DataFrame(out);s.to_csv(O/'tables/module_member_sensitivity_summary.tsv',sep='\t',index=False)
(O/'logs/module_member_sensitivity_complete.json').write_text(json.dumps({'status':'success','tests':len(a),'relationships':len(plan)},indent=2),encoding='utf-8')
print(s[(s.variant=='Raw')&(s.model=='JointStratum')].to_string(index=False))
