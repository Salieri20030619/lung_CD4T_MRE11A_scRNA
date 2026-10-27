from pathlib import Path
import sys,pandas as pd,numpy as np
O=Path(sys.argv[1]);a=pd.read_csv(O/'direction/mouse_candidate_cell_inputs.tsv.gz',sep='\t');out=[];coverage=[]
for (ds,day),z in a[a.scope=='ActCD4'].groupby(['dataset','day']):
 z=z.copy();v=np.log1p(z.library.to_numpy());w=z['sample'].map(1/z.groupby('sample').size()).to_numpy();ix=np.argsort(v,kind='stable');cuts=np.unique(np.interp([1/3,2/3],np.cumsum(w[ix])/w.sum(),v[ix]));z['stratum']=np.searchsorted(cuts,v,side='left');ct=pd.crosstab(z['sample'],z.stratum);common=ct.columns[(ct>0).all()];prop=ct[common].div(ct[common].sum(axis=1),axis=0);groups=z.groupby('sample').group.first();weights=(prop.loc[groups=='Saline'].mean()+prop.loc[groups=='Silica'].mean())/2
 for sample,b in z.groupby('sample'):coverage.append(dict(dataset=ds,day=day,sample=sample,group=b.group.iloc[0],mean_library=b.library.mean(),median_library=b.library.median(),cells=len(b),coverage=b.stratum.isin(common).mean(),common_bins=';'.join(map(str,common))))
 for g in ['AHNAK','MDFIC','TTC39C','ST6GALNAC3']:
  vals=[]
  for sample,b in z.groupby('sample'):
   mu=b.groupby('stratum')[g+'_CP10k'].mean();vals.append((b.group.iloc[0],sum(weights[j]*mu[j] for j in common)))
  vals=pd.DataFrame(vals,columns=['group','value']);out.append(dict(dataset=ds,day=day,gene=g,depth_standardized_difference=vals[vals.group=='Silica'].value.mean()-vals[vals.group=='Saline'].value.mean(),method='equal-library depth tertiles, common support, fixed common weights; diagnostic only'))
pd.DataFrame(out).to_csv(O/'direction/ActCD4_depth_standardization.tsv',sep='\t',index=False,encoding='utf-8');pd.DataFrame(coverage).to_csv(O/'direction/ActCD4_depth_coverage.tsv',sep='\t',index=False,encoding='utf-8')
print(pd.DataFrame(out).query("gene in ['MDFIC','TTC39C']").to_string(index=False));print(pd.DataFrame(coverage).to_string(index=False))