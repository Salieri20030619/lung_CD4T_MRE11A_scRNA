from pathlib import Path
import json,datetime,hashlib
import numpy as np,pandas as pd
from scipy import stats
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/breakpoint_TNF_IL18_20261004_v1');R=O/'results'
assert not (O/'provenance/missing_LOO_complete.json').exists()
path=R/'aligned_chain_models.tsv';m=pd.read_csv(path,sep='\t');inp=pd.read_csv(R/'aligned_model_inputs.tsv.gz',sep='\t');q=m[m.LOO_min.isna()]
assert len(q)==12 and q.signature.nunique()==12 and q.calculation.eq('reused_historical').all()
before=m[['spec_id','rho','P','FDR']].copy();(O/'provenance/aligned_models_before_LOO.tsv').write_bytes(path.read_bytes());rows=[]
for i,v in q.iterrows():
 z=inp[inp.spec_id.eq(v.spec_id)].sort_values('sample');covs=v.covariates.split('+');vals=[]
 for omit in z['sample']:
  d=z[z['sample'].ne(omit)];C=[np.ones(len(d))]
  for c in covs:
   if c=='G':C.append(d.group.eq('CWP').to_numpy(float))
   elif c=='S':C.extend([d.series.eq(s).to_numpy(float) for s in sorted(d.series.unique())[1:]])
   else:
    x=stats.rankdata(d[c]);C.append((x-x.mean())/x.std() if x.std()>0 else np.zeros(len(d)))
  C=np.column_stack(C);rank=np.linalg.matrix_rank(C);df=len(d)-rank-1;rho=np.nan;status='FIT'
  if rank<C.shape[1] or df<3:status='INSUFFICIENT_RANK_OR_DF'
  else:
   x=stats.rankdata(d.x_value);y=stats.rankdata(d.y_value);x-=C@np.linalg.lstsq(C,x,rcond=None)[0];y-=C@np.linalg.lstsq(C,y,rcond=None)[0]
   if min(np.linalg.norm(x),np.linalg.norm(y))<1e-10:status='CONSTANT_RESIDUAL'
   else:rho=float(np.corrcoef(x,y)[0,1]);vals.append(rho)
  rows.append(dict(spec_id=v.spec_id,signature=v.signature,omitted_sample=omit,n=len(d),df=df,rho=rho,status=status))
 m.loc[i,['LOO_valid','LOO_min','LOO_max','LOO_same_direction']]=[len(vals),min(vals) if vals else np.nan,max(vals) if vals else np.nan,bool(vals and all(np.sign(x)==np.sign(v.rho) for x in vals))]
 m.loc[i,'LOO_detail']='new_LOO_only_for_missing_historical_diagnostic'
assert before.equals(m[['spec_id','rho','P','FDR']])
m.to_csv(path,sep='\t',index=False,na_rep='NA');pd.DataFrame(rows).to_csv(R/'historical_missing_LOO_completed.tsv',sep='\t',index=False,na_rep='NA')
# Refresh the view table; no model is fit here.
m[(m.policy=='aligned_steps')&(m.covariates=='G+S+receiver_library+sender_library')].to_csv(R/'common_final_chain_edges.tsv',sep='\t',index=False,na_rep='NA')
receipt=dict(time=datetime.datetime.now().isoformat(),historical_models_with_LOO_added=12,omission_diagnostics=len(rows),main_rho_P_FDR_unchanged=True,remaining_missing_LOO=int(m.LOO_min.isna().sum()))
(O/'provenance/missing_LOO_complete.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt,indent=2))
print(m[m.spec_id.isin(q.spec_id)][['chain','modality','policy','covariates','n','rho','P','FDR','LOO_valid','LOO_min','LOO_max']].to_string(index=False))
