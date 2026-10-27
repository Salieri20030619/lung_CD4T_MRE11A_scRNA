from pathlib import Path
import pandas as pd,numpy as np,json
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/death_priority_20261018_v2')
d=pd.read_csv(O/'tables/module_member_sensitivity.tsv',sep='\t');rows=[]
for fam,z in d.groupby('family'):
 z=z[z.P.notna()];p=z.P.to_numpy();n=len(p);ix=np.argsort(p);v=np.minimum.accumulate((p[ix]*n/np.arange(1,n+1))[::-1])[::-1];q=np.empty(n);q[ix]=np.minimum(v,1);err=float(np.max(abs(q-z.FDR)));assert err<1e-12;assert z.n.eq(10).all();rows.append(dict(family=fam,tests=n,FDR_error=err,passed=True))
pd.DataFrame(rows).to_csv(O/'tables/module_member_QA.tsv',sep='\t',index=False)
(O/'logs/module_member_QA_complete.json').write_text(json.dumps({'status':'success','families':len(rows),'tests':len(d)},indent=2),encoding='utf-8')
print('MEMBER_QA',len(rows),len(d))
