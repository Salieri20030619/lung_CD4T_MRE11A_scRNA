from pathlib import Path
import pandas as pd,numpy as np,json
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports');O=P/'front_half_20261004_v1'
n=pd.read_csv(O/'results/receiver_audit_inputs.tsv',sep='\t')
o=pd.read_csv(P/'missing_links_20261003_v1/results/receiver_paired_modalities.tsv',sep='\t');o=o[o.species.eq('Human')]
d=n.merge(o,on=['species','scope','sample'],suffixes=('_new','_old'),validate='one_to_one');rows=[]
for a,b in [('State_Treg','TregNoTF_NoRegUnion'),('State_Activation','Activation_NoRegUnion')]:
 v=d[a]-d[b+'_old'];rows.append(dict(new_panel=a,old_panel=b,n=len(v),max_abs_error=float(np.max(np.abs(v))),identical=np.allclose(d[a],d[b+'_old'],rtol=1e-12,atol=1e-12)))
pd.DataFrame(rows).to_csv(O/'results/duplicate_score_audit.tsv',sep='\t',index=False)
print(json.dumps(rows))
