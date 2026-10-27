from pathlib import Path
import pandas as pd,numpy as np,itertools,json
from scipy.stats import rankdata
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/four_axis_screen_20261017_v15');S=['Saline_1','Saline_2','Saline_3','Silica_1','Silica_2','Silica_3'];rd=lambda f:pd.read_csv(f,sep='\t');identity=rd(O/'06_novelty/DUSP5_link_identity_counts.tsv');links=rd(O/'05_intersection/DUSP5_Treg_primary_links.tsv');y=identity[(identity.scope=='Treg')&(identity.gene=='DUSP5')].set_index('sample').loc[S].CP10k.to_numpy();rows=[]
for l in links.itertuples():
 x=identity[(identity.scope==l.sender)&(identity.gene==l.ligand)].set_index('sample').loc[S].CP10k.to_numpy()
 for adj in [False,True]:
  xx=rankdata(x);yy=rankdata(y)
  for k in ([np.arange(3),np.arange(3,6)] if adj else [np.arange(6)]):xx[k]-=xx[k].mean();yy[k]-=yy[k].mean()
  den=np.linalg.norm(xx)*np.linalg.norm(yy);rho=xx@yy/den;perms=(u+v for u in itertools.permutations(range(3)) for v in itertools.permutations(range(3,6))) if adj else itertools.permutations(range(6));null=[xx[list(p)]@yy/den for p in perms];p=sum(abs(z)>=abs(rho)-1e-12 for z in null)/len(null);refp=l.adj_P if adj else l.P;refr=l.adj_rho if adj else l.rho;assert abs(p-refp)<1e-10 and abs(rho-refr)<1e-10;rows.append(dict(ligand=l.ligand,sender=l.sender,scope='Treg',gene='DUSP5',model='Adjusted' if adj else 'Pooled',rho=rho,P=p,status='PASS'))
pd.DataFrame(rows).to_csv(O/'provenance/focal_source_exact_audit.tsv',sep='\t',index=False);print('SOURCE_EXACT_PASS',len(rows))