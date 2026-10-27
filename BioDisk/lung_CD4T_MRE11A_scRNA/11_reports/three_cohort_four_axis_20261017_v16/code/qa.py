from pathlib import Path
import pandas as pd,numpy as np,json
from scipy.stats import false_discovery_control
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/three_cohort_four_axis_20261017_v16');d=pd.read_csv(O/'tables/DE_all_models.tsv',sep='\t',low_memory=False);checks=[]
for fam,z in d.groupby('family'):
 assert z.gene.is_unique;assert z.PValue.between(0,1).all();fdr=false_discovery_control(z.PValue.to_numpy());delta=np.max(np.abs(fdr-z.FDR));assert delta<1e-10;checks.append(dict(check='BH_family_recalculation',family=fam,n=len(z),max_abs_error=delta,status='PASS'))
inv=pd.read_csv(O/'tables/sample_inventory.tsv',sep='\t');assert not inv.duplicated(['dataset','scope','sample']).any();assert inv[(inv.dataset=='GSE274372')&(inv.scope=='TLike')].cells.sum()==494;assert len(inv[(inv.dataset=='GSE183682')&(inv.scope=='Core')])==4
for f in (O/'inputs').glob('*.counts.tsv.gz'):
 x=pd.read_csv(f,sep='\t',index_col=0);assert x.index.is_unique;assert (x.values>=0).all() and np.allclose(x.values,np.rint(x.values));stem=f.name.replace('.counts.tsv.gz','');m=pd.read_csv(O/'inputs'/f'{stem}.samples.tsv',sep='\t');assert list(x.columns)==list(m['sample']);assert (x.sum().values<=m.library.values+1e-6).all();checks.append(dict(check='pseudobulk_integrity',family=stem,n=x.shape[1],max_abs_error=0,status='PASS'))
pd.DataFrame(checks).to_csv(O/'tables/QA_statistical_integrity.tsv',sep='\t',index=False);print('PASS',len(checks),'families',d.family.nunique())