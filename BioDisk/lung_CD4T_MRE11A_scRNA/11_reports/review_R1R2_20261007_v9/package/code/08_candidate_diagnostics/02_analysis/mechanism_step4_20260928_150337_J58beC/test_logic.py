import importlib.util,tempfile
from pathlib import Path
import numpy as np,pandas as pd
from scipy import sparse
from scipy.io import mmwrite
sp=importlib.util.spec_from_file_location('audit',Path(__file__).with_name('analyze.py'));a=importlib.util.module_from_spec(sp);sp.loader.exec_module(a)
with tempfile.TemporaryDirectory() as td:
 p=Path(td)
 # Same marginal receptor prevalence but no same-cell co-expression for sample A.
 x=sparse.csr_matrix([[1,0,2,1,0],[0,1,0,0,1],[1,1,1,0,0],[0,0,1,0,0]],dtype=int)
 mmwrite(p/'x.mtx',x.T)
 pd.DataFrame({'symbol':['Il2rb','Il2rg','Il6ra','Il6r','Il1r2']}).to_csv(p/'g.tsv',sep='\t',index=False)
 m=pd.DataFrame({'cell':['a','b','c','d'],'sample':['A','A','B','B'],'group':['Saline','Saline','Silica','Silica'],'subtype':['ActCD4']*4,'library':np.asarray(x.sum(axis=1)).ravel()})
 m.to_csv(p/'m.tsv',sep='\t',index=False)
 y,g,md,mp=a.load(p/'x.mtx',p/'g.tsv',p/'m.tsv','symbol')
 assert len(g)==4 and 'IL6R' in g and 'IL6RA' not in g
 assert np.allclose(y[:,list(g).index('IL6R')].toarray().ravel(),[3,0,1,1])
 co,ss=a.audit_scope(y,g,md,'test','ActCD4','All','CD4');d=pd.DataFrame(co)
 assert d[(d.axis=='IL15_TP')&(d['sample']=='A')].codetected.iloc[0]==0
 assert d[(d.axis=='IL15_TP')&(d['sample']=='B')].codetected.iloc[0]==1
 assert not d[d.axis=='IL18'].all_features_present.any() and d[d.axis=='IL18'].pct.isna().all()
 r=a.partial([1,3,2,4,6,5],[1,3,2,4,6,5],['A']*3+['B']*3);assert abs(r-1)<1e-10
 assert np.isnan(a.partial([1]*6,range(6),['A']*3+['B']*3))
 assert np.isnan(a.partial(range(6),range(6),['A']+['B']*5))
print('ANALYSIS_LOGIC_PASS alias_merge, complex_codetection, missing_vs_zero, sample_correlation')

