from pathlib import Path
import os,ast,json
os.environ['OPENBLAS_NUM_THREADS']='1'
import pandas as pd,numpy as np
from scipy import stats
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/human_alignment_20261005_v1';R=O/'results'
tree=ast.parse((O/'code/align_compute.py').read_text());ns=globals();exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['ols','bh','wr']],type_ignores=[]),'<reused finite functions>','exec'),ns)
d=pd.read_csv(R/'TF_readout_donors.tsv.gz',sep='\t');des=pd.read_csv(R/'stage1_donor_design.tsv',sep='\t');q=d.merge(des[['sample','group','series']],on='sample',suffixes=('_readout','_design'));bad=q[(q.group_readout!=q.group_design)|(q.series_readout!=q.series_design)];print('METADATA_MISMATCH',len(bad));print(bad[['sample','scope','readout','group_readout','group_design','series_readout','series_design']].drop_duplicates().to_string(index=False))
wr(bad,'TF_metadata_mismatch_audit.tsv.gz')
if len(bad):raise SystemExit('METADATA_MISMATCH_STOP')
rows=[];loo=[]
for (scope,tf),g in d[(d.readout=='RNA')&(d['filter']=='All')].groupby(['scope','tf']):
 for roster,cut in [('PositiveCells',1),('Cells20',20)]:
  v=g[g.cells>=cut].copy();v['cells']=np.maximum(v.cells,20)
  for model in ['Base','Source','SourceDep']:
   tag=dict(scope=scope,tf=tf,roster=roster,model=model);a=ols(v,model);rows.append(dict(tag,**a))
   if scope=='CD4-GZMA':
    for s in v['sample']:loo.append(dict(tag,omitted=s,full_effect=a['effect'],**ols(v[v['sample']!=s],model)))
rows=bh(pd.DataFrame(rows),['roster','model']);wr(rows,'TF_RNA_roster_models.tsv');wr(loo,'GZMA_TF_RNA_donor_influence.tsv')
print('ROSTER_GZMA',rows[rows.scope=='CD4-GZMA'].to_string(index=False))
print('LOO_GZMA',pd.DataFrame(loo).groupby(['tf','roster','model']).agg(computed=('P','count'),min_effect=('effect','min'),max_effect=('effect','max'),max_P=('P','max')).to_string())
de=pd.read_csv(R/'sender_disease_fixed_candidates.tsv',sep='\t');lig='IL15 TNF CLSTN1 ADAM15 PLAU DSC2 CXCL2 ICAM2 VTN VEGFB GSTP1 RTN4 IL15RA'.split();grid=pd.MultiIndex.from_product([lig,['MacAM','MonoMac','MacWide'],['Base','Source']],names=['gene','sender','model']).to_frame(index=False).merge(de,on=['gene','sender','model'],how='left');grid['status']=np.where(grid.PValue.notna(),'TESTED','FILTERED_NOT_TESTED');wr(grid,'sender_disease_candidates_with_status.tsv')
