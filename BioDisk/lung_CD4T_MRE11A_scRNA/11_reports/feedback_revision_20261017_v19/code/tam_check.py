from pathlib import Path
import pandas as pd,numpy as np,scipy.sparse as sp,json
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports');B=P/'three_cohort_four_axis_20261017_v16';O=P/'feedback_revision_20261017_v19';rows=[]
for ds,core in [('GSE314943','Core'),('GSE274372','TLike')]:
 m=pd.read_csv(B/'inputs'/f'{ds}__{core}.cells.tsv.gz',sep='\t');g=pd.read_csv(B/'inputs'/f'{ds}__{core}.genes.txt',header=None)[0].tolist();idx={v:i for i,v in enumerate(g)};x=sp.load_npz(B/'inputs'/f'{ds}__{core}.cellcounts.npz');tidx=[idx[v] for v in ['CD3D','CD3E','CD3G','TRAC','CD247'] if v in idx];tp=np.asarray(x[:,tidx].sum(1)).ravel()>0
 for scope in ([core,'Naive','ActCD4','Treg'] if core=='Core' else [core]):
  take=np.ones(len(m),bool) if scope==core else m.subtype.eq(scope).to_numpy()
  for sample,z0 in m.loc[take].groupby('sample'):
   for variant in ['Raw','BG1']:
    z=z0 if variant=='Raw' else z0[z0.background_fraction<=.01];ids=z.index.to_numpy()
    for gene in ['AXL','TYRO3','MERTK']:
     present=gene in idx;v=x[ids,idx[gene]].toarray().ravel() if present else np.zeros(len(z));pos=v>0
     rows.append(dict(dataset=ds,scope=scope,sample=sample,group=z0.group.iloc[0],variant=variant,gene=gene,cells=len(z),feature_present=present,UMI=v.sum(),positive_cells=int(pos.sum()),fraction=float(pos.mean()) if len(z) else np.nan,Tmarker_positive=int((pos&tp[ids]).sum()),multiUMI=int((v>=2).sum()),threshold_eligible=bool(len(z)>=10 and present),sample_expression_pass=bool(len(z)>=10 and present and pos.mean()>=.05)))
d=pd.DataFrame(rows);d.to_csv(O/'tables/TAM_receiver_sample_coverage.tsv',sep='\t',index=False);qual=[]
for (ds,scope,gr,variant,gene),z in d.groupby(['dataset','scope','group','variant','gene']):
 n=int(z.threshold_eligible.sum());passed=int(z.sample_expression_pass.sum());qual.append(dict(dataset=ds,scope=scope,group=gr,variant=variant,gene=gene,available_samples=len(z),eligible_samples=n,passing_samples=passed,group_expression_pass=bool(n and passed/n>=.5)))
q=pd.DataFrame(qual);q.to_csv(O/'tables/TAM_receiver_group_qualification.tsv',sep='\t',index=False)
print('ACTCD4_RAW',d[(d.dataset=='GSE314943')&(d.scope=='ActCD4')&(d.variant=='Raw')][['sample','gene','cells','positive_cells','fraction']].to_string(index=False));print('GROUPS',q[(q.dataset=='GSE314943')&(q.variant=='Raw')].to_string(index=False))