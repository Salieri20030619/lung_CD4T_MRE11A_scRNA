from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,anndata as ad
from scipy.sparse.csgraph import connected_components
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');V=P/'11_reports/death_priority_20261018_v31';W=P/'11_reports/death_priority_20261018_v2'
O=Path('/mnt/BioData/bio/wangmingyuan/lung_CD4T_MRE11A_scRNA/20_v11_targeted_review_20261019');O.mkdir(exist_ok=True)
def wr(d,n):pd.DataFrame(d).to_csv(O/n,sep='\t',index=False)
before={};qs=[];bins=[];leverage=[];graphs=[];hvg=[];roots=[]
qa=pd.read_csv(V/'tables/graph_QA.tsv',sep='\t')
for track in ['Naive_ActCD4','Treg']:
 f=V/'objects'/f'{track}_geometry.h5ad';before[str(f)]=hashlib.sha256(f.read_bytes()).hexdigest()
 a=ad.read_h5ad(f);d=pd.read_csv(V/'tables'/f'{track}_cell_trajectory.tsv.gz',sep='\t')
 assert list(a.obs_names)==d.uid.tolist()
 n,cl=connected_components(a.obsp['connectivities'],directed=False)
 for v in ['K15_R1','K15_R2','K30_R1','K30_R2']:
  rr=qa[(qa.track==track)&(qa.variant==v)].iloc[0];root=rr.root_cell
  roots.append(dict(track=track,variant=v,root=root,root_subtype=d.loc[d.uid==root,'subtype'].iloc[0],saved_object_iroot=str(a.obs_names[int(a.uns['iroot'])]),saved_object_matches_variant=np.allclose(a.obs.dpt_pseudotime.to_numpy(),d[v].to_numpy(),equal_nan=True)))
  for scope,z in d.groupby('subtype'):
   raw=z[v].to_numpy();finite=np.isfinite(raw);val=raw[finite];scaled=(val-val.min())/(val.max()-val.min())
   r=dict(track=track,scope=scope,variant=v,cells=len(z),finite=int(finite.sum()),nan=int(np.isnan(raw).sum()),inf=int(np.isinf(raw).sum()))
   for lab,xx in [('raw',val),('scaled',scaled)]:
    for q in [0,.01,.25,.5,.75,.9,.95,.99,1]:r[lab+'_q'+str(q)]=float(np.quantile(xx,q))
   qs.append(r)
   if v=='K15_R1':
    zz=z.loc[finite].copy();zz['scaled']=scaled;zz['bin']=np.minimum((scaled*4).astype(int),3)
    for sam,t in zz.groupby('sample_id'):
     for b in range(4):bins.append(dict(track=track,scope=scope,sample_id=sam,bin=b,cells=int((t.bin==b).sum())))
     tt=t.scaled.to_numpy();dep=np.log1p(t.library.to_numpy());C=np.column_stack([np.ones(len(t)),dep]);res=tt-C@np.linalg.lstsq(C,tt,rcond=None)[0];hl=res**2/sum(res**2);i=int(hl.argmax())
     leverage.append(dict(track=track,scope=scope,sample_id=sam,cells=len(t),nonempty_bins=t.bin.nunique(),DPT_min=tt.min(),DPT_max=tt.max(),DPT_span=np.ptp(tt),max_partial_leverage=hl.max(),max_leverage_cell=t.uid.iloc[i],max_leverage_DPT=tt[i],top1pct_leverage_sum=np.sort(hl)[-max(1,int(np.ceil(len(hl)*.01))):].sum()))
 for comp in range(n):
  z=d[cl==comp];graphs.append(dict(track=track,k=15,component=comp,cells=len(z),samples=z.sample_id.nunique(),subtypes=';'.join(sorted(z.subtype.unique())),root1_in_component=qa[(qa.track==track)&(qa.variant=='K15_R1')].root_cell.iloc[0] in set(z.uid)))
 for g in ['HIF1A','IL2RA','ARL6IP1']:hvg.append(dict(track=track,gene=g,in_geometry=g in a.var_names,n_geometry_features=a.n_vars))
wr(qs,'DPT_quantiles.tsv');wr(bins,'DPT_donor_bin_coverage.tsv');wr(leverage,'DPT_partial_leverage.tsv');wr(graphs,'graph_components.tsv');wr(hvg,'candidate_geometry_membership.tsv');wr(roots,'root_and_saved_object_audit.tsv')
cur=pd.read_csv(V/'tables/module_trajectory_curves.tsv',sep='\t');ss=pd.DataFrame(bins)
z=cur[cur.module==cur.module.iloc[0]].merge(ss,on=['scope','sample_id','bin'],suffixes=('_old','_check'))
assert (z.cells_old==z.cells_check).all()
lr=pd.read_csv(W/'tables/joint_LR_eligibility.tsv',sep='\t');lr=lr[(lr.scope=='Treg')&lr.ligand_parts.eq('LGALS9')];wr(lr,'LGALS9_Treg_receptor_eligibility.tsv')
re=pd.read_csv(W/'tables/receiver_expression.tsv',sep='\t');se=pd.read_csv(W/'tables/sender_expression.tsv',sep='\t')
print('RECEIVER_COLUMNS',re.columns.tolist());print('SENDER_COLUMNS',se.columns.tolist())
genes=set(';'.join(lr.receptor_parts).split(';'));wr(re[(re.scope=='Treg')&re.gene.isin(genes)],'LGALS9_receptor_per_sample.tsv');wr(se[se.gene.eq('LGALS9')],'LGALS9_sender_per_sample.tsv')
for f,h in before.items():assert hashlib.sha256(Path(f).read_bytes()).hexdigest()==h
(O/'source_snapshot.json').write_text(json.dumps(before,indent=2),encoding='utf-8')
summary=dict(status='read_only_checks_complete',DPT_quantiles=qs,graph_components=graphs,candidate_geometry=hvg,source_hashes_unchanged=True)
(O/'audit_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print('GRAPH',graphs);print('HVG',hvg);print('BINS',ss.groupby('scope').apply(lambda x:x.groupby('sample_id').cells.apply(lambda a:int((a>0).sum())).to_dict(),include_groups=False).to_dict());print('LEVERAGE',pd.DataFrame(leverage).groupby('scope').max_partial_leverage.agg(['min','median','max']).to_dict());print('LR',lr.to_dict('records'))
