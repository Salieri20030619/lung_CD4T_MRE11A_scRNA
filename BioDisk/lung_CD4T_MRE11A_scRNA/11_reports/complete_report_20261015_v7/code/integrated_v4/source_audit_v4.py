from pathlib import Path
import sys,pandas as pd,numpy as np,json
from scipy.io import mmread
O=Path(sys.argv[1]);P=O.parents[1];M=P/'02_analysis/mechanism_step1_20260927_152218_s4V5eF';act=pd.read_csv(O/'02_sender/ligand_activity_all_states.tsv',sep='\t');focal=sorted(set(act.loc[(act.direction=='Up')&(act['rank']<=5),'test_ligand'])|{'IL15','IL18','TNF'})
md=pd.read_csv(M/'comm_input/Mac_cells.tsv',sep='\t');gn=pd.read_csv(M/'comm_input/Mac_genes.tsv',sep='\t');gi=pd.Index(gn.canonical);x=mmread(M/'comm_input/Mac_counts.mtx.gz').T.tocsr();panels={'myeloid':'LST1 TYROBP FCER1G CSF1R ADGRE1 C1QA C1QB C1QC MRC1'.split(),'T_lineage':'CD3D CD3E TRAC LCK'.split(),'epithelial':'SFTPC EPCAM KRT8 KRT18 SCGB1A1'.split()};flags={}
for name,gs in panels.items():
 ix=gi.get_indexer(gs);ix=ix[ix>=0];flags[name]=np.asarray((x[:,ix]>0).sum(axis=1)).ravel()>=2
rows=[]
for source in ['AMlike','AllMac']:
 for sample,b in md.groupby('sample'):
  take=md['sample'].eq(sample).to_numpy()&(np.ones(len(md),bool) if source=='AllMac' else md.subtype.eq(source).to_numpy())
  for g in focal:
   if g not in gi:continue
   pos=(x[:,gi.get_loc(g)].toarray().ravel()>0)&take;n=int(pos.sum());r=dict(sender=source,sample=sample,group=b.group.iloc[0],ligand=g,source_cells=int(take.sum()),positive_cells=n,pct_positive=100*n/take.sum())
   for name,f in flags.items():r[name+'_positive_cells']=int((pos&f).sum());r[name+'_fraction_among_positive']=float(f[pos].mean()) if n else np.nan
   rows.append(r)
a=pd.DataFrame(rows);a.to_csv(O/'02_sender/lead_ligand_source_identity_audit.tsv',sep='\t',index=False,na_rep='NA')
# Source-specific ranks retain the same activity score; coverage restriction is a separate sensitivity.
el=pd.read_csv(O/'02_sender/ligand_scope_eligibility.tsv',sep='\t');rob=pd.read_csv(O/'02_sender/LR_legacy10pct_robustness.tsv.gz',sep='\t').groupby(['sender','scope','ligand']).legacy10pct_supported.any().reset_index();ranks=[]
for (source,scope),v in el.groupby(['sender','scope']):
 z=act[act.scope==scope].merge(v[['ligand','expression_supported','supported_receptors']],left_on='test_ligand',right_on='ligand').merge(rob[(rob.sender==source)&(rob.scope==scope)][['ligand','legacy10pct_supported']],on='ligand');z=z[z.expression_supported].copy();z['sender']=source
 z['source_specific_rank']=z.groupby('direction').aupr_corrected.rank(method='min',ascending=False);z['coverage10_rank']=np.nan
 for direction in ['Up','Down']:
  ix=z.index[(z.direction==direction)&z.legacy10pct_supported];z.loc[ix,'coverage10_rank']=z.loc[ix,'aupr_corrected'].rank(method='min',ascending=False)
 ranks.append(z)
ranks=pd.concat(ranks);ranks.to_csv(O/'02_sender/source_specific_activity_ranks.tsv',sep='\t',index=False,na_rep='NA')
print('COVERAGE_RANKS');print(ranks[(ranks.sender=='AMlike')&(ranks.test_ligand.isin(['IL15','IL18','TNF']))&(ranks.direction=='Up')][['scope','test_ligand','rank','source_specific_rank','legacy10pct_supported','coverage10_rank']].to_string(index=False));print('IFNG_SOURCE_CASES');print(a[(a.ligand=='IFNG')&(a.group=='Silica')].to_string(index=False))