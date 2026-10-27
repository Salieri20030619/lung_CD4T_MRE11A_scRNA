from pathlib import Path
import pandas as pd,numpy as np,json,re
from scipy.stats import false_discovery_control
O=Path(__file__).parent;H=O.parent;T=O/'tables';U=O/'server_tables';V=H/'lung_append_analysis_20261017_v18'
mk=pd.read_csv(U/'candidate_positive_identity.tsv',sep='\t');z=mk.groupby('gene')[['positive','positive_Tmarker','positive_BG1']].sum();z['Tmarker_fraction']=z.positive_Tmarker/z.positive;z['BG1_retained_fraction']=z.positive_BG1/z.positive;z.to_csv(T/'identity_summary.tsv',sep='\t')
# Annotate gene content in displayed enrichment; do not reinterpret named diseases literally.
rows=[]
for num in ['F07','F08','F09']:
 d=pd.read_csv(O/'plot_data'/(num+'.tsv'),sep='\t');col='hit_genes' if 'hit_genes' in d else 'leadingEdge'
 for r in d.itertuples():
  genes=str(getattr(r,col)).split(';');igs=[g for g in genes if re.match(r'^(IGH|IGK|IGL)',g)];ribo=[g for g in genes if re.match(r'^(RPL|RPS|MRPL|MRPS)',g)];rows.append(dict(figure=num,scope=r.scope,database=r.database,pathway=r.pathway,description=r.description,member_type=col,members=len(genes),immunoglobulin_members=';'.join(igs),ribosomal_members=';'.join(ribo),immunoglobulin_fraction=len(igs)/len(genes),ribosomal_fraction=len(ribo)/len(genes)))
pd.DataFrame(rows).to_csv(T/'displayed_pathway_member_context.tsv',sep='\t',index=False)
# Recheck whole trajectory program BH families from saved full tables, not display subset.
tr=O/'sources/11_reports/customer_sequence_20261015_v2/03A_state_pseudotime';d=pd.read_csv(tr/'program_trajectory_tests.tsv',sep='\t');q=[]
for f,z in d.groupby('family'):
 e=np.max(abs(false_discovery_control(z.P)-z.BH));assert e<1e-10;q.append(dict(check='trajectory_program_full_family_BH',family=f,error=e))
pd.DataFrame(q).to_csv(T/'trajectory_BH_QA.tsv',sep='\t',index=False)
g=pd.read_csv(U/'GSEA_all.tsv',sep='\t');g['P_at_floor']=g.pval<=1e-8;g['precision_status']=np.where(g.P_at_floor,'P<=1e-8 (algorithm floor)',np.where(g.log2err.isna(),'log2err unavailable; possible conservative P estimation','finite log2err'));g.to_csv(T/'GSEA_with_precision_status.tsv',sep='\t',index=False)
# Candidate selection membership unchanged; newly added pathway evidence does not feed back into selection.
a=pd.read_csv(H/'lung_feedback_revision_20261017_v19/tables/R05_current_all90.tsv',sep='\t');assert [int(a.A_status_current.str.startswith(x).sum()) for x in ['A1','A2','A3']]==[52,37,1]
print('Identity',z.shape,'pathway_context',len(rows),'GSEA_floor',int(g.P_at_floor.sum()),'precision_NA_notfloor',int((g.log2err.isna()&~g.P_at_floor).sum()))