from pathlib import Path
import os,json,runpy
os.environ['OPENBLAS_NUM_THREADS']='1'
import pandas as pd,numpy as np
from scipy import stats
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/human_alignment_20261005_v1';R=O/'results';H=P/'11_reports/human_first_20261004_v1/results'
def rd(n):return pd.read_csv(R/n,sep='\t',low_memory=False)
def wr(d,n):pd.DataFrame(d).to_csv(R/n,sep='\t',index=False,na_rep='NA')
g=rd('alignment_gene_donors.tsv.gz');des=rd('stage1_donor_design.tsv');s=rd('IL15_sender_joint_coverage.tsv.gz').merge(des[['sample','group','series']],on='sample');r=rd('IL15_receiver_complex.tsv.gz');recgenes=g[g.gene.isin(['IL2RB','IL2RG'])].pivot(index=['scope','filter','sample'],columns='gene',values='pct').reset_index().rename(columns={'scope':'receiver'})
r=r.merge(recgenes,on=['receiver','filter','sample']);pairs=[]
for sender,dd in s.groupby('sender'):
 d=r.merge(dd,on='sample',suffixes=('_receiver','_sender'));d=d[(d.cells_receiver>=20)&(d.cells_sender>=20)];d['sender']=sender
 for cut in [.1,.05]:
  q=d.copy();q['gate']=cut;q['sender_marginal']=(q.IL15>=cut)&(q.IL15RA>=cut);q['receiver_marginal']=(q.IL2RB>=cut)&(q.IL2RG>=cut);q['both_marginal']=q.sender_marginal&q.receiver_marginal;q['both_joint']=(q.joint>=cut)&(q.pct>=cut);pairs.append(q)
pairs=pd.concat(pairs,ignore_index=True);wr(pairs,'IL15_matched_donor_coverage.tsv.gz');summary=pairs.groupby(['sender','receiver','filter','group','gate']).agg(n=('sample','size'),sender_marginal_n=('sender_marginal','sum'),receiver_marginal_n=('receiver_marginal','sum'),both_marginal_n=('both_marginal','sum'),both_joint_n=('both_joint','sum'),sender_joint_median=('joint','median'),receiver_joint_median=('pct','median')).reset_index();wr(summary,'IL15_matched_coverage_summary.tsv')
lt=rd('fixed_ligand_target_complete_members.tsv.gz');net=rd('existing_TF_edges_complete.tsv.gz');z=lt.merge(net[['filter','seed','tf','target','weight','retained_DropOverlap']],on=['filter','target'],suffixes=('_LT','_SCENIC'));z=z[z.target!=z.tf];wr(z,'ligand_disease_target_TF_edge_alignment.tsv.gz')
edgecount=net[net.target!=net.tf].groupby(['filter','seed','tf']).agg(edges=('target','size'),positive_edges=('weight',lambda x:sum(x>0))).reset_index();wr(edgecount,'existing_TF_edge_counts.tsv')
pm=rd('stage1_complete_program_members.tsv');nn=net[net.target!=net.tf].merge(pm[['program','gene']],left_on='target',right_on='gene');wr(nn,'TF_edges_to_existing_program_genes.tsv.gz')
# New observed disease target profiles for all fixed candidates and exact primary contexts.
ov=rd('TF_disease_target_overlap.tsv');q=ov[(ov['filter']=='All')&(ov.model=='Source')&(ov.direction=='Up')&(ov.test=='FDR05')&(ov.receiver=='CD4-GZMA')];print('GZMA_FDR_TF_OVERLAP',q[['tf','seed','variant','target_n','regulon_background_n','overlap_n','overlap_genes','P','FDR']].to_string(index=False))
d=rd('TF_disease_models.tsv');print('RNA_GZMA',d[(d.readout=='RNA')&(d['filter']=='All')&(d.scope=='CD4-GZMA')].to_string(index=False))
print('IL15_MATCHED_CWP',summary[(summary['filter']=='All')&(summary.group=='CWP')&(summary.gate==.1)].to_string(index=False))
print('TF_EDGES',edgecount.to_string(index=False))
print('SENDER_SOURCE',rd('sender_disease_fixed_candidates.tsv').query("model == 'Source'").to_string(index=False))
(O/'provenance/supplement_complete.txt').write_text('complete\n')
