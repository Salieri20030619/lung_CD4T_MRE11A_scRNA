from pathlib import Path
import os,json,sys
os.environ['OPENBLAS_NUM_THREADS']='1'
import pandas as pd,numpy as np
from scipy.io import mmread
from scipy import sparse
from threadpoolctl import threadpool_limits
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');A=P/'02_analysis/mechanism_step2_20260927_154919_EsN4gk';O=P/'11_reports/human_gap_diagnostics_20261005_v1';R=O/'results';H=P/'11_reports/human_first_20261004_v1/results'
sys.path.insert(0,str(A));from shared import canonical
rd=lambda p:pd.read_csv(p,sep='\t',low_memory=False)
def wr(d,n):pd.DataFrame(d).to_csv(R/n,sep='\t',index=False,na_rep='NA')
m=rd(A/'input/human_cells.tsv');g=rd(A/'input/human_genes.tsv').symbol.astype(str);canon=np.array([canonical(v) for v in g]);ctrl=rd(A/'tables/human_score_controls.tsv.gz');cg=ctrl[ctrl.module=='TReg'].gene.tolist();members=rd(R/'independence_program_members.tsv.gz')
with threadpool_limits(limits=2):x=mmread(A/'input/human_counts.mtx').tocsr()
need=sorted(set(cg)|set(members[members.program.str.startswith('TReg')].gene));assert set(need)<=set(canon)
logs={}
for gene in need:logs[gene]=np.log1p(np.asarray(x[canon==gene].sum(axis=0)).ravel()/m.library*1e4)
control=np.mean([logs[s] for s in cg],axis=0);old=rd(A/'objects/human_cell_scores.tsv.gz');m=m.merge(old[['cell','bg','mt','nfeature','score_Naive','score_TAct','score_TReg']],on='cell',validate='one_to_one');out=[]
for name,q in members[members.program.str.startswith('TReg')&(members.type=='target')].groupby('program'):
 m[name]=np.mean([logs[s] for s in q.gene],axis=0)-control
 if name=='TReg_original':assert np.max(abs(m[name]-m.score_TReg))<1e-10
 for sample,v in m[m.subtype=='CD4-GZMA'].groupby('sample'):out.append(dict(sample=sample,group=v.group.iloc[0],series=v.series.iloc[0],cells=len(v),library=v.library.sum(),mean_depth=v.library.mean(),eligible20=len(v)>=20,scope='CD4-GZMA',filter='All',lineage='CD4',program=name,members_n=len(q),score=v[name].mean()))
ind=rd(R/'independence_program_donors.tsv.gz');draft=ind[ind.scope=='CD4-GZMA'];wr(draft,'INVALID_mapped_space_TReg_reconstruction.tsv.gz');ind=pd.concat([ind[ind.scope!='CD4-GZMA'],pd.DataFrame(out)],ignore_index=True);wr(ind,'independence_program_donors.tsv.gz')
# Per-cell scores for the frozen disease gene sets in the exact mapped NicheNet space.
mp=rd(H/'stage2_CD4_mapping.tsv');assert len(mp)==len(g) and np.array_equal(mp.original.astype(str),g)
need=sorted(set(members[members.program.str.startswith('IL7R')&(members.type=='target')].gene));ilogs={}
for gene in need:ilogs[gene]=np.log1p(np.asarray(x[mp.mapped.eq(gene).to_numpy()].sum(axis=0)).ravel()/m.library*1e4)
for name,q in members[members.program.str.startswith('IL7R')&(members.type=='target')].groupby('program'):m[name]=np.mean([ilogs[s] for s in q.gene],axis=0)
for gene in ['IKZF2','FOXP3','CTLA4','TIGIT','IL2RA','ENTPD1','LRRC32','IKZF4']:
 m['count_'+gene]=np.asarray(x[canon==gene].sum(axis=0)).ravel()
wr(m,'independent_program_cell_scores_SERVER_ONLY.tsv.gz')
check=[]
for name in [c for c in m if c.startswith('IL7R_nominal')]:
 a=m[m.subtype=='CD4-IL7R'].groupby('sample')[name].mean();b=ind[ind.program==name].set_index('sample').score;idx=a.index.intersection(b.index);err=float(np.max(abs(a.loc[idx]-b.loc[idx])));assert err<1e-10;check.append(dict(program=name,max_error=err))
mapped=set(mp.mapped);audit=pd.DataFrame({'control_gene':cg,'present_in_nichenet_mapped': [s in mapped for s in cg]});wr(audit,'TReg_control_mapping_mismatch.tsv')
receipt={'original_TReg_cell_reconstruction_error':float(np.max(abs(m.TReg_original-m.score_TReg))),'frozen_controls':len(cg),'controls_absent_from_mapped_space':audit.loc[~audit.present_in_nichenet_mapped,'control_gene'].tolist(),'IL7R_donor_reconstruction':check,'resolution':'Reconstruct TReg in original M2 canonical raw space; IL7R disease targets remain NicheNet mapped space; no change to old published scores.'};(O/'provenance/module_reconstruction.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt,indent=2))
