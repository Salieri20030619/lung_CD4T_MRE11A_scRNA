from pathlib import Path
import os,json
os.environ['OPENBLAS_NUM_THREADS']='1'
import pandas as pd,numpy as np
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/independent_candidates_revision_20261005_v1';D=O/'results/source_evidence';R=O/'results'
pool='IL15 TNF IKZF2 STAT1 CLSTN1 ADAM15 PLAU DSC2 CXCL2 ICAM2 VTN VEGFB GSTP1 RTN4 IL18'.split()
def rd(tag):return pd.read_csv(next(D.glob(tag+'.tsv*')),sep='\t',low_memory=False)
def show(tag,q,cols=None,n=80):
 print('\nTABLE',tag,'n',len(q),flush=True);print((q[[c for c in cols if c in q]] if cols else q).head(n).to_string(index=False),flush=True)
for tag,col in [('human_ligand_evidence','ligand'),('human_ligand_program','ligand'),('human_ligand_activity','test_ligand'),('human_sender_DE','gene')]:
 q=rd(tag);q=q[q[col].str.upper().isin(pool)];q.to_csv(R/(tag+'_fixed15.tsv.gz'),sep='\t',index=False)
 if tag=='human_sender_DE':show(tag,q,['gene','sender','model','logFC','PValue','FDR'])
 if tag=='human_ligand_program':
  q=q[(q['filter']=='All')&(q.sender!='TNK')];show('ALL_MAC_ASSOC_NOMINAL',q[q.P<.05],['ligand','sender','receiver','model','test','direction','n','rho','P','FDR'],120)
 if tag=='human_ligand_evidence':
  z=q[(q['filter']=='All')&(q.sender!='TNK')];a=z.groupby('ligand').agg(contexts=('status','size'),LR=('LR_expression_supported','sum'),activity=('aupr_corrected',lambda s:s.notna().sum()),positive_activity=('aupr_corrected',lambda s:(s>0).sum()),assoc=('association_P',lambda s:s.notna().sum()),nominal=('association_P',lambda s:(s<.05).sum()),FDR05=('association_FDR',lambda s:(s<.05).sum()));a.to_csv(R/'candidate_context_counts.tsv',sep='\t');show(tag,a)
for tag in ['mouse_ligand_activity','mouse_ligand_context','independent_G296','old_M1','human_program_disease']:
 q=rd(tag);show(tag,q.head(3))
q=rd('IL18_TNF_direct_phenotype');show('DIRECT_EDGES',q[['edge','x','y','covariates']].drop_duplicates(),n=45)
q=q[(q.P<.05)&q.edge.isin(['ligand_phenotype','ligand_program'])];show('IL18_TNF_DIRECT_NOMINAL',q,['axis','sender','scope','edge','x','y','covariates','n','rho','P','FDR'],60)
q=rd('IL7R_IKZF2_excess');show('EXCESS_ALL',q,['variant','seed','endpoint','n','rho','P'],24)
q=rd('old_TF_program');q=q[(q.species=='Human')&(q['filter']=='All')&(q.tf=='STAT1')&(q.variant=='DropOverlap')];show('STAT1_DROPOVERLAP_NOMINAL',q[q.p<.05],['scope','seed','feature','model','n','rho','p','FDR'],50)
q=rd('human_program_disease');print('PROGRAM_VALUES', {c:q[c].dropna().astype(str).unique()[:30].tolist() for c in q if any(x in c for x in ['scope','model','program','feature','endpoint','kind'])})
print('DONE',flush=True)
