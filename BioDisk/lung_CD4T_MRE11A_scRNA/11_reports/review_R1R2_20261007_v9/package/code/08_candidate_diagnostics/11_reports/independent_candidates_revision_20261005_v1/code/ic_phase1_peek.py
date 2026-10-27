from pathlib import Path
import os,json,hashlib
os.environ['OPENBLAS_NUM_THREADS']='1'
import pandas as pd,numpy as np
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/independent_candidates_revision_20261005_v1';D=O/'results/source_evidence';R=O/'results';pool='IL15 TNF IKZF2 STAT1 CLSTN1 ADAM15 PLAU DSC2 CXCL2 ICAM2 VTN VEGFB GSTP1 RTN4 IL18'.split()
def rd(tag):
 f=next(D.glob(tag+'.tsv*'));return pd.read_csv(f,sep='\t',low_memory=False)
def sh(tag,df,cols=None,n=35):
 print('\nFOCUS',tag,'rows',len(df),'COLS',list(df.columns),flush=True);print((df[[x for x in cols if x in df]] if cols else df).head(n).to_string(index=False),flush=True)
manifest=[];de=[]
for f in (P/'11_reports/human_first_20261004_v1/results').glob('*DE*tsv*'):
 try:q=pd.read_csv(f,sep='\t')
 except Exception:continue
 if 'gene' in q:
  x=q[q.gene.str.upper().isin(pool)].copy();x['source_file']=str(f);de.append(x);manifest.append(dict(path=str(f),rows=len(q),selected=len(x),sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
if de:
 q=pd.concat(de,ignore_index=True);q.to_csv(R/'human_existing_DE_fixed15.tsv.gz',sep='\t',index=False);sh('TF_HUMAN_DE',q[q.gene.isin(['IKZF2','STAT1'])],['source_file','receiver','scope','filter','model','gene','logFC','PValue','FDR'],80)
(O/'provenance/additional_DE_source_manifest.json').write_text(json.dumps(manifest,indent=2))
q=rd('human_TF_disease');sh('human_TF_disease',q[(q['filter']=='All')&(q.readout=='RNA')],['scope','tf','readout','model','n','n_control','n_case','effect','P','FDR','status'],100)
q=rd('human_TF_roster');sh('human_TF_roster',q,['scope','tf','roster','model','n','n_control','n_case','effect','P','FDR','status'],90)
for tag in ['human_TF_program_bridge','human_ligand_program','human_ligand_evidence','human_ligand_activity','human_sender_DE','mouse_sender','mouse_candidate_program','IL18_TNF_direct_phenotype','IL18_TNF_proximal','STAT1_IFN_context','old_baseline','network_edges','network_targets','IL7R_IKZF2_excess','IL7R_IKZF2_influence','GZMA_TReg_independence']:
 q=rd(tag);sh(tag,q,n=2)
q=rd('old_TF_group');q=q[(q.species=='Human')&(q['filter']=='All')];sh('OLD_HUMAN_AUC_GROUP_NOMINAL',q[q.p<.05],['scope','tf','variant','seed','model','n_control','n_case','coef','p','FDR','df'],45)
q=rd('old_TF_program');q=q[(q.species=='Human')&(q['filter']=='All')];sh('OLD_STAT1_PROGRAM_NOMINAL',q[(q.tf=='STAT1')&(q.p<.05)],['scope','tf','variant','seed','feature','model','n','rho','p','FDR','df'],40)
print('PEEK_COMPLETE',flush=True)
