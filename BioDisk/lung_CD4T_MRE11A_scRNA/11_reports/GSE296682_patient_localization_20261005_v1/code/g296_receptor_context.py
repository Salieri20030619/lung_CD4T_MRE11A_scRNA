from pathlib import Path
import json
import pandas as pd,numpy as np
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/GSE296682_patient_localization_20261005_v1';R=O/'results'
(O/'provenance/finite_receptor_context_addendum.txt').write_text('Final finite annotation diagnostic prompted by sparse/heterogeneous author CD4-effector markers: inspect already-fixed IL2RB+IL2RG positive cells using already-fixed TCR2/CD4/CD8/NK/MatureMac3 markers and existing metadata count/depth. No new gene, threshold search, cell relabeling, DE or inferential model. Source comparison uses existing CWP rows with original >=20-cell qualification, donor summaries only, not cross-cohort tests.\n')
z=np.load(R/'official_fixed_expression_SERVER_ONLY.npz',allow_pickle=True);m=pd.read_csv(R/'patient_cell_marker_context_SERVER_ONLY.tsv.gz',sep='\t',index_col=0).loc[z['cells']];genes=list(z['genes']);X=z['values']
def g(s):return X[:,genes.index(s)]
m['IL2RB_positive']=g('IL2RB')>0;m['IL2RG_positive']=g('IL2RG')>0;m['RB_RG_joint']=m.IL2RB_positive&m.IL2RG_positive;m['metadata_counts']=np.expm1(m.log1p_total_counts);m['metadata_features']=np.expm1(m.log1p_nFeaturess_RNA_by_counts)
q=m[m.RB_RG_joint].copy()
for s in 'IL2RB IL2RG CD3D CD3E TRAC CD4 CD8A CD8B NKG7 GNLY KLRD1 PRF1 C1QA C1QC APOE MRC1'.split():q[s]=g(s)[m.RB_RG_joint.to_numpy()]
q.to_csv(R/'receptor_joint_cells_SERVER_ONLY.tsv.gz',sep='\t');rows=[];depth=[]
for (sample,scope),u in m.groupby(['Sample','Cell_Subtype'],observed=True):
 v=u[u.RB_RG_joint];rows.append(dict(sample=sample,scope=scope,N=len(u),RB_RG_positive=len(v),with_TCR2=int(v.TCR2.sum()),with_CD4=int(v.CD4_detected.sum()),with_TCR2_CD4=int((v.TCR2&v.CD4_detected).sum()),with_CD8AB=int(v.CD8AB.sum()),with_NKcyto2=int(v.NKcyto2.sum()),with_MatureMac3=int(v.MatureMac3.sum())))
 depth.append(dict(sample=sample,scope=scope,N=len(u),counts_median=u.metadata_counts.median(),counts_q25=u.metadata_counts.quantile(.25),counts_q75=u.metadata_counts.quantile(.75),features_median=u.metadata_features.median(),MT_pct_median=u.pct_counts_mt.median()))
pd.DataFrame(rows).to_csv(R/'receptor_joint_lineage_summary.tsv',sep='\t',index=False);pd.DataFrame(depth).to_csv(R/'author_scope_depth_context.tsv',sep='\t',index=False)
prior=pd.read_csv(R/'prior_main_fixed_expression_REUSED_SERVER_ONLY.tsv.gz',sep='\t');prior=prior[(prior.group=='CWP')&(prior.cells>=20)&prior.key.isin(['MacAM','MonoMac','MacWide','Core_All','CD4-IL7R_All','CD4-GZMA_All','Treg_All'])];s=prior.groupby(['key','gene']).agg(donors=('sample','nunique'),pct_min=('pct','min'),pct_median=('pct','median'),pct_max=('pct','max'),donors_detected=('pct',lambda x:int((x>0).sum())),donors5=('pct',lambda x:int((x>=.05).sum())),donors10=('pct',lambda x:int((x>=.10).sum()))).reset_index();s['interpretation']='EXISTING_CWP_DONOR_SUMMARY_NO_CROSS_COHORT_TEST';s.to_csv(R/'prior_main_fixed_CWP_coverage_summary_REUSED.tsv',sep='\t',index=False)
focus=['CD4 effector T cell','CD4 helper T cell','Regulatory T cell','CD8 effector T cell','NK cell'];print('JOINT_CONTEXT',pd.DataFrame(rows)[pd.DataFrame(rows).scope.isin(focus)].to_string(index=False));print('DEPTH',pd.DataFrame(depth)[pd.DataFrame(depth).scope.isin(focus)].to_string(index=False));print('MAIN_PRIOR',s[s.gene.isin(['IL15','IL15RA','IL2RB','IL2RG','TNF','TNFRSF1A','TNFRSF1B'])].to_string(index=False))
