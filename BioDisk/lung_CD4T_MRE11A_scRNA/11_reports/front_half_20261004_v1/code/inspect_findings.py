from pathlib import Path
import pandas as pd
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/front_half_20261004_v1');R=O/'results'
def rd(n):return pd.read_csv(R/n,sep='\t',low_memory=False)
m=rd('new_link_models.tsv');g=rd('group_context_models.tsv');f=rd('fixed_final_models.tsv')
print('ALL_SPEC_SUMMARY');print(m.groupby(['FDR_family','sample_policy']).agg(rows=('model_id','size'),nominal=('nominal_P_lt_05','sum'),FDR05=('FDR_lt_05','sum'),min_n=('n','min')).to_string())
print('GROUP_BEST');print(g.sort_values('P').head(20)[['kind','scope','y','covariates','n','effect','P','FDR','status']].to_string(index=False))
print('SOURCE_DESIGN');print(rd('source_group_design.tsv').to_string(index=False))
print('ORIGINAL_HUMAN');print(rd('reused_human_original_TF_RNA.tsv').to_string(index=False))
keys=[('TF','receiver','IKZF2_AUC_17','State_Effector','CD4-IL7R'),('TF','receiver','STAT1_AUC_17','State_Effector','Core'),('ligand','MacAM','L_TNF','State_Effector','Core'),('sender_state','MonoMac','SPP1Mac','State_Activation','Core'),('ligand','MacWide','L_TNF','State_Activation','CD4-IL7R'),('ligand','MonoMac','L_CXCL16','State_Inhibitory','Core')]
for kind,sender,x,y,scope in keys:
 print('TRACE',kind,sender,x,y,scope)
 print(m[(m.kind==kind)&(m.sender==sender)&(m.x==x)&(m.y==y)&(m.scope==scope)][['model_id','covariates','sample_policy','n','rho','P','FDR','LOO_min','LOO_max','LOO_all_subjects_valid','excluded_ids']].to_string(index=False))
print('NICHENET_FROZEN');nn=rd('reused_step3_nichenet_ranks.tsv');print(nn.to_string(index=False))
print('MOUSE_PROGRAM_OVERLAP');print(rd('mouse_target_program_overlap.tsv').to_string(index=False))
