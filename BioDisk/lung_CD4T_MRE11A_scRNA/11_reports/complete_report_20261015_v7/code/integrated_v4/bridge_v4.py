from pathlib import Path
import sys,json
import numpy as np,pandas as pd
O=Path(sys.argv[1]);rd=lambda n:pd.read_csv(O/n,sep='\t')
pairs=rd('04_chains/all_sender_receiver_associations.tsv.gz');tf=rd('03_receiver/all_TF_program_associations.tsv.gz');tf=tf[tf.seed==17].copy();tf['endpoint']='AUC::'+tf.x;rna=rd('03_receiver/candidate_RNA_program_associations.tsv');rna['endpoint']='RNA::'+rna.x;prog=pd.concat([tf[tf.x.str.split('(').str[0].str.upper().isin(['IKZF2','STAT1'])],rna],ignore_index=True)
D=rd('03_receiver/candidate_disease_DE.tsv');elig=rd('02_sender/ligand_scope_eligibility.tsv');ranks=rd('02_sender/source_specific_activity_ranks.tsv');prior=rd('inputs/ligand_candidate_prior_weights.tsv');cov=rd('03_receiver/TF_regulon_coverage.tsv');loc=rd('05_human/patient_localization_existing.tsv');phen=pd.concat([rd('01_phenotype/01_CD4_death_group_tests.tsv'),rd('01_phenotype/02_mito_autophagy_group_tests.tsv')]);phen=phen[phen['filter']=='All']
rows=[]
for sender in ['AMlike','AllMac']:
 de=rd(f'02_sender/{sender}_DE_all_genes.tsv.gz').set_index('canonical')
 for scope in ['Core','Naive','ActCD4','Treg']:
  for ligand in ['IL15','IL18','TNF']:
   el=elig[(elig.sender==sender)&(elig.scope==scope)&(elig.ligand==ligand)].iloc[0]
   la=ranks[(ranks.sender==sender)&(ranks.scope==scope)&(ranks.test_ligand==ligand)]
   px=pairs[(pairs.sender==sender)&(pairs.scope==scope)&(pairs.x==ligand)]
   for endpoint in sorted(prog.endpoint.unique()):
    pp=px[px.y==endpoint];gene=endpoint.split('::')[1].split('(')[0].upper();dd=D[(D.scope==scope)&(D.gene==gene)].iloc[0];aa=prog[(prog.scope==scope)&(prog.endpoint==endpoint)]
    for module in sorted(aa.y.unique()):
     r=dict(sender=sender,scope=scope,ligand=ligand,receiver_gene=gene,readout=endpoint.split('::')[0],endpoint=endpoint,module=module,expression_supported=el.expression_supported,supported_receptors=el.supported_receptors,source_DE_logFC=de.loc[ligand,'logFC'] if ligand in de.index else np.nan,source_DE_P=de.loc[ligand,'PValue'] if ligand in de.index else np.nan,source_DE_BH=de.loc[ligand,'FDR'] if ligand in de.index else np.nan,receiver_DE_logFC=dd.logFC,receiver_DE_P=dd.P,receiver_DE_BH=dd.BH)
     for direction in ['Up','Down']:
      rr=la[la.direction==direction];r['activity_'+direction+'_rank']=rr.source_specific_rank.iloc[0] if len(rr) else np.nan;r['activity_'+direction+'_AUPR']=rr.aupr_corrected.iloc[0] if len(rr) else np.nan;r['legacy10pct_supported']=rr.legacy10pct_supported.iloc[0] if len(rr) else False
     for model,prefix in [('Pooled','pooled'),('ConditionAdjusted','condition')]:
      b=pp[pp.model==model];c=aa[(aa.model==model)&(aa.y==module)]
      for k in ['rho','P','BH','loo_same_sign','loo_valid']:
       r[prefix+'_sender_'+k]=b[k].iloc[0] if len(b) else np.nan;r[prefix+'_program_'+k]=c[k].iloc[0] if len(c) else np.nan
     ph=phen[(phen.scope==scope)&(phen.module==module)].iloc[0];r['phenotype_group_delta']=ph.difference;r['phenotype_group_P']=ph.P;r['phenotype_group_BH']=ph.BH
     r['both_edges_nominal']=(r['pooled_sender_P']<.05)&(r['pooled_program_P']<.05);r['both_edges_condition_nominal']=(r['condition_sender_P']<.05)&(r['condition_program_P']<.05)
     r['interpretation']='observational/predictive evidence segments; no causal-chain P'
     rows.append(r)
bridge=pd.DataFrame(rows);bridge.to_csv(O/'04_chains/old_ligand_candidate_evidence_segments.tsv.gz',sep='\t',index=False,na_rep='NA')
# Original TF focal cards remain readable and keep the same RNA/AUC readout in both edges.
f=bridge[(bridge.sender=='AMlike')&(bridge.readout=='AUC')&(bridge.module.isin(['PyroExec','ApoExec','MitoROS']))]
f.to_csv(O/'04_chains/old_TF_chain_cards.tsv',sep='\t',index=False,na_rep='NA')
# Role-specific candidate matrix; no trajectory or time gate.
mat=[]
for dd in D.itertuples():
 if dd.gene in ['IKZF2','STAT1']:
  a=tf[(tf.scope==dd.scope)&(tf.x.str.split('(').str[0].str.upper()==dd.gene)&(tf.model=='Pooled')]
 else:a=rna[(rna.scope==dd.scope)&(rna.x==dd.gene)&(rna.model=='Pooled')]
 hits=a[a.P<.05]
 mat.append(dict(gene=dd.gene,scope=dd.scope,role='regulatory_candidate' if dd.gene in ['IKZF2','STAT1'] else 'response_candidate',DE_tested=dd.DE_tested,DE_logFC=dd.logFC,DE_P=dd.P,DE_BH=dd.BH,nominal_programs=';'.join(hits.y),program_readout='DropOverlap_AUC' if dd.gene in ['IKZF2','STAT1'] else 'RNA',trajectory_required=False,time_required=False,interpretation='retain evidence by role; no automatic final-target declaration'))
pd.DataFrame(mat).to_csv(O/'04_chains/role_specific_candidate_matrix.tsv',sep='\t',index=False,na_rep='NA')
# New high-scoring combinations are audit-only, not promoted based on nominal P or novelty assumptions.
primary_pairs=pairs[pairs.model=='Pooled'];primary_prog=prog[(prog.model=='Pooled')&(prog.P<.05)][['scope','endpoint','y','rho','P','BH']].rename(columns={'y':'module','rho':'program_rho','P':'program_P','BH':'program_BH'});new=primary_pairs[primary_pairs.P<.05].merge(primary_prog,left_on=['scope','y'],right_on=['scope','endpoint']);new['interpretation']='audit-only pooled nominal co-associations; source specificity, prior support and condition confounding require review';new.to_csv(O/'04_chains/all_pooled_two_segment_hits_AUDIT_ONLY.tsv.gz',sep='\t',index=False)
summary={'old_pool_cards':len(bridge),'old_AUC_pyro_double_nominal':int(bridge[(bridge.readout=='AUC')&(bridge.module=='PyroExec')].both_edges_nominal.sum()),'old_AUC_any_double_nominal':int(bridge[bridge.readout=='AUC'].both_edges_nominal.sum()),'old_pool_condition_double_nominal':int(bridge.both_edges_condition_nominal.sum()),'new_nominal_combinations_audit_only':len(new),'patient_scope_labels':sorted(loc.scope.unique().tolist())}
(O/'provenance/bridge_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(json.dumps(summary));print('OLD_DOUBLE_NOMINAL');print(bridge[(bridge.sender=='AMlike')&bridge.both_edges_nominal][['ligand','scope','readout','receiver_gene','module','pooled_sender_rho','pooled_sender_P','pooled_program_rho','pooled_program_P','condition_sender_P','condition_program_P']].to_string(index=False))