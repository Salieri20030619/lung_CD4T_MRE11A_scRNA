from pathlib import Path
import sys,pandas as pd,numpy as np,json
O=Path(sys.argv[1]);prefix=(O/'code/mouse_stats_v4.py').read_text(encoding='utf-8').split("phen=rd('01_phenotype/sample_program_scores.tsv')")[0];ns={};exec(compile(prefix,str(O/'code/mouse_stats_v4.py'),'exec'),ns);assoc=ns['assoc'];samples=ns['SAMPLES'];mods=ns['MODS'];cand=ns['CAND']
rec=pd.read_csv(O/'03_receiver/receptor_candidate_sample_expression.tsv.gz',sep='\t');weights=pd.read_csv(O/'06_support/Core_fixed_composition_weights.tsv',sep='\t').set_index('state').fixed_weight;Y=pd.read_csv(O/'06_support/Core_fixed_composition_program_inputs.tsv',sep='\t',index_col=0).loc[samples];X=None
for state,w in weights.items():
 a=rec[(rec.scope==state)&rec.gene.isin(cand)].pivot(index='sample',columns='gene',values='CP10k').reindex(index=samples,columns=cand);X=a*w if X is None else X+a*w
z=assoc(X,Y,'Core_fixedcomposition','none',0,'RNA_program');z.to_csv(O/'06_support/Core_fixed_composition_RNA_program.tsv',sep='\t',index=False,na_rep='NA')
act=pd.read_csv(O/'03_receiver/TF_all_regulon_sample_activity.tsv.gz',sep='\t');send=pd.read_csv(O/'02_sender/all_ligand_sample_expression.tsv.gz',sep='\t');elig=pd.read_csv(O/'02_sender/ligand_scope_eligibility.tsv',sep='\t');rows=[]
for scope in ['Core','Naive','ActCD4','Treg']:
 for seed in [29,43]:
  a=act[(act.scope==scope)&(act.seed==seed)];A=a.pivot(index='sample',columns='regulon',values='mean_auc').loc[samples];A=A[[c for c in A if c.rsplit('(',1)[0].upper() in ['IKZF2','STAT1']]];A.columns=['AUC::'+c for c in A]
  for source in ['AMlike','AllMac']:
   lig=elig[(elig.sender==source)&(elig.scope==scope)&elig.expression_supported].ligand.to_list();X=send[(send.scope==source)&send.gene.isin(lig)].pivot(index='sample',columns='gene',values='CP10k').reindex(index=samples,columns=lig);rows.append(assoc(X,A,scope,source,seed,'sender_TF_seed_sensitivity'))
z=pd.concat(rows,ignore_index=True);z.to_csv(O/'06_support/sender_TF_seed_sensitivity.tsv.gz',sep='\t',index=False,na_rep='NA');print('SENSITIVITY_COMPLETE',len(z));print(z[(z.x=='IL18')&(z.y=='AUC::Ikzf2(+)')&(z.scope=='ActCD4')&(z.sender=='AMlike')][['seed','model','rho','P','BH','loo_same_sign']].to_string(index=False))