import json,pandas as pd
F=json.load(open('manuscript_delivery/added_evidence_aggregate.json',encoding='utf8'))
E=json.load(open('revision_delivery/report_evidence_aggregate.json',encoding='utf8'))
def show(k,q=None,cols=None,src=F):
 d=pd.DataFrame(src[k]);d=d.query(q) if q else d
 print('\n'+k+'\n'+(d[cols] if cols else d).to_string(index=False))
show('human_GSEA_summary');show('human_GSVA_summary')
show('human_UCell_module_group_tests',cols=['analysis','model','module','estimate','p_value','fdr_bh'])
show('human_UCell_subtype_tests',"p_value<0.05")
show('human_UCell_module_gene_coverage',cols=['module','expected_n','present_n'])
show('human_mitochondrial_pathway_stability',cols=['analysis','pathway_short','NES_Full','NES_NoRP+MT','p.adjust_NoRP+MT','all_fdr_lt_005'])
show('mouse_UCell_score_group_tests')
show('all_TF_program_associations',"species=='Human' and filter=='All' and variant=='DropOverlap' and scope=='Core' and model=='SourceDep'",['tf','seed','feature','rho','p','FDR'])
show('mouse_GSEA_focused',"specification.str.contains('Full_Standard') and pathway in ['GO:0006119','GO:0097194','GO:0042110','mmu04210','mmu00190']",['specification','pathway','description','NES','pval','padj'])
for k in ['network_display_edges','bootstrap_effect_summary','human_NN_main','mouse_NN','scope_coverage','stage1_models','human_TF_disease','old_TF_group','S4_direction_replication']:
 d=pd.DataFrame(E[k]);print('\nSCHEMA '+k+'\n'+d.head(1).to_string(index=False))
