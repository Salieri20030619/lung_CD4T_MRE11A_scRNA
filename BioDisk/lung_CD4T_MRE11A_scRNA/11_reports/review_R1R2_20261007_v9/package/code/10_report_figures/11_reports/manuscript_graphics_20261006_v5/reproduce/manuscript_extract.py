from pathlib import Path
import pandas as pd,json,hashlib,re,shutil
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');B=P/'11_reports';O=B/'manuscript_integrated_20261006_v3';T=B/'technical_delivery_20261002_ikzf2_stat1_v1/revision_v2/final_release';out={};manifest=[]
def rd(f):
 f=Path(f);manifest.append(dict(path=str(f.relative_to(P)),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()));return pd.read_csv(f,sep='\t',low_memory=False)
def put(key,d):
 # Export only aggregate scientific results, never individual donor/cell rows.
 d=d.drop(columns=[x for x in ['sample','samples','sample_ids','donor','Subject','cell','omitted','most_changed_donor'] if x in d],errors='ignore');out[key]=json.loads(d.to_json(orient='records',double_precision=15))
fixed=rd(B/'human_first_20261004_v1/results/stage1_legacy_504_repaired_BH.tsv');put('human_all_12_programs',fixed)
put('module_definitions',rd(T/'source_tables/M1_module_definition.tsv'));put('human_module_coverage',rd(T/'source_tables/M2_human_module_coverage.tsv'));put('multi_cohort_module_coverage',rd(T/'source_tables/M1_module_coverage.tsv'))
d=rd(T/'source_tables/M1_group_effects.tsv');put('matched_score_group_effects',d[d.metric.str.startswith('score_')])
d=rd(T/'data/S2_TF_associations.tsv');put('all_TF_program_associations',d)
put('IKZF2_death_audit',rd(T/'data/IKZF2_Core_death_module_direction_audit.tsv'))
h=P/'10_tables/CD4_final';rx=r'OXIDATIVE_PHOSPHORYLATION|ELECTRON_TRANSPORT|MITOCHONDRION_ORGANIZATION|MITOPHAGY|AUTOPHAGY_OF_MITOCHONDRION|OXIDATIVE_STRESS|APOPTOTIC|APOPTOSIS|PYROPT|INFLAMMASOME|T_CELL_ACTIVATION|T_CELL_RECEPTOR|JAK_STAT|ATP_SYNTHESIS|PROTON_TRANSMEMBRANE'
g=rd(h/'pathway_v1/GO_KEGG_GSEA.tsv');put('human_GSEA_focused',g[g.ID.str.contains(rx,case=False,regex=True)])
put('human_GSEA_summary',g.groupby(['analysis','database']).agg(tested=('ID','size'),nominal=('pvalue',lambda x:int(x.lt(.05).sum())),BH05=('p.adjust',lambda x:int(x.lt(.05).sum()))).reset_index())
put('human_GSEA_top_display',g.sort_values(['p.adjust','pvalue']).groupby(['analysis','database']).head(12))
gs=rd(h/'pathway_v1/GO_KEGG_GSVA_tests.tsv');put('human_GSVA_focused',gs[gs.pathway_id.str.contains(rx,case=False,regex=True)])
put('human_GSVA_summary',gs.groupby(['analysis','database']).agg(tested=('pathway_id','size'),nominal=('P.Value',lambda x:int(x.lt(.05).sum())),BH05=('adj.P.Val',lambda x:int(x.lt(.05).sum()))).reset_index())
for name in ['mitochondrial_pathway_sensitivity','mitochondrial_pathway_stability','mitochondrial_leading_edge_summary','sensitivity_GSEA_summary','gene_exclusion_summary']:
 put('human_'+name,rd(h/'pathway_sensitivity_v1'/f'{name}.tsv'))
for name in ['module_group_tests','module_within_G227','module_gene_coverage']:
 put('human_UCell_'+name,rd(h/'module_audit_v1'/f'{name}.tsv'))
put('human_UCell_subtype_tests',rd(h/'subtype_module_audit_v1/subtype_module_tests.tsv'))
# MRE11 associations deliberately remain distinctly named historical evidence.
put('MRE11_human_historical',rd(h/'module_audit_v1/module_MRE11_association.tsv'))
m=P/'02_analysis/GSE314943/stage7_20260905_234200_aiQYL0/tables';parts=[];top=[]
for f in m.glob('*_GSEA.tsv'):
 d=rd(f);d['specification']=f.stem.replace('_GSEA','')
 keep=d.description.str.contains(r'oxidative phosphorylation|electron transport|mitochondrion organization|mitophagy|oxidative stress|apopt|pyropt|inflammasome|T cell activation|T cell receptor|JAK.STAT|ATP synthesis|proton transmembrane',case=False,regex=True,na=False)
 parts.append(d[keep]);top.append(d.sort_values(['padj','pval']).head(8))
put('mouse_GSEA_focused',pd.concat(parts));put('mouse_GSEA_top_display',pd.concat(top));put('mouse_GSEA_summary',rd(m/'GSEA_summary.tsv'))
m=P/'02_analysis/GSE314943/stage8_20260906_003153_6mFKmg/tables'
for name in ['score_group_tests','module_definitions','module_coverage']:
 put('mouse_UCell_'+name,rd(m/f'{name}.tsv'))
put('MRE11_mouse_historical',rd(m/'Mre11a_correlations.tsv'))
# Fixed original counts for background-only atlas; no disease inference from cell totals.
put('lineage_summary',rd(T/'source_tables/Human_preliminary_lineage_summary.tsv'))
(O/'provenance/added_evidence_aggregate.json').write_text(json.dumps(out,ensure_ascii=False,indent=2));(O/'provenance/added_source_manifest.json').write_text(json.dumps(manifest,indent=2))
print('AGGREGATE_EXPORTED',{k:len(v) for k,v in out.items()})
print('FOCAL_HUMAN_PROGRAMS',fixed[(fixed['filter']=='All')&fixed.metric.isin(['score_OXPHOS','score_MitoISR','score_MitoROS','score_Mitophagy','score_ApoExec','score_PyroExec','score_Inflam'])&(fixed.p_nominal<.05)].to_string(index=False))
for tag in ['human_GSEA_focused','human_GSVA_focused','mouse_GSEA_focused']:
 d=pd.DataFrame(out[tag]);c='pathway_id' if tag=='human_GSVA_focused' else 'ID' if tag=='human_GSEA_focused' else 'description'
 q=d[d[c].str.contains('OXIDATIVE_PHOSPHORYLATION|APOPTOSIS|oxidative phosphorylation|^Apoptosis$|T_CELL_RECEPTOR|T cell receptor',case=False,regex=True)]
 print(tag,q.drop(columns=['Description','core_enrichment','leadingEdge','name','leading_edge'],errors='ignore').head(60).to_string(index=False))
