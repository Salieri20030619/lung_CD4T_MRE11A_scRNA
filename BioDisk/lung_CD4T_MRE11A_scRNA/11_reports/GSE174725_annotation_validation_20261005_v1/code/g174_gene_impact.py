from pathlib import Path
import json,re,collections
import pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/GSE174725_annotation_validation_20261005_v1';R=O/'results';A=P/'11_reports/GSE174725_reaudit_20261005_v1/results';H=P/'11_reports/human_first_20261004_v1/results';G=P/'11_reports/human_gap_diagnostics_20261005_v1/results';V=P/'11_reports/human_alignment_20261005_v1/results'
def rd(p):return pd.read_csv(p,sep='\t')
mp=rd(A/'unambiguous_gene_mapping.tsv.gz');samples=sorted(mp['sample'].unique());alias=rd(H/'stage2_CD4_mapping.tsv');alias=alias[alias.original.notna()&alias.mapped.notna()];inv={}
for _,v in alias.iterrows():inv.setdefault(v.mapped,set()).add(v.original)
L='IL15 TNF CLSTN1 ADAM15 PLAU DSC2 CXCL2 ICAM2 VTN VEGFB GSTP1 RTN4'.split();lr=rd(V/'all_fixed_candidate_LR_pairs.tsv.gz');lr=lr[lr['from'].isin(L)][['from','to']].drop_duplicates();lr.to_csv(R/'fixed_candidate_receptor_pairs.tsv',sep='\t',index=False)
sets=[('Fixed12Ligands',g,'candidate') for g in L]+[('TF_RNA',g,'TF') for g in ['IKZF2','STAT1']]+[('IL15_complex',g,'receptor_component') for g in ['IL15','IL15RA','IL2RB','IL2RG']]
for _,v in lr.iterrows():
 for g in re.split('[;&_+]',v['to']):sets.append(('LR_'+v['from'],g,'receptor'))
d=rd(H/'stage2_targetset_members.tsv.gz');d=d[(d['filter']=='All')&(d.model=='Source')&(d.direction=='Up')]
for (rec,test),q in d.groupby(['receiver','test']):sets.extend([('Disease_'+rec+'_'+test,g,'disease_target') for g in q.gene])
d=rd(H/'stage1_complete_program_members.tsv');d=d[d.program.isin(['TAct','TReg','Inhib','Naive'])]
for program,q in d.groupby('program'):sets.extend([('Module_'+program,g,'functional_program') for g in q.gene])
ctrl=rd(P/'02_analysis/mechanism_step2_20260927_154919_EsN4gk/tables/human_score_controls.tsv.gz')
for program,q in ctrl[ctrl.module.isin(['TAct','TReg','Inhib','Naive'])].groupby('module'):sets.extend([('Module_'+program+'_background',g,'fixed_functional_background') for g in q.gene])
d=rd(G/'independence_program_members.tsv.gz')
for (program,kind),q in d.groupby(['program','type']):sets.extend([('Prior_'+program+'_'+kind,g,'independence_program') for g in q.gene])
d=rd(G/'matched_control_members.tsv.gz')
for group,q in d.groupby('set'):sets.extend([(group,g,'fixed_matched_background') for g in q.control])
spec=pd.DataFrame(sets,columns=['set','gene','kind']).drop_duplicates();rows=[];unmapped=[]
for gsm in samples:
 q=mp[mp['sample']==gsm];available=set(q.pooled_symbol);ix=q.set_index('pooled_symbol');raw=rd(A/(gsm+'_unmapped_gene_rows.tsv'));raw['sample']=gsm;raw['date_corrupted']=raw.gene.str.match(r'^\d{1,2}-(Mar|Sep|Dec|Oct|Jan|Feb|Apr|May|Jun|Jul|Aug|Nov)$');unmapped.append(raw);uppermap={}
 for symbol in available:uppermap.setdefault(symbol.upper(),set()).add(symbol)
 for _,v in spec.iterrows():
  upper=uppermap.get(v.gene.upper(),set());candidates={v.gene} if v.gene in available else upper if len(upper)==1 else inv.get(v.gene,set())&available
  if len(candidates)==1:
   g=next(iter(candidates));match=ix.loc[g];status='EXACT' if g==v.gene else 'CANONICAL_CASE' if g.upper()==v.gene.upper() else 'EXISTING_VERIFIED_ALIAS';symbol=g;ens=match.ensembl;csv=match.csv_symbol
  else:status=('CSV_SYMBOL_WITH_NONUNIQUE_FEATURE_ID' if raw.gene.str.upper().eq(v.gene.upper()).any() else 'ABSENT_OR_UNRESOLVED_ALIAS') if len(candidates)==0 else 'AMBIGUOUS_ALIAS';symbol=ens=csv=''
  rows.append(dict(sample=gsm,set=v['set'],gene=v.gene,kind=v.kind,status=status,pooled_symbol=symbol,ensembl=ens,csv_symbol=csv))
out=pd.DataFrame(rows);out.to_csv(R/'candidate_receptor_program_feature_audit.tsv.gz',sep='\t',index=False);pd.concat(unmapped).to_csv(R/'all_unresolved_csv_gene_rows.tsv',sep='\t',index=False)
covered=out.status.isin(['EXACT','CANONICAL_CASE','EXISTING_VERIFIED_ALIAS']);summary=out.assign(covered=covered).groupby(['set','sample']).agg(requested=('gene','size'),covered=('covered','sum')).reset_index();summary['fraction']=summary.covered/summary.requested;summary.to_csv(R/'feature_coverage_by_program_donor.tsv',sep='\t',index=False)
common=out.assign(covered=covered).groupby(['set','gene','kind']).agg(donors_covered=('covered','sum'),symbols=('pooled_symbol',lambda x:';'.join(sorted(set(x)-{''})))).reset_index();common.to_csv(R/'program_members_and_common_coverage.tsv',sep='\t',index=False)
focus=out[out.kind.isin(['candidate','TF','receptor_component','receptor'])];print('PROGRAM_COVERAGE',summary.groupby('set').agg(requested=('requested','max'),min_covered=('covered','min'),max_covered=('covered','max')).to_string());print('CASE_REPAIRS',sorted(out[out.status=='CANONICAL_CASE'].gene.unique()));(O/'provenance/gene_impact_complete.json').write_text(json.dumps({'status':'COMPLETE','rows':len(out),'sets':spec['set'].nunique(),'focus_unresolved':int((~focus.status.isin(['EXACT','CANONICAL_CASE','EXISTING_VERIFIED_ALIAS'])).sum()),'mm10_note':'Human ENSG feature annotation verified in preceding audit; true alignment reference cannot be confirmed without author reference/build records. No special correction for selected genes.'},indent=2))
