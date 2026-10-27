from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,shutil
O=Path(__file__).resolve().parents[1];W=O.parent
inputs={
 'receiver':W/'stepwise_diagnostic/inputs/receiver_all_positive_cell_samples.tsv',
 'sender':W/'stepwise_diagnostic/inputs/sender_IL18_MacAM_all_samples.tsv',
 'human_scores':W/'revision_v2/technical_delivery/source_tables/M2_human_sample_scores.tsv',
 'mouse_scores':W/'revision_v2/technical_delivery/source_tables/M1_sample_scores.tsv',
 'legacy_activity':W/'revision_v2/technical_delivery/data/S2_TF_sample_activity.tsv',
 'legacy_associations':W/'evidence_story_20261002/existing_S2_TF_to_six_candidates.tsv'}
tables={};manifest=[]
for k,p in inputs.items():
 tables[k]=pd.read_csv(p,sep='\t');dest=O/'inputs'/p.name;shutil.copy2(p,dest)
 manifest.append(dict(name=k,source=str(p),copy=str(dest.relative_to(O)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
six=['IL1R2','PENK','TTC39C','FOXP1','TNFRSF9','LAG3'];joined=[];audits=[]
for (sp,arm,scope),z in tables['receiver'].groupby(['species','filter','scope']):
 q=tables['human_scores' if sp=='Human' else 'mouse_scores']
 if sp=='Mouse':q=q[q.dataset=='GSE314943']
 qscope=('All' if scope=='Core' else scope) if sp=='Mouse' else ('Core' if scope=='Core' else 'Subtype:'+scope)
 q=q[(q['filter']==arm)&(q.scope==qscope)].copy()
 cols=['sample','cells','library']+[v+'_'+g for g in six for v in ['expr','pct']]
 a=z.merge(q[cols],on='sample',how='left',validate='many_to_one',indicator=True)
 assert a['_merge'].eq('both').all() and a.cells.eq(a.receiver_cells).all()
 assert np.allclose(a.library,a.receiver_library)
 a=a.drop(columns=['_merge','cells','library'])
 if sp=='Human':
  s=tables['sender'];s=s[s.sender_filter=='All']
  a=a.merge(s.drop(columns=['group','series']),on='sample',how='left',validate='many_to_one')
  for row in a[(a.seed==17)&(a.tf=='IKZF2')].itertuples():
   assert s[s['sample']==row.sample].group.iloc[0]==row.group
 else:
  for c in ['sender_cells','sender_mean_CP10k','sender_library','sender_epi']:a[c]=np.nan
 joined.append(a)
d=pd.concat(joined,ignore_index=True);d.to_csv(O/'inputs/aligned_receiver_candidates_sender.tsv',sep='\t',index=False)
rosters=[];transitions=[];summary=[]
for (sp,arm,scope),z in d[(d.seed==17)&(d.tf=='IKZF2')].groupby(['species','filter','scope']):
 specs=[('Receiver',t) for t in [1,5,10,20]]+([('Paired',t) for t in [1,5,10,20]] if sp=='Human' else [])
 for branch,t in specs:
  keep=(z.receiver_cells>=t)&(z.sender_cells>=t if branch=='Paired' else True)
  q=z[keep];key=f'{sp}/{arm}/{scope}/{branch}{t}'
  summary.append(dict(roster_id=key,species=sp,filter=arm,scope=scope,branch=branch,threshold=t,n=len(q),n_control=int(q.group.isin(['Control','Saline']).sum()),n_case=int(q.group.isin(['CWP','Silica']).sum()),sample_ids=';'.join(sorted(q['sample']))))
  for _,r in z.iterrows():
   reasons=[]
   if r.receiver_cells<t:reasons.append(f'receiver_cells={r.receiver_cells}<{t}')
   if branch=='Paired' and (pd.isna(r.sender_cells) or r.sender_cells<t):reasons.append(f'sender_cells={r.sender_cells}<{t}' if pd.notna(r.sender_cells) else 'sender_missing')
   rosters.append(dict(roster_id=key,species=sp,filter=arm,scope=scope,branch=branch,threshold=t,sample=r['sample'],group=r.group,series=r.series,receiver_cells=r.receiver_cells,sender_cells=r.sender_cells,included=not reasons,exclusion_reason=';'.join(reasons)))
 for branch in ['Receiver']+(['Paired'] if sp=='Human' else []):
  previous=set(z['sample'])
  for t in [1,5,10,20]:
   rr=[x for x in rosters if x['roster_id']==f'{sp}/{arm}/{scope}/{branch}{t}'];now={x['sample'] for x in rr if x['included']}
   for x in rr:
    if x['sample'] in previous-now:transitions.append(dict(**x,transition=f'previous_to_{t}'))
   previous=now
pd.DataFrame(rosters).to_csv(O/'results/roster_membership_and_reasons.tsv',sep='\t',index=False)
pd.DataFrame(summary).to_csv(O/'results/roster_summary.tsv',sep='\t',index=False)
pd.DataFrame(transitions).to_csv(O/'results/sequential_exclusions.tsv',sep='\t',index=False)
pd.DataFrame(manifest).to_csv(O/'provenance/input_manifest.tsv',sep='\t',index=False)
focus=pd.DataFrame(summary);focus=focus[(focus.species=='Human')&(focus['filter']=='All')&(focus.scope=='CD4-IL7R')]
print(focus.to_string(index=False))
print('ALL_INPUT_CELL_COUNTS_AND_DEPTHS_MATCH; original edgeR unchanged; no scientific models run in this stage')
