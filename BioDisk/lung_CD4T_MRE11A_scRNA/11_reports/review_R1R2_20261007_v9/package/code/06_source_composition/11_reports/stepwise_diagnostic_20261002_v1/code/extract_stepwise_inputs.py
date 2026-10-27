from pathlib import Path
import pandas as pd,numpy as np,json,hashlib

ROOT=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
OUT=ROOT/'11_reports/stepwise_diagnostic_20261002_v1'
for d in ['inputs','code','results','figures','provenance']:(OUT/d).mkdir(parents=True,exist_ok=True)
S1=ROOT/'02_analysis/scenic1_20260929_103513_xRoEnr'
S2=ROOT/'02_analysis/scenic2_20261001_092655_ZKNWD6'
S3=ROOT/'02_analysis/scenic3_20261001_123452_Kvp9h5'
TECH=ROOT/'11_reports/technical_delivery_20261002_ikzf2_stat1_v1/revision_v2/final_release'
manifest=[]
def rd(p,**kw):
 manifest.append({'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size})
 return pd.read_csv(p,**kw)
human=rd(TECH/'source_tables/M2_human_sample_scores.tsv',sep='\t')
mouse=rd(TECH/'source_tables/M1_sample_scores.tsv',sep='\t');mouse=mouse[mouse.dataset=='GSE314943']
rna=rd(TECH/'data/RNA_sample_summary.tsv',sep='\t')
sender=rd(TECH/'prior_review_repaired/scenic3_review/results/reviewed_sender_expression.tsv',sep='\t')
sender=sender[(sender.sender=='MacAM')&(sender.ligand=='IL18')]
sender.to_csv(OUT/'inputs/sender_IL18_MacAM_all_samples.tsv',sep='\t',index=False)
identity=rd(S3/'tables/mouse_identity_sample.tsv',sep='\t')
allrows=[]
for sp,scopes in [('Human',['Core','CD4-IL7R']),('Mouse',['Core','ActCD4'])]:
 for arm in ['All','BG1']:
  cells=rd(S1/sp/arm/'cells.tsv',sep='\t')
  assert cells.cell.is_unique
  for seed in [17,29,43]:
   cols=['Cell']+(['IKZF2(+)','STAT1(+)'] if sp=='Human' else ['Ikzf2(+)','Stat1(+)'])
   auc=rd(S2/'objects'/f'{sp}_{arm}_{seed}_DropOverlap_auc.csv.gz',usecols=cols).rename(columns={'Cell':'cell',cols[1]:'IKZF2',cols[2]:'STAT1'})
   m=cells.merge(auc,on='cell',validate='one_to_one');assert len(m)==len(cells)
   for scope in scopes:
    z=m if scope=='Core' else m[m.subtype==scope]
    a=z.groupby('sample',observed=True).agg(group=('group','first'),series=('series','first'),receiver_cells=('cell','size'),receiver_library=('library','mean'),receiver_bg=('bg','mean'),IKZF2=('IKZF2','mean'),STAT1=('STAT1','mean')).reset_index()
    qscope=('All' if scope=='Core' else scope) if sp=='Mouse' else ('Core' if scope=='Core' else 'Subtype:'+scope)
    q=(mouse if sp=='Mouse' else human)
    q=q[(q['filter']==arm)&(q.scope==qscope)]
    qcols=['sample','cells','library']+[c for c in q if c.startswith('score_') or c.startswith('detect_')]
    a=a.merge(q[qcols],on='sample',validate='one_to_one')
    assert (a.receiver_cells==a.cells).all(),(sp,arm,scope,'cell mismatch')
    assert np.allclose(a.receiver_library,a.library),(sp,arm,scope,'depth mismatch')
    a=a.drop(columns=['cells','library'])
    if sp=='Mouse' and scope=='ActCD4':
     ident=identity[(identity['filter']==arm)&(identity.seed==seed)]
     a=a.merge(ident[['sample','TregNoTF','FOXP3_det','FOXP3_IL2RA_det']],on='sample',validate='one_to_one')
    for tf in ['IKZF2','STAT1']:
     rr=rna[(rna.species==sp)&(rna['filter']==arm)&(rna.scope==scope)&(rna.tf==tf)][['sample','cells','mean_CP10k','pct','present']]
     v=a.merge(rr,on='sample',validate='one_to_one')
     assert (v.cells==v.receiver_cells).all()
     v['mean_auc']=v[tf];v['tf']=tf;v['species']=sp;v['filter']=arm;v['scope']=scope;v['seed']=seed
     v=v.drop(columns=['IKZF2','STAT1','cells']).rename(columns={'mean_CP10k':'tf_mean_CP10k','pct':'tf_pct','present':'tf_present'})
     allrows.append(v)
   print('FROZEN_AUC_AGGREGATED',sp,arm,seed,flush=True)
data=pd.concat(allrows,ignore_index=True)
data.to_csv(OUT/'inputs/receiver_all_positive_cell_samples.tsv',sep='\t',index=False)
keys=['species','filter','scope','seed','tf','sample']
assert not data.duplicated(keys).any()
old=rd(TECH/'data/S2_TF_sample_activity.tsv',sep='\t');old=old[old.variant=='DropOverlap']
cmp=data.merge(old,on=keys,suffixes=('','_old'),validate='one_to_one')
assert np.allclose(cmp.mean_auc,cmp.mean_auc_old,atol=1e-12)
manifest=pd.DataFrame(manifest).drop_duplicates('path')
manifest.to_csv(OUT/'provenance/input_manifest.tsv',sep='\t',index=False)
summary={'rows':len(data),'legacy_eligible_auc_means_reproduced':len(cmp),'max_absolute_auc_difference':float(np.max(np.abs(cmp.mean_auc-cmp.mean_auc_old))),'scenic_recomputed':False,'minimum_receiver_cells_included':int(data.receiver_cells.min()),'output':str(OUT)}
(OUT/'provenance/input_verification.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
