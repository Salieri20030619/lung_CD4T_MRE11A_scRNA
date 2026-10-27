from pathlib import Path
import pandas as pd,numpy as np,json
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/three_cohort_four_axis_20261017_v16')
def rd(n):return pd.read_csv(O/n,sep='\t',low_memory=False)
d=rd('tables/DE_all_models.tsv');f=rd('tables/function_membership_extended.tsv');lr=rd('tables/expressed_LR_with_decoy_literature.tsv');nn=rd('tables/NicheNet_top250_targets.tsv');d=d[d.model.isin(['314_disease','joint_stratified','274_time'])].copy();d['direction']=np.where(d.logFC>0,'up','down');nom=d[d.PValue<.05];rows=[];edges=[]
cp_cache={sc:rd(f'tables/joint_stratified__{sc}__logCPM.tsv').set_index('gene') for sc in ['Core','Naive','ActCD4','Treg']}
for r in nom.itertuples():
 ds='GSE274372' if r.model=='274_time' else 'GSE314943';direct=lr[(lr.dataset==ds)&(lr.scope==r.scope)&(lr.sender=='AllMac')&lr.receptor_parts.str.split(';').map(lambda z:r.gene in z)]
 n=nn[(nn.gene==r.gene)&(nn.scope==r.scope)&(nn.family==r.family)&(nn.direction==r.direction)] if r.contrast!='time_omnibus' else nn[(nn.gene==r.gene)&(nn.model==r.model)&(nn.scope==r.scope)]
 ff=f[f.gene==r.gene]
 a=r._asdict();a.pop('Index');a['function_axis']=len(ff)>0;a['function_mechanisms']=';'.join(sorted(ff.mechanism.unique()));a['function_sources']=';'.join(sorted(ff.source.unique()));a['function_only_literature']=len(ff)>0 and (ff.source=='primary_literature_extension').all();a['direct_receptor']=len(direct)>0;a['predicted_downstream']=len(n)>0;a['mac_axis']=len(direct)>0 or len(n)>0;a['three_axis']=a['function_axis'] and a['mac_axis'];a['ligands']=';'.join(sorted(set(direct.ligand_parts)|set(n.ligand)));a['max_aupr_corrected']=n.aupr_corrected.max() if len(n) else np.nan;a['best_target_rank']=n.target_rank.min() if len(n) else np.nan
 if ds=='GSE314943':
  cp=cp_cache[r.scope]
  for day,c,t in [(7,'GSM5567501','GSM5567502'),(56,'GSM5567503','GSM5567504')]:a[f'183_D{day}_descriptive_delta']=cp.loc[r.gene,t]-cp.loc[r.gene,c] if r.gene in cp.index else np.nan
  orig=d[(d.model=='314_disease')&(d.scope==r.scope)&(d.gene==r.gene)];a['314_logFC']=orig.logFC.iloc[0] if len(orig) else np.nan;a['314_P']=orig.PValue.iloc[0] if len(orig) else np.nan
  joint=d[(d.model=='joint_stratified')&(d.scope==r.scope)&(d.gene==r.gene)];a['joint_logFC']=joint.logFC.iloc[0] if len(joint) else np.nan;a['joint_P']=joint.PValue.iloc[0] if len(joint) else np.nan;a['joint_FDR']=joint.FDR.iloc[0] if len(joint) else np.nan
 rows.append(a)
 if a['three_axis']:
  for e in direct.itertuples():edges.append(dict(gene=r.gene,family=r.family,scope=r.scope,route='direct_receptor',ligand=e.ligand_parts,receptor=e.receptor_parts,interaction=e.interaction,group=e.group,evidence=e.evidence,annotation=e.annotation,target_rank=np.nan,activity_rank=np.nan))
  for e in n.itertuples():
   matching=lr[(lr.dataset==ds)&(lr.scope==r.scope)&(lr.sender=='AllMac')&(lr.ligand_parts==e.ligand)]
   for le in matching.itertuples():edges.append(dict(gene=r.gene,family=r.family,scope=r.scope,route='NicheNet_predicted_target',ligand=e.ligand,receptor=le.receptor_parts,interaction=le.interaction,group=le.group,evidence=le.evidence,annotation=le.annotation,target_rank=e.target_rank,activity_rank=e.activity_rank))
a=pd.DataFrame(rows);a.to_csv(O/'tables/nominal_DE_axis_membership.tsv',sep='\t',index=False);s=a[a.three_axis].copy();s.to_csv(O/'tables/three_axis_candidates_all_contrasts.tsv',sep='\t',index=False);pd.DataFrame(edges).drop_duplicates().to_csv(O/'tables/candidate_complete_routes.tsv',sep='\t',index=False)
print('UNION_COUNTS',json.dumps({'DE':a.gene.nunique(),'DE_function':a.loc[a.function_axis,'gene'].nunique(),'DE_function_mac':s.gene.nunique(),'disease':s.loc[s.model!='274_time','gene'].nunique(),'temporal':s.loc[s.model=='274_time','gene'].nunique()}))
print('CANDIDATES',';'.join(sorted(s.gene.unique())))
print(s.sort_values(['direct_receptor','joint_FDR'],ascending=[False,True]).drop_duplicates(['gene','scope'])[['gene','scope','model','logFC','PValue','FDR','direct_receptor','predicted_downstream','183_D7_descriptive_delta','183_D56_descriptive_delta']].to_string(index=False))