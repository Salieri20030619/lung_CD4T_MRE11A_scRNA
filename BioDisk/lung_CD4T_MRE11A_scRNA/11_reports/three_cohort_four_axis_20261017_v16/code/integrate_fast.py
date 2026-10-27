from pathlib import Path
import pandas as pd,numpy as np,json,os,signal
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/three_cohort_four_axis_20261017_v16')
for p in Path('/proc').glob('[0-9]*/cmdline'):
 try:
  if p.read_bytes().endswith(str(O/'code/intersect.py').encode()+b'\x00'):os.kill(int(p.parent.name),signal.SIGTERM)
 except (FileNotFoundError,PermissionError,ProcessLookupError):pass
rd=lambda n:pd.read_csv(O/n,sep='\t',low_memory=False)
d=rd('tables/DE_all_models.tsv');d=d[d.model.isin(['314_disease','joint_stratified','274_time'])].copy();d['direction']=np.where(d.logFC>0,'up','down');a=d[d.PValue<.05].copy();a['dataset']=np.where(a.model=='274_time','GSE274372','GSE314943')
f=rd('tables/function_membership_extended.tsv');fs=f.groupby('gene').agg(function_mechanisms=('mechanism',lambda z:';'.join(sorted(set(z)))),function_sources=('source',lambda z:';'.join(sorted(set(z)))),function_only_literature=('source',lambda z:(z=='primary_literature_extension').all())).reset_index();a=a.merge(fs,how='left',on='gene');a['function_axis']=a.function_sources.notna();a['function_only_literature']=a.function_only_literature.fillna(False).astype(bool)
lr=rd('tables/expressed_LR_with_decoy_literature.tsv');lr=lr[lr.sender=='AllMac'];dr=lr.assign(gene=lr.receptor_parts.str.split(';')).explode('gene');agg=lambda z:';'.join(sorted(set(z)))
dg=dr.groupby(['dataset','scope','gene']).agg(direct_ligands=('ligand_parts',agg)).reset_index();a=a.merge(dg,on=['dataset','scope','gene'],how='left');a['direct_receptor']=a.direct_ligands.notna()
nn=rd('tables/NicheNet_top250_targets.tsv');ng=nn.groupby(['gene','family','direction']).agg(predicted_ligands=('ligand',agg),max_aupr_corrected=('aupr_corrected','max'),best_target_rank=('target_rank','min')).reset_index();a=a.merge(ng,on=['gene','family','direction'],how='left');om=nn[nn.model=='274_time'].groupby('gene').agg(predicted_ligands=('ligand',agg),max_aupr_corrected=('aupr_corrected','max'),best_target_rank=('target_rank','min'))
for i in a.index[a.contrast=='time_omnibus']:
 g=a.at[i,'gene']
 if g in om.index:
  for col in om.columns:a.at[i,col]=om.at[g,col]
a['predicted_downstream']=a.predicted_ligands.notna();a['mac_axis']=a.direct_receptor|a.predicted_downstream;a['three_axis']=a.function_axis&a.mac_axis;a['ligands']=[ ';'.join(sorted(set((str(x)+';'+str(y)).split(';'))-{'nan',''})) for x,y in zip(a.direct_ligands,a.predicted_ligands)]
for model,prefix in [('314_disease','314'),('joint_stratified','joint')]:
 z=d[d.model==model][['gene','scope','logFC','PValue','FDR']].rename(columns={'logFC':prefix+'_logFC','PValue':prefix+'_P','FDR':prefix+'_FDR'});a=a.merge(z,on=['gene','scope'],how='left')
cps=[]
for sc in ['Core','Naive','ActCD4','Treg']:
 c=rd(f'tables/joint_stratified__{sc}__logCPM.tsv');c['scope']=sc;c['183_D7_descriptive_delta']=c.GSM5567502-c.GSM5567501;c['183_D56_descriptive_delta']=c.GSM5567504-c.GSM5567503;cps.append(c[['gene','scope','183_D7_descriptive_delta','183_D56_descriptive_delta']])
a=a.merge(pd.concat(cps),on=['gene','scope'],how='left');a.to_csv(O/'tables/nominal_DE_axis_membership.tsv',sep='\t',index=False);s=a[a.three_axis].copy();s.to_csv(O/'tables/three_axis_candidates_all_contrasts.tsv',sep='\t',index=False)
# Matched full ligand-receptor-target routes; keep omnibus evidence explicitly separate.
keys=s[['gene','family','scope','dataset','model','contrast','direction']];direct=keys.merge(dr,on=['gene','scope','dataset']);direct['route']='direct_receptor';direct['target_rank']=np.nan;direct['activity_rank']=np.nan;direct['ligand']=direct.ligand_parts;direct['receptor']=direct.receptor_parts
nt=keys[keys.contrast!='time_omnibus'].merge(nn,on=['gene','family','scope','dataset','model','direction']);ot=keys[keys.contrast=='time_omnibus'].merge(nn.drop(columns=['family','direction']),on=['gene','scope','dataset','model']);nt=pd.concat([nt,ot],ignore_index=True);net=nt.merge(lr.rename(columns={'ligand':'original_ligand','ligand_parts':'ligand'}),on=['dataset','scope','ligand']);net['route']='NicheNet_predicted_target';net['receptor']=net.receptor_parts
cols=['gene','family','scope','route','ligand','receptor','interaction','group','evidence','annotation','target_rank','activity_rank'];pd.concat([direct[cols],net[cols]]).drop_duplicates().to_csv(O/'tables/candidate_complete_routes.tsv',sep='\t',index=False)
counts={'DE':int(a.gene.nunique()),'DE_function':int(a.loc[a.function_axis,'gene'].nunique()),'DE_function_mac':int(s.gene.nunique()),'disease':int(s.loc[s.model!='274_time','gene'].nunique()),'temporal':int(s.loc[s.model=='274_time','gene'].nunique())};(O/'tables/intersection_counts.json').write_text(json.dumps(counts,indent=2),encoding='utf-8');print('COUNTS',counts);print('CANDIDATES',';'.join(sorted(s.gene.unique())));print(s.sort_values(['direct_receptor','joint_FDR'],ascending=[False,True]).drop_duplicates(['gene','scope'])[['gene','scope','model','logFC','PValue','FDR','direct_receptor','predicted_downstream','183_D7_descriptive_delta','183_D56_descriptive_delta']].to_string(index=False));(O/'logs/intersection_complete.txt').write_text('success',encoding='utf-8')