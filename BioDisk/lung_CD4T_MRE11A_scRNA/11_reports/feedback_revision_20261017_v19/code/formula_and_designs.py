from pathlib import Path
import pandas as pd,numpy as np,json,math,itertools
O=Path(r'C:\Users\31937\lung_feedback_revision_20261017_v19');V=Path(r'C:\Users\31937\lung_append_analysis_20261017_v18');B=Path(r'C:\Users\31937\lung_review_20261017_v14\revision_v14');T=O/'tables'
s=pd.read_csv(V/'tables/sample_gene_and_program_scores.tsv',sep='\t');e=pd.read_csv(V/'tables/sample_selected_gene_expression.tsv.gz',sep='\t');c=pd.read_csv(V/'tables/module_coverage.tsv',sep='\t');keys=['dataset','variant','scope','sample'];x=e.merge(s[keys+['library']],on=keys,how='left',validate='many_to_one');assert x.library.notna().all();x['reconstructed_CP10k']=x.UMI/x.library*10000;assert np.allclose(x.CP10k,x.reconstructed_CP10k,equal_nan=True);x.to_csv(T/'R03_gene_counts_with_library.tsv.gz',sep='\t',index=False)
checks=[];reconstructed={}
for _,row in s[s.variant.isin(['Raw','BG1'])].iterrows():
 key=tuple(row[k] for k in keys);z=x[(x.dataset==row.dataset)&(x.variant==row.variant)&(x.scope==row.scope)&x['sample'].eq(row['sample'])].set_index('gene');rr={}
 for col in s.columns:
  if col.startswith('G_'):g=col[2:];rr[col]=float(np.log1p(z.at[g,'reconstructed_CP10k'])) if g in z.index else np.nan
  if col.startswith('M_'):
   cov=c[(c.dataset==row.dataset)&c.module.eq(col[2:])].iloc[0];genes=cov.mapped.split(';');rr[col]=float(np.log1p(z.loc[genes,'reconstructed_CP10k']).mean()) if cov.status=='SCORABLE' else np.nan
 for col,val in rr.items():
  err=abs(val-row[col]) if np.isfinite(val) and np.isfinite(row[col]) else 0 if pd.isna(val) and pd.isna(row[col]) else np.inf;assert err<1e-12;checks.append(dict(dataset=row.dataset,variant=row.variant,scope=row.scope,sample=row['sample'],variable=col,error=err,status='PASS'))
 reconstructed[key]=rr
weights=[]
for _,row in s[s.variant.isin(['FixedComp','BG1FixedComp'])].iterrows():
 base='Raw' if row.variant=='FixedComp' else 'BG1';z=s[(s.dataset==row.dataset)&s.variant.eq(base)&s.scope.isin(['Naive','ActCD4','Treg'])&s['sample'].eq(row['sample'])];assert len(z)==3
 for state in ['Naive','ActCD4','Treg']:
  zz=z[z.scope==state].iloc[0];weights.append(dict(scheme='v18_equal_state_standardization',dataset=row.dataset,variant=row.variant,sample=row['sample'],state=state,weight=1/3,included_cells=int(zz.cells),source_variant=base))
 for col in [cc for cc in s if cc.startswith(('G_','M_'))]:
  vals=[reconstructed[(row.dataset,base,state,row['sample'])][col] for state in ['Naive','ActCD4','Treg']];val=float(np.log1p(np.mean(np.expm1(vals)))) if col.startswith('G_') else float(np.mean(vals));err=abs(val-row[col]);assert err<1e-12;checks.append(dict(dataset=row.dataset,variant=row.variant,scope=row.scope,sample=row['sample'],variable=col,error=err,status='PASS'))
pd.DataFrame(checks).to_csv(T/'R03_score_reconstruction_QA.tsv',sep='\t',index=False);pd.DataFrame(weights).to_csv(T/'R03_equal_state_weights_actual.tsv',sep='\t',index=False)
old=pd.read_csv(B/'supplement/S18_fixed_composition_weights.tsv',sep='\t');old['scheme']='v14_mean_sample_proportion_not_equal';d=s[(s.dataset=='GSE314943')&s.variant.eq('Raw')&s.scope.isin(['Naive','ActCD4','Treg'])].copy();d['fraction']=d.cells/d.groupby('sample').cells.transform('sum');calc=d.groupby('scope').fraction.mean();old['reconstructed_from_six_samples']=old.state.map(calc);assert np.allclose(old.fixed_weight,old.reconstructed_from_six_samples);old.to_csv(T/'R03_old_and_new_weight_distinction.tsv',sep='\t',index=False)
# Explicit design matrices and exchangeability blocks, using the frozen analysis encoding.
a=pd.read_csv(V/'tables/all_associations.tsv',sep='\t');designs=[];blocks=[];rosters=[]
for (co,var,scope,model),z0 in a.groupby(['cohort','variant','scope','model']):
 ds=['GSE314943'] if co=='314' else ['GSE314943','GSE183682'] if co=='Joint' else ['GSE274372'];d=s[s.dataset.isin(ds)&s.variant.eq(var)&s.scope.eq(scope)].sort_values(['dataset','sample']).reset_index(drop=True)
 if co=='274_min10':d=d[d.cells>=10].reset_index(drop=True)
 design_id='|'.join([str(co),var,scope,model]);C={'Intercept':np.ones(len(d))};lab=(d.dataset.astype(str)+'_'+d.day.astype(str)) if model.startswith('Joint') else d.day.astype(str) if model=='TimeDayAdjusted' else pd.Series(['all']*len(d))
 if model.startswith('Joint') or model=='TimeDayAdjusted':
  dum=pd.get_dummies(lab,drop_first=True)
  for col in dum:C['stratum_'+str(col)]=dum[col].to_numpy(float)
 if 'Condition' in model:C['Silica']=d.group.eq('Silica').to_numpy(float)
 block=lab+'_'+d.group if 'Condition' in model else lab
 counts=block.value_counts().to_dict();nperm=math.prod(math.factorial(v) for v in counts.values());assert len(d)==int(z0.n.iloc[0]);assert nperm==int(z0.permutations_total.iloc[0]);assert ';'.join(d['sample'])==z0.sample_list.iloc[0]
 mat=np.column_stack(list(C.values()));rank=np.linalg.matrix_rank(mat);assert rank==mat.shape[1]
 for i,rr in d.iterrows():
  rosters.append(dict(design_id=design_id,cohort=co,variant=var,variant_label={'FixedComp':'等权状态标准化','BG1FixedComp':'BG1加等权状态标准化'}.get(var,var),scope=scope,model=model,dataset=rr.dataset,day=rr.day,sample=rr['sample'],group=rr.group,cells=rr.cells,library=rr.library,design_rank=rank,residual_df=len(d)-rank-1))
  blocks.append(dict(design_id=design_id,sample=rr['sample'],permuted_variable='rank(y), then residualize with same C',block=block.iloc[i],block_n=counts[block.iloc[i]],permutations_total=nperm,evaluated=min(nperm,10000),method='full enumeration' if nperm<=10000 else 'MonteCarlo10000_seed20261017'))
  for col,val in C.items():designs.append(dict(design_id=design_id,sample=rr['sample'],column=col,value=float(val[i])))
pd.DataFrame(designs).to_csv(T/'R04_design_matrices.tsv',sep='\t',index=False);pd.DataFrame(blocks).to_csv(T/'R04_permutation_blocks.tsv',sep='\t',index=False);pd.DataFrame(rosters).to_csv(T/'R04_actual_sample_rosters.tsv',sep='\t',index=False)
# All actual group-effect input rosters; score table also provides the complete effect-model variables.
eff=pd.read_csv(V/'tables/programme_group_effects.tsv',sep='\t');eff[['cohort','variant','scope','sample_list','method','n','P','family']].to_csv(T/'R04_group_effect_model_index.tsv',sep='\t',index=False)
summary={'reconstructed_score_values':len(checks),'max_error':float(pd.DataFrame(checks).error.max()),'designs':len(set(r['design_id'] for r in rosters)),'all_designs_match_saved_n_permutations_and_sample_order':True,'old_weights_reproduced':True};(O/'reconstruction_QA.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(json.dumps(summary))