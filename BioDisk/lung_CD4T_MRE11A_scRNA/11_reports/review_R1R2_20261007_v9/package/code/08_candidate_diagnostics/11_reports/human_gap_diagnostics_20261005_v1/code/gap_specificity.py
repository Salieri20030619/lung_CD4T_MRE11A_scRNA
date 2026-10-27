from pathlib import Path
import os,ast,json
os.environ['OPENBLAS_NUM_THREADS']='1'
import pandas as pd,numpy as np
from scipy import stats
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/human_gap_diagnostics_20261005_v1';R=O/'results';V=P/'11_reports/human_alignment_20261005_v1/results';A=P/'02_analysis'
tree=ast.parse((O/'code/gap_compute.py').read_text());exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef)],type_ignores=[]),'<frozen functions>','exec'),globals())
assert (O/'code/SPECIFICITY_FOLLOWUP_SPEC.txt').exists();(O/'provenance/SPECIFICITY_FOLLOWUP_SPEC.txt').write_bytes((O/'code/SPECIFICITY_FOLLOWUP_SPEC.txt').read_bytes())
ind=rd(R/'independence_program_donors.tsv.gz');ct=rd(R/'matched_control_donor_scores.tsv.gz');tf=rd(V/'TF_readout_donors.tsv.gz');tf=tf[(tf.scope=='CD4-IL7R')&(tf['filter']=='All')&(tf.tf=='IKZF2')&(tf.readout=='AUC')];endpoints=pd.concat([ct,ind[ind.program.isin(['IL7R_nominal_original','IL7R_nominal_nonRiboMito'])].rename(columns={'program':'endpoint'})],ignore_index=True);rows=[]
for (variant,seed),d in tf.groupby(['variant','seed']):
 base=d[['sample','group','series','cells','mean_depth','value']].rename(columns={'cells':'receiver_cells','mean_depth':'receiver_depth','value':'x'})
 for endpoint,q in endpoints.groupby('endpoint'):
  dd=complete(base.merge(q[['sample','score','members_n']],on='sample').rename(columns={'score':'y'}));rows.append(dict(variant=variant,seed=seed,endpoint=endpoint,members_n=q.members_n.iloc[0],**cor(dd,'SourceReceiverDepth')))
rows=pd.DataFrame(rows);test=bh(rows[rows.endpoint!='MatchedTargetSubset'].reset_index(drop=True),['variant','seed']);subset=rows[rows.endpoint=='MatchedTargetSubset'].copy();subset['FDR']=np.nan;subset['BH_family']='DESCRIPTIVE_MATCHING_RESTRICTED_TARGET_SUBSET';wr(pd.concat([test,subset]),'matched_control_specificity_models.tsv')
# Existing cells, existing AUC, donor-balanced descriptions: no inferential cell p-values.
auc=rd(R/'existing_TF_auc_selected_cell_context_SERVER_ONLY.tsv.gz');sc=rd(R/'independent_program_cell_scores_SERVER_ONLY.tsv.gz');features=['TReg_original','TReg_minus_IKZF2','TReg_minus_IKZF2_network','IL7R_nominal_original','IL7R_nominal_nonRiboMito','IL7R_nominal_minus_IKZF2_network','nfeature'];auc=auc.merge(sc[['cell']+features],on='cell',validate='many_to_one');auc=auc[(auc.subtype=='CD4-IL7R')&(auc.tf=='IKZF2')];within=[];bins=[]
ref=auc.drop_duplicates('cell').copy();sz=ref.groupby('sample').size();ref=ref[ref['sample'].map(sz)>=20];ref['cell_weight']=1/ref['sample'].map(sz);cuts={}
for axis in ['score_Naive','library']:
 rr=ref.sort_values(axis);w=rr.cell_weight.cumsum()/rr.cell_weight.sum();cuts[axis]=np.interp([.25,.5,.75],w,rr[axis]);assert len(set(cuts[axis]))==3
(O/'provenance/common_state_cutoffs.json').write_text(json.dumps({k:v.tolist() for k,v in cuts.items()},indent=2))
for (variant,seed,sample),d in auc.groupby(['variant','seed','sample']):
 tag=dict(variant=variant,seed=seed,sample=sample,group=d.group.iloc[0],series=d.series.iloc[0],cells=len(d),mean_depth=d.library.mean())
 if len(d)<20:continue
 X=np.column_stack([np.ones(len(d)),stats.rankdata(d.library),stats.rankdata(d.nfeature)])
 for endpoint in ['IL7R_nominal_original','IL7R_nominal_nonRiboMito']:
  x=stats.rankdata(d.auc);y=stats.rankdata(d[endpoint]);x-=X@np.linalg.lstsq(X,x,rcond=None)[0];y-=X@np.linalg.lstsq(X,y,rcond=None)[0];rho=np.corrcoef(x,y)[0,1] if min(np.std(x),np.std(y))>1e-10 else np.nan;within.append(dict(tag,endpoint=endpoint,rho=rho,interpretation='WITHIN_DONOR_DESCRIPTIVE_NO_CELL_P'))
 for axis in ['score_Naive','library']:
  cut=np.digitize(d[axis],cuts[axis])
  for quart,g in d.groupby(cut):bins.append(dict(tag,axis=axis,quartile=int(quart)+1,bin_cells=len(g),auc=g.auc.mean(),original=g.IL7R_nominal_original.mean(),nonRiboMito=g.IL7R_nominal_nonRiboMito.mean(),axis_mean=g[axis].mean(),raw_auc=d.auc.mean(),raw_original=d.IL7R_nominal_original.mean(),raw_nonRiboMito=d.IL7R_nominal_nonRiboMito.mean()))
within=pd.DataFrame(within);bins=pd.DataFrame(bins);wr(within,'within_donor_IKZF2_program_correlations.tsv');wr(bins,'within_donor_state_quartile_scores.tsv.gz')
keys=['variant','seed','sample','group','series','cells','mean_depth','axis'];standard=bins.groupby(keys)[['auc','original','nonRiboMito','raw_auc','raw_original','raw_nonRiboMito']].mean().reset_index();support=bins.groupby(keys).agg(observed_bins=('quartile','nunique'),min_bin_cells=('bin_cells','min')).reset_index();standard=standard.merge(support,on=keys);wr(standard,'common_state_support_donors.tsv');standard=standard[standard.observed_bins==4];std=[]
for (variant,seed,axis),d in standard.groupby(['variant','seed','axis']):
 for endpoint in ['original','nonRiboMito']:
  for mode in ['MeasuredCommonRoster','EqualBinWeights']:
   xc='raw_auc' if mode=='MeasuredCommonRoster' else 'auc';yc='raw_'+endpoint if mode=='MeasuredCommonRoster' else endpoint;q=d.rename(columns={xc:'x',yc:'y','cells':'receiver_cells','mean_depth':'receiver_depth'});std.append(dict(variant=variant,seed=seed,axis=axis,endpoint=endpoint,mode=mode,**cor(q,'SourceReceiverDepth')))
std=bh(pd.DataFrame(std),['variant','seed','axis','mode']);wr(std,'state_quartile_standardized_donor_models.tsv')
# GZMA rare-cell composition: observed marker counts, no annotation changes.
gz=sc[sc.subtype=='CD4-GZMA'];comp=[]
for sample,d in gz.groupby('sample'):
 pos=d.count_IKZF2>0;z=d[pos];v=dict(sample=sample,group=d.group.iloc[0],series=d.series.iloc[0],cells=len(d),IKZF2_positive=int(pos.sum()),IKZF2_pct=pos.mean(),IKZF2_umi=d.count_IKZF2.sum(),mean_depth=d.library.mean(),positive_mean_depth=z.library.mean(),negative_mean_depth=d[~pos].library.mean(),positive_TReg_noTF=z.TReg_minus_IKZF2.mean(),negative_TReg_noTF=d[~pos].TReg_minus_IKZF2.mean())
 for g in ['FOXP3','IL2RA','CTLA4','TIGIT']:v[g+'_codetected']=int((z['count_'+g]>0).sum())
 comp.append(v)
comp=pd.DataFrame(comp);wr(comp,'GZMA_IKZF2_positive_cell_composition.tsv');group=comp.groupby(['series','group']).agg(donors=('sample','size'),cells=('cells','sum'),IKZF2_positive=('IKZF2_positive','sum'),IKZF2_umi=('IKZF2_umi','sum'),donors_positive=('IKZF2_positive',lambda x:sum(x>0)),median_pct=('IKZF2_pct','median'),median_depth=('mean_depth','median')).reset_index();wr(group,'GZMA_IKZF2_source_composition_summary.tsv')
summary={'matched_control_models':len(rows),'within_donor_rows':len(within),'state_standardized_models':len(std),'GZMA_cells':len(gz),'GZMA_IKZF2_positive':int((gz.count_IKZF2>0).sum())};(O/'provenance/specificity_complete.json').write_text(json.dumps(summary,indent=2));print('CONTROL_TESTS',test[['variant','seed','endpoint','members_n','n','df','rho','P','FDR']].to_string(index=False));print('WITHIN_DONOR',within.groupby(['variant','seed','endpoint']).agg(donors=('sample','size'),negative=('rho',lambda x:sum(x<0)),median=('rho','median'),min=('rho','min'),max=('rho','max')).to_string());print('GZMA_COMPOSITION',group.to_string(index=False));print('STANDARDIZED',std[['variant','seed','axis','endpoint','rho','P','FDR']].to_string(index=False))
