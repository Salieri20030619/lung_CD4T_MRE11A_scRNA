from pathlib import Path
import sys,json,itertools,hashlib
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from scipy.io import mmread
from scipy.stats import permutation_test,false_discovery_control
O=Path(sys.argv[1]);P=O.parents[1];M=P/'02_analysis/mechanism_step1_20260927_152218_s4V5eF'
cands=pd.concat([pd.read_csv(O/'04_candidate_intersection'/f'{s}_primary_intersection.tsv',sep='\t') for s in ['Naive','ActCD4']],ignore_index=True)
genes=sorted(cands.canonical.unique());manifest=[]
def rec(p):
 p=Path(p);manifest.append({'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()});return p
def wr(a,n):a.to_csv(O/'05_supporting_evidence'/n,sep='\t',index=False,encoding='utf-8',na_rep='NA')
def bh(v):
 a=np.asarray(v,float);q=np.full(len(a),np.nan);ok=np.isfinite(a)
 if ok.any():q[ok]=false_discovery_control(a[ok])
 return q
def summarize(values,control,case):
 c=values[values.group==control].CP10k.to_numpy();s=values[values.group==case].CP10k.to_numpy();p=np.nan
 if len(c)>=2 and len(s)>=2 and np.ptp(np.r_[c,s])>0:p=permutation_test((c,s),lambda a,b:abs(b.mean()-a.mean()),vectorized=False,n_resamples=np.inf,alternative='greater').pvalue
 return dict(n_control=len(c),n_case=len(s),control_mean=c.mean() if len(c) else np.nan,case_mean=s.mean() if len(s) else np.nan,difference=s.mean()-c.mean() if len(c) and len(s) else np.nan,P=p)
# Limited independent mouse direction, with no inferential tests.
f=rec(P/'02_analysis/screen2_CD4_target_wz1lfx/GSE183682/CD4_core.h5ad');a=ad.read_h5ad(f)
raw=a.layers['counts'] if 'counts' in a.layers else a.X
raw=sparse.csr_matrix(raw);assert np.allclose(raw.data,np.rint(raw.data));lib=np.asarray(raw.sum(axis=1)).ravel()
sy=a.var['symbol'].astype(str).str.upper().to_numpy() if 'symbol' in a.var else np.array(a.var_names.str.upper());md=a.obs.copy();md['scope']=md.F08.astype(str).map({'0':'Naive','3':'ActCD4','5':'Treg'})
meta={'GSM5567501':('NS',7),'GSM5567502':('Silica',7),'GSM5567503':('NS',56),'GSM5567504':('Silica',56)}
rows=[]
for scope in ['Core','Naive','ActCD4','Treg']:
 take=md.scope.notna().to_numpy() if scope=='Core' else (md.scope==scope).to_numpy()
 for sample in sorted(md['sample'].astype(str).unique()):
  mask=take&(md['sample'].astype(str).values==sample);n=int(mask.sum())
  for gene in genes:
   ix=np.flatnonzero(sy==gene);v=np.asarray(raw[mask][:,ix].sum(axis=1)).ravel() if len(ix) else np.full(n,np.nan)
   rows.append(dict(dataset='GSE183682',scope=scope,sample=sample,group=meta[sample][0],day=meta[sample][1],gene=gene,cells=n,feature_present=bool(len(ix)),CP10k=np.mean(v/lib[mask]*1e4) if n and len(ix) else np.nan,pct_positive=float(np.mean(v>0)*100) if n and len(ix) else np.nan))
r=pd.DataFrame(rows);wr(r,'GSE183682_candidate_sample_expression.tsv')
dirs=[]
for (scope,day,gene),v in r.groupby(['scope','day','gene']):
 c=v[v.group=='NS'].CP10k.iloc[0];s=v[v.group=='Silica'].CP10k.iloc[0];dirs.append(dict(scope=scope,day=day,gene=gene,control=c,case=s,difference=s-c,interpretation='one library per condition; descriptive only'))
wr(pd.DataFrame(dirs),'GSE183682_candidate_directions.tsv');del a,raw
# Reviewed human CD4 labels and pre-recorded flags, no cell-count gate.
H=P/'11_reports/GSE174725_annotation_validation_20261005_v1/results';a=ad.read_h5ad(rec(H/'broad_annotation_review_SERVER_ONLY.h5ad'));md=pd.read_csv(rec(H/'frozen_annotation_cells_SERVER_ONLY.tsv.gz'),sep='\t',index_col=0).loc[a.obs_names]
raw=sparse.csr_matrix(a.layers['counts']);assert np.allclose(raw.data,np.rint(raw.data));lib=md.total_counts.to_numpy();rows=[]
for filt in ['ReviewedAll','ReviewedFlagExcluded']:
 allowed=np.ones(len(md),bool) if filt=='ReviewedAll' else ~md.sensitivity_flag.astype(bool).to_numpy()
 for scope in ['Core','MemCD4','NaiveCD4','Treg']:
  take=(md.old_CD4.ne('NotInCore') if scope=='Core' else md.old_CD4.eq(scope)).to_numpy()&allowed
  for sample in sorted(md['sample'].unique()):
   mask=take&(md['sample'].values==sample);n=int(mask.sum());group=md.loc[md['sample']==sample,'group'].iloc[0]
   for gene in genes:
    ix=np.flatnonzero(a.var_names.str.upper()==gene);v=np.asarray(raw[mask][:,ix].sum(axis=1)).ravel() if len(ix) else np.full(n,np.nan)
    rows.append(dict(dataset='GSE174725',filter=filt,scope=scope,sample=sample,group=group,gene=gene,cells=n,feature_present=bool(len(ix)),CP10k=np.mean(v/lib[mask]*1e4) if n and len(ix) else np.nan,pct_positive=float(np.mean(v>0)*100) if n and len(ix) else np.nan))
h=pd.DataFrame(rows);wr(h,'GSE174725_candidate_donor_expression.tsv')
hr=[]
for (filt,scope,gene),v in h.groupby(['filter','scope','gene']):
 r=summarize(v.dropna(subset=['CP10k']),'Exposed','Silicosis');r.update(filter=filt,scope=scope,gene=gene,family=filt+'_'+('Core' if scope=='Core' else 'Subtypes')+'_four_selected_candidates',interpretation='limited BALF support; exposed controls; no human-mouse subtype equivalence assumed');hr.append(r)
hr=pd.DataFrame(hr);hr['BH']=hr.groupby('family').P.transform(bh);wr(hr,'GSE174725_candidate_group_comparisons.tsv');del a,raw
# Existing patient-localization extraction only; absent extraction never means no expression.
f=rec(P/'11_reports/GSE296682_patient_localization_20261005_v1/results/all_fixed_and_marker_localization.tsv.gz');loc=pd.read_csv(f,sep='\t');subset=loc[loc.gene.isin(genes)].copy();wr(subset,'GSE296682_existing_candidate_localization.tsv')
wr(pd.DataFrame([dict(gene=g,status='existing extraction available' if g in set(loc.gene) else 'not in legacy fixed-gene extraction; expression unknown',role='patient localization only, not disease validation') for g in genes]),'GSE296682_candidate_coverage_status.tsv')
# Legacy restricted communication and regulatory analyses are supplemental, not selection gates.
for rel,name,col in [('02_analysis/mechanism_step3_20260928_103057_lIRL38/tables/candidate_ligand_prior.tsv','legacy_candidate_ligand_prior.tsv','canonical'),('02_analysis/mechanism_step6_20260928_163015_VoRQqA/tables/TF_target_edges.tsv','legacy_candidate_TF_edges.tsv','target')]:
 a=pd.read_csv(rec(P/rel),sep='\t');c=col if col in a else next((c for c in a if 'target' in c.lower()),a.columns[0]);b=a[a[c].astype(str).str.upper().isin(genes)].copy();b['reuse_limit']='legacy restricted candidates/thresholds; not a new no-cell-gate inference';wr(b,name)
# Candidate-autophagy association: within-sample depth AND pseudotime residual correlation.
# Remove each candidate from the GO signature/control pool before evaluating that candidate.
ce=pd.read_csv(rec(M/'comm_input/CD4_cells.tsv'),sep='\t');gn=pd.read_csv(rec(M/'comm_input/CD4_genes.tsv'),sep='\t');raw=mmread(rec(M/'comm_input/CD4_counts.mtx.gz')).T.tocsr();y=raw.multiply((1e4/ce.library.to_numpy())[:,None]).tocsr();y.data=np.log1p(y.data);gidx=pd.Index(gn.canonical)
members=pd.read_csv(O/'02_mito_autophagy/autophagy_measured_members.tsv',sep='\t').gene.to_list();controls=pd.read_csv(O/'02_mito_autophagy/autophagy_control_members.tsv',sep='\t').gene.to_list();pt=pd.read_csv(O/'03A_state_pseudotime/pseudotime_cell_membership.tsv.gz',sep='\t',index_col=0)
assocs=[];summary=[]
for row in cands.itertuples():
 gene=row.canonical;state=row.state;mi=gidx.get_indexer([g for g in members if g!=gene]);ci=gidx.get_indexer([g for g in controls if g!=gene]);score=np.asarray(y[:,mi].mean(axis=1)).ravel()-np.asarray(y[:,ci].mean(axis=1)).ravel();ex=y[:,gidx.get_loc(gene)].toarray().ravel();rr=[]
 for sample in sorted(ce['sample'].unique()):
  mask=(ce.subtype==state)&(ce['sample']==sample);cc=ce.loc[mask,'cell'];Z=np.column_stack([np.ones(mask.sum()),np.log1p(ce.loc[mask,'library']),pt.loc[cc,'K15_R1']]);a=ex[mask];b=score[mask];a=a-Z@np.linalg.lstsq(Z,a,rcond=None)[0];b=b-Z@np.linalg.lstsq(Z,b,rcond=None)[0];rho=float(a@b/np.sqrt((a@a)*(b@b))) if a@a>1e-15 and b@b>1e-15 else np.nan;rr.append(rho);assocs.append(dict(state=state,gene=gene,sample=sample,cells=int(mask.sum()),partial_rho=rho,removed_from_signature=gene in members,removed_from_controls=gene in controls))
 rr=np.array(rr);ok=np.isfinite(rr);v=rr[ok];p=np.nan
 if len(v)>=2:
  signs=np.array(list(itertools.product([-1,1],repeat=len(v))));p=float((np.abs(signs@v/len(v))>=abs(v.mean())-1e-12).mean())
 summary.append(dict(state=state,gene=gene,n_libraries=len(v),mean_partial_rho=v.mean(),positive_libraries=int((v>0).sum()),negative_libraries=int((v<0).sum()),P=p,in_GO_autophagy=gene in members,interpretation='post-selection exploratory association; candidate removed from score; not causal'))
summary=pd.DataFrame(summary);summary['BH']=bh(summary.P);wr(pd.DataFrame(assocs),'candidate_autophagy_sample_associations.tsv');wr(summary,'candidate_autophagy_association_summary.tsv')
wr(pd.DataFrame(manifest).drop_duplicates('path'),'input_manifest_support.tsv')
print('G183_DIRECTIONS');print(pd.DataFrame(dirs).query("scope in ['Naive','ActCD4']").to_string(index=False))
print('HUMAN_CORE');print(hr[(hr['filter']=='ReviewedAll')&(hr.scope=='Core')][['gene','difference','P','BH','n_control','n_case']].to_string(index=False))
print('AUTOPHAGY');print(summary[['state','gene','mean_partial_rho','P','BH','in_GO_autophagy']].to_string(index=False))
print('SUPPORT_COMPLETE')