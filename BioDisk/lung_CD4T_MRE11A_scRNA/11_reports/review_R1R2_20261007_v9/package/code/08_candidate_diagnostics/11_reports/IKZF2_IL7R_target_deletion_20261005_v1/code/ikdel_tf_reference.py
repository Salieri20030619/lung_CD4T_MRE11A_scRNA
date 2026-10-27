from pathlib import Path
import os,json,hashlib,fcntl
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
import numpy as np,pandas as pd,h5py
from scipy import stats
from pyscenic.cli.utils import load_signatures
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/IKZF2_IL7R_target_deletion_20261005_v1';R=O/'results';N=P/'02_analysis/scenic1_20260929_103513_xRoEnr/Human/All';V=P/'11_reports/human_alignment_20261005_v1/results';G=P/'11_reports/human_gap_diagnostics_20261005_v1/results';B=P/'02_analysis/scenic2_20261001_092655_ZKNWD6'
lock=(O/'provenance/tf_reference.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);assert not (O/'provenance/tf_reference_complete.json').exists()
spec=O/'code/IKDEL_TF_REFERENCE_SPEC.txt';(O/'provenance/IKDEL_TF_REFERENCE_SPEC.txt').write_bytes(spec.read_bytes())
def rd(p):return pd.read_csv(p,sep='\t')
def wr(x,n):pd.DataFrame(x).to_csv(R/n,sep='\t',index=False,na_rep='NA')
def canon(g):
 a={'MRE11A':'MRE11','ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','PARK2':'PRKN'};s=str(g).upper();return a.get(s,s)
def dec(a):return np.array([x.decode() if isinstance(x,bytes) else str(x) for x in a])
roster=rd(R/'roster_17_to_16_audit.tsv');samples=sorted(roster.loc[roster.included16,'sample']);md=rd(N/'cells.tsv');md=md[(md.subtype=='CD4-IL7R')&md['sample'].isin(samples)].copy()
with h5py.File(N/'score.loom','r') as h:
 genes=dec(h['row_attrs/Gene'][:]);cells=dec(h['col_attrs/CellID'][:]);ix=np.flatnonzero(np.isin(cells,md.cell));X=h['matrix'][:,ix].T
e=pd.DataFrame(X,index=cells[ix],columns=genes).loc[md.cell];md=md.set_index('cell').loc[e.index];means=[];det=[]
for s in samples:
 z=e.loc[md['sample']==s].values;means.append(z.mean(axis=0));det.append((z>0).mean(axis=0))
f=pd.DataFrame({'gene':genes,'expression_rank':stats.rankdata(np.mean(means,axis=0))/len(genes),'detection_rank':stats.rankdata(np.mean(det,axis=0))/len(genes)}).set_index('gene');whole=rd(N/'regulon_edges.tsv');whole=whole[whole.tf.str.upper()!=whole.target.str.upper()]
mods=rd(P/'02_analysis/mechanism_step1_20260927_152218_s4V5eF/tables/module_definition.tsv');excluded={canon(g) for g in mods.gene}|set('IL1R2 PENK TTC39C FOXP1 TNFRSF9 LAG3 MRE11'.split());rows=[];members=[];sets={}
for seed in [17,29,43]:
 deg=whole[whole.seed==seed].groupby('target').tf.nunique();ff=f.copy();ff['log_degree']=np.log2(1+ff.index.map(deg).fillna(0));regs=load_signatures(str(N/f'seed{seed}/regulons.csv'))
 for reg in regs:
  tf=canon(reg.transcription_factor);gs=sorted(g for g in reg.gene2weight if canon(g) not in (excluded|{tf}));assert set(gs)<=set(genes)
  if not gs:continue
  sets[seed,tf]=set(gs);v=ff.loc[gs];rows.append(dict(seed=seed,tf=tf,regulon=reg.name,n_targets=len(gs),expression_rank=v.expression_rank.mean(),detection_rank=v.detection_rank.mean(),log_degree=v.log_degree.mean()));members.extend(dict(seed=seed,tf=tf,target=g) for g in gs)
features=pd.DataFrame(rows);assert not features.duplicated(['seed','tf']).any();base=features[features.tf=='IKZF2'].set_index('seed');alltf=set.intersection(*[set(features[features.seed==s].tf) for s in [17,29,43]])-set('IKZF2 STAT1 IL1R2 PENK TTC39C FOXP1 TNFRSF9 LAG3 MRE11'.split());match=[]
for tf in sorted(alltf):
 z=features[features.tf==tf].set_index('seed').loc[[17,29,43]];size=np.log2(z.n_targets/base.n_targets);ex=z.expression_rank-base.expression_rank;de=z.detection_rank-base.detection_rank;dg=z.log_degree-base.log_degree;ok=bool((size.abs()<=1).all()&(ex.abs()<=.15).all()&(de.abs()<=.15).all()&(dg.abs()<=1).all());distance=float(np.mean(size**2+(ex/.15)**2+(de/.15)**2+dg**2));match.append(dict(tf=tf,eligible=ok,distance=distance,max_abs_log2_size_ratio=size.abs().max(),max_abs_expression_diff=ex.abs().max(),max_abs_detection_diff=de.abs().max(),max_abs_log_degree_diff=dg.abs().max()))
match=pd.DataFrame(match).sort_values(['distance','tf']);chosen=match[match.eligible].head(3).tf.tolist();assert chosen;match['chosen']=match.tf.isin(chosen);wr(match,'reference_TF_matching_all_candidates.tsv');wr(features[features.tf.isin(chosen+['IKZF2'])],'reference_TF_selected_features.tsv');wr(pd.DataFrame(members)[lambda d:d.tf.isin(chosen+['IKZF2'])],'reference_TF_frozen_members.tsv.gz');(O/'provenance/reference_TFs_FROZEN_BEFORE_OUTCOME.json').write_text(json.dumps(dict(chosen=chosen,selection='matching features only; no AUC/program association read yet',spec_sha256=hashlib.sha256(spec.read_bytes()).hexdigest()),indent=2));print('FROZEN_REFERENCES',chosen,flush=True)
# Only now read selected existing AUC columns and fit finite comparison readouts.
old=rd(V/'TF_readout_donors.tsv.gz');old=old[(old.readout=='AUC')&(old.tf=='IKZF2')&(old.scope=='CD4-IL7R')&(old['filter']=='All')&(old.variant=='DropOverlap')];meta=old[old.seed==17].set_index('sample').loc[samples][['group','series','mean_depth']];pg=rd(G/'IL7R_fixed_reference_contrast_donors.tsv').set_index('sample').loc[samples];cov=np.column_stack([np.ones(16),meta.group.eq('CWP'),pd.get_dummies(meta.series,drop_first=True,dtype=float),stats.rankdata(meta.mean_depth)]);assert np.linalg.matrix_rank(cov)==5;b=stats.rankdata(pg.excess);b-=cov@np.linalg.lstsq(cov,b,rcond=None)[0];out=[];don=[]
for seed in [17,29,43]:
 path=B/f'objects/Human_All_{seed}_DropOverlap_auc.csv.gz';head=pd.read_csv(path,nrows=0);rr=features[(features.seed==seed)&features.tf.isin(chosen+['IKZF2'])];cols=rr.regulon.tolist();ac=pd.read_csv(path,usecols=[head.columns[0]]+cols,index_col=0).loc[e.index];dc=ac.groupby(md['sample']).mean().loc[samples]
 for _,v in rr.iterrows():
  x=stats.rankdata(dc[v.regulon]);x-=cov@np.linalg.lstsq(cov,x,rcond=None)[0];rho=float(np.corrcoef(x,b)[0,1]);pv=2*stats.t.sf(abs(rho)*np.sqrt(10/(1-rho*rho)),10);ci=np.tanh(np.arctanh(rho)+np.array([-1,1])*stats.norm.ppf(.975)/3);inter=len(sets[seed,v.tf]&sets[seed,'IKZF2']);union=len(sets[seed,v.tf]|sets[seed,'IKZF2']);out.append(dict(seed=seed,tf=v.tf,n=16,df=10,n_targets=v.n_targets,rho=rho,P=pv,CI_low=ci[0],CI_high=ci[1],overlap_with_IKZF2=inter,Jaccard_with_IKZF2=inter/union,kind='Reused_IKZF2' if v.tf=='IKZF2' else 'Matched_reference'))
  for s in samples:don.append(dict(seed=seed,tf=v.tf,sample=s,AUC=dc.loc[s,v.regulon],excess=pg.loc[s,'excess']))
q=pd.DataFrame(out);q['FDR']=np.nan;ix=q.index[q.kind=='Matched_reference'];p=q.loc[ix,'P'].values;order=np.argsort(p);fq=np.empty(len(p));fq[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);q.loc[ix,'FDR']=fq;q['BH_family']=np.where(q.kind=='Matched_reference',f'fixed_reference_TFs_x_3seeds;m={len(ix)}','REUSED_BASELINE_FDR_NOT_COMPUTED');wr(q,'reference_TF_finite_associations.tsv');wr(don,'reference_TF_donor_scores_SERVER_ONLY.tsv');info=dict(status='COMPLETE',chosen=chosen,models=len(q),new_inferential_models=len(ix),no_new_AUC=True);(O/'provenance/tf_reference_complete.json').write_text(json.dumps(info,indent=2));print(q.to_string(index=False));print(json.dumps(info),flush=True)
