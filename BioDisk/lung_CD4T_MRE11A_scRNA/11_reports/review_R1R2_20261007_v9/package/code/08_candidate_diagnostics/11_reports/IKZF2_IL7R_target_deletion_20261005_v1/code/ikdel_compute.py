from pathlib import Path
import os,json,hashlib,math,fcntl,inspect,sys
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
import numpy as np,pandas as pd,h5py
from scipy import stats
from scipy.optimize import linear_sum_assignment
from pyscenic.aucell import create_rankings,aucell4r
from ctxcore.genesig import GeneSignature
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/IKZF2_IL7R_target_deletion_20261005_v1';R=O/'results';N=P/'02_analysis/scenic1_20260929_103513_xRoEnr/Human/All';V=P/'11_reports/human_alignment_20261005_v1/results';G=P/'11_reports/human_gap_diagnostics_20261005_v1/results';T=P/'11_reports/technical_delivery_20261002_ikzf2_stat1_v1/revision_v2/final_release'
lock=(O/'provenance/compute.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
assert not (O/'provenance/compute_complete.json').exists(),'Already complete; do not duplicate'
spec=O/'code/IKDEL_FROZEN_SPEC.txt';assert spec.exists();(O/'provenance/IKDEL_FROZEN_SPEC.txt').write_bytes(spec.read_bytes())
def rd(p):return pd.read_csv(p,sep='\t',low_memory=False)
def wr(x,n):pd.DataFrame(x).to_csv(R/n,sep='\t',index=False,na_rep='NA')
def dec(a):return np.array([x.decode() if isinstance(x,bytes) else str(x) for x in a])
roster=rd(R/'roster_17_to_16_audit.tsv');samples=sorted(roster.loc[roster.included16,'sample']);assert len(samples)==16
md=rd(N/'cells.tsv');md=md[(md.subtype=='CD4-IL7R')&md['sample'].isin(samples)].copy();assert len(md)==1361
pg=rd(G/'IL7R_fixed_reference_contrast_donors.tsv').set_index('sample').loc[samples]
baseold=rd(V/'TF_readout_donors.tsv.gz');baseold=baseold[(baseold.readout=='AUC')&(baseold.scope=='CD4-IL7R')&(baseold.tf=='IKZF2')&(baseold['filter']=='All')&(baseold.variant=='DropOverlap')]
meta=baseold[baseold.seed==17].set_index('sample').loc[samples][['group','series','cells','mean_depth']].copy();meta['sample']=samples;meta['y']=pg.excess;assert meta.cells.sum()==1361 and meta.y.notna().all()
members=rd(P/'11_reports/human_first_20261004_v1/results/stage2_targetset_members.tsv.gz');mm=members[(members.receiver=='CD4-IL7R')&(members['filter']=='All')&(members.model=='Source')&(members.test=='Nominal05')&(members.direction=='Up')];assert 'IKZF2' not in set(mm.gene)
with h5py.File(N/'score.loom','r') as h:
 genes=dec(h['row_attrs/Gene'][:]);cells=dec(h['col_attrs/CellID'][:]);ix=np.flatnonzero(np.isin(cells,md.cell));X=h['matrix'][:,ix].T
expression=pd.DataFrame(X,index=cells[ix],columns=genes).loc[md.cell];assert expression.shape==(1361,8938);assert np.isfinite(expression.values).all()
md=md.set_index('cell').loc[expression.index];means=[];detect=[]
for s in samples:
 z=expression.loc[md['sample']==s].values;means.append(z.mean(axis=0));detect.append((z>0).mean(axis=0))
features=pd.DataFrame({'gene':genes,'mean_log':np.mean(means,axis=0),'detection':np.mean(detect,axis=0)}).set_index('gene');features['expression_rank']=stats.rankdata(features.mean_log)/len(features);features['detection_rank']=stats.rankdata(features.detection)/len(features)
net=rd(T/'data/TF_network_all_edges.tsv');net=net[(net.species=='Human')&(net['filter']=='All')&(net.tf=='IKZF2')&net.retained_DropOverlap].copy();assert not net.target.eq('IKZF2').any();whole=rd(N/'regulon_edges.tsv');whole=whole[whole.tf.str.upper()!=whole.target.str.upper()].copy()
# Technical delivery canonicalized symbols to upper case; restore exact original-universe spelling only.
case_map={g.upper():g for g in genes};assert len(case_map)==len(genes)
net['target_canonical']=net.target;net['target']=net.target.map(lambda x:case_map.get(x.upper(),x));whole['target']=whole.target.map(lambda x:case_map.get(x.upper(),x));wr(net[['seed','target_canonical','target']].drop_duplicates(),'exact_case_symbol_mapping.tsv')
(O/'provenance/case_mapping_addendum.txt').write_text('Before scoring, exact one-to-one case-only mapping restores original score.loom/regulons.csv spellings (C18ORF25 -> C18orf25; C12ORF73 -> C12orf73). No gene added/dropped, no alias guess, and all planned sizes/thresholds unchanged. Original per-cell baseline must reproduce.\n')
def cor(x,d=meta):
 n=len(d);cov=np.column_stack([np.ones(n),d.group.eq('CWP').to_numpy(float),pd.get_dummies(d.series,drop_first=True,dtype=float),stats.rankdata(d.mean_depth)]);rank=np.linalg.matrix_rank(cov);df=n-rank-1
 if rank<cov.shape[1] or df<3:return dict(rho=np.nan,P=np.nan,CI_low=np.nan,CI_high=np.nan,df=df,status='NOT_IDENTIFIABLE')
 a=stats.rankdata(x);b=stats.rankdata(d.y);a-=cov@np.linalg.lstsq(cov,a,rcond=None)[0];b-=cov@np.linalg.lstsq(cov,b,rcond=None)[0]
 if min(a.std(),b.std())<1e-12:return dict(rho=np.nan,P=np.nan,CI_low=np.nan,CI_high=np.nan,df=df,status='CONSTANT')
 rho=float(np.clip(np.corrcoef(a,b)[0,1],-1,1));p=2*stats.t.sf(abs(rho)*np.sqrt(df/max(1e-15,1-rho*rho)),df);ci=np.tanh(np.arctanh(np.clip(rho,-1+1e-12,1-1e-12))+np.array([-1,1])*stats.norm.ppf(.975)/np.sqrt(df-1));return dict(rho=rho,P=p,CI_low=ci[0],CI_high=ci[1],df=df,status='EXPLORATORY_APPROX_FISHER_CI')
def bh(p):
 p=np.asarray(p);ix=np.argsort(p);z=p[ix]*len(p)/np.arange(1,len(p)+1);out=np.empty(len(p));out[ix]=np.minimum(1,np.minimum.accumulate(z[::-1])[::-1]);return out
all_models=[];all_members=[];all_pairs=[];all_status=[];donors=[];checks=[];mechanical=[];all_features=[];loo=[]
for seed in [17,29,43]:
 q=net[net.seed==seed].sort_values(['weight','target'],ascending=[False,True]);assert len(q)=={17:68,29:169,43:156}[seed] and q.target.is_unique;target=q.target.tolist();weights=q.set_index('target').weight.to_dict();assert set(target)<=set(genes)
 f=features.copy();deg=whole[whole.seed==seed].groupby('target').tf.nunique();f['indegree']=f.index.map(deg).fillna(0).astype(int);f['log_degree']=np.log2(1+f.indegree);sd=max(f.log_degree.std(),1e-12);f['seed']=seed;all_features.append(f.loc[target].reset_index())
 kmax=math.ceil(.2*len(target));high=target[:kmax];pool=target[kmax:];M=f.loc[high,['expression_rank','detection_rank','log_degree']].values;C=f.loc[pool,['expression_rank','detection_rank','log_degree']].values;diff=M[:,None,:]-C[None,:,:];valid=(abs(diff[:,:,0])<=.35)&(abs(diff[:,:,1])<=.35)&(abs(diff[:,:,2])<=1);cost=np.sum((diff/np.array([1,1,sd]))**2,axis=2)
 specs=[dict(name='Baseline',kind='Baseline',fraction=0.,control_seed=0,deleted=[],retained=target)]
 for frac in [.1,.2]:
  k=math.ceil(frac*len(target));specs.append(dict(name=f'High_{frac}',kind='HighWeight',fraction=frac,control_seed=0,deleted=target[:k],retained=target[k:]))
 for control_seed in range(2026100501,2026100521):
  rng=np.random.default_rng(np.random.SeedSequence([seed,control_seed]));co=np.where(valid,cost+rng.uniform(0,.2,cost.shape),1e9);ii,jj=linear_sum_assignment(co);ok=len(ii)==kmax and valid[ii,jj].all();all_status.append(dict(seed=seed,control_seed=control_seed,status='MATCHED' if ok else 'DATA_INSUFFICIENT_MATCHING',n_pairs=int(valid[ii,jj].sum()),needed=kmax))
  if not ok:continue
  match=[pool[j] for j in jj];assert np.array_equal(ii,np.arange(kmax)) and len(set(match))==kmax
  for i,j in zip(ii,jj):all_pairs.append(dict(seed=seed,control_seed=control_seed,high_target=high[i],control_target=pool[j],high_order=int(i+1),expression_rank_diff=float(diff[i,j,0]),detection_rank_diff=float(diff[i,j,1]),log2_degree_diff=float(diff[i,j,2]),distance=float(cost[i,j])))
  for frac in [.1,.2]:
   k=math.ceil(frac*len(target));deleted=match[:k];specs.append(dict(name=f'Control_{control_seed}_{frac}',kind='MatchedControl',fraction=frac,control_seed=control_seed,deleted=deleted,retained=[t for t in target if t not in set(deleted)]))
 for sp in specs:
  for gene in target:all_members.append(dict(seed=seed,spec=sp['name'],kind=sp['kind'],fraction=sp['fraction'],control_seed=sp['control_seed'],target=gene,weight=weights[gene],deleted=gene in sp['deleted']))
 # One ranking per original seed; signatures change, expression and gene universe never change.
 ranks=create_rankings(expression,seed=seed);signatures=[GeneSignature(name=sp['name'],gene2weight={g:weights[g] for g in sp['retained']}) for sp in specs];auc=aucell4r(ranks,signatures,auc_threshold=.05,noweights=True,normalize=False,num_workers=1);auc=auc.loc[expression.index]
 oldfile=P/f'02_analysis/scenic2_20261001_092655_ZKNWD6/objects/Human_All_{seed}_DropOverlap_auc.csv.gz';head=pd.read_csv(oldfile,nrows=0);oldcell=pd.read_csv(oldfile,usecols=[head.columns[0],'IKZF2(+)'],index_col=0).loc[expression.index,'IKZF2(+)'];cellerr=float(np.max(abs(auc.Baseline-oldcell)));assert cellerr<1e-9,(seed,cellerr)
 dauc=auc.groupby(md['sample']).mean().loc[samples];oldmeans=baseold[baseold.seed==seed].set_index('sample').loc[samples].value;derr=float(np.max(abs(dauc.Baseline-oldmeans)));assert derr<1e-9
 baseline=cor(dauc.Baseline);oldcor=rd(G/'IKZF2_matched_reference_contrast_models.tsv');oldcor=oldcor[(oldcor.variant=='DropOverlap')&(oldcor.seed==seed)&(oldcor.endpoint=='excess')].iloc[0];assert abs(baseline['rho']-oldcor.rho)<1e-9 and abs(baseline['P']-oldcor.P)<1e-9
 checks.append(dict(seed=seed,cells=len(md),genes=len(genes),targets=len(target),max_cell_AUC_error=cellerr,max_donor_AUC_error=derr,rho=baseline['rho'],P=baseline['P'],rho_error=baseline['rho']-oldcor.rho,P_error=baseline['P']-oldcor.P))
 auc.to_csv(R/f'fixed_AUC_seed{seed}_SERVER_ONLY.tsv.gz',sep='\t')
 for sp in specs:
  name=sp['name'];nret=len(sp['retained']);x=dauc[name].values;fit=cor(x);delta=fit['rho']-baseline['rho'];fixed=auc[name]*nret/len(target);assert np.max(fixed-auc.Baseline)<1e-10
  jk=[]
  if sp['kind']=='HighWeight':
   for j,s in enumerate(samples):
    mask=np.arange(16)!=j;ff=cor(x[mask],meta.iloc[mask]);bb=cor(dauc.Baseline.values[mask],meta.iloc[mask]);v=ff['rho']-bb['rho'];jk.append(v);loo.append(dict(seed=seed,fraction=sp['fraction'],omitted=s,rho=ff['rho'],baseline_rho=bb['rho'],delta_rho=v,status=ff['status']))
  sejk=float(np.sqrt(15/16*np.sum((np.asarray(jk)-np.mean(jk))**2))) if len(jk)==16 and np.isfinite(jk).all() else np.nan
  z=dict(seed=seed,spec=name,kind=sp['kind'],fraction=sp['fraction'],control_seed=sp['control_seed'],n=16,n_control=7,n_case=9,targets_before=len(target),deleted=len(sp['deleted']),retained=nret,deleted_set=';'.join(sorted(sp['deleted'])),**fit,delta_rho=delta,delta_abs_rho=abs(fit['rho'])-abs(baseline['rho']),delta_JK_low=delta-stats.t.ppf(.975,15)*sejk,delta_JK_high=delta+stats.t.ppf(.975,15)*sejk,LOO_valid=int(np.isfinite(jk).sum()) if jk else 0,AUC_donor_rank_agreement=stats.spearmanr(dauc.Baseline,x).statistic);all_models.append(z)
  for j,s in enumerate(samples):donors.append(dict(seed=seed,spec=name,kind=sp['kind'],fraction=sp['fraction'],control_seed=sp['control_seed'],sample=s,group=meta.group.iloc[j],series=meta.series.iloc[j],cells=int(meta.cells.iloc[j]),mean_depth=meta.mean_depth.iloc[j],excess=meta.y.iloc[j],AUC=x[j],baseline_AUC=dauc.Baseline.iloc[j],delta_AUC=x[j]-dauc.Baseline.iloc[j],baseline_denominator_AUC=x[j]*nret/len(target),retained_fraction=nret/len(target)))
  b=dauc.Baseline.values;raw=x*nret/len(target);mechanical.append(dict(seed=seed,spec=name,kind=sp['kind'],fraction=sp['fraction'],control_seed=sp['control_seed'],deleted_fraction=len(sp['deleted'])/len(target),mean_donor_delta_normalized_AUC=float(np.mean(x-b)),mean_donor_raw_contribution_loss=float(np.mean(b-raw)),mean_donor_fraction_contribution_lost=float(np.mean((b-raw)/b)),size_normalization_factor=len(target)/nret,mean_donor_rank_change=float(np.mean(abs(stats.rankdata(x)-stats.rankdata(b))))))
 print('SEED_COMPLETE',seed,'signatures',len(specs),'baseline',baseline,'controls',sum(z['kind']=='MatchedControl' for z in specs),flush=True)
models=pd.DataFrame(all_models);models['FDR']=np.nan;models['BH_family']='NOT_COMPUTED_DESCRIPTIVE_CONTROL_OR_REUSED_BASELINE';sel=models.kind=='HighWeight';assert sel.sum()==6;models.loc[sel,'FDR']=bh(models.loc[sel,'P']);models.loc[sel,'BH_family']='six_predeclared_high_weight_deletion_associations';wr(models,'all_fixed_association_models.tsv');wr(all_members,'all_target_deletion_members.tsv.gz');wr(all_pairs,'control_pair_matching.tsv.gz');wr(all_status,'control_matching_status.tsv');wr(donors,'donor_fixed_readouts_SERVER_ONLY.tsv.gz');wr(checks,'baseline_reproduction.tsv');wr(mechanical,'mechanical_scoring_diagnostics.tsv');wr(pd.concat(all_features),'target_matching_features.tsv');wr(loo,'high_weight_donor_delete_intervals.tsv')
summ=[]
for (seed,frac),q in models[models.kind=='MatchedControl'].groupby(['seed','fraction']):
 h=models[(models.seed==seed)&(models.fraction==frac)&(models.kind=='HighWeight')].iloc[0];summ.append(dict(seed=seed,fraction=frac,controls=len(q),distinct_deleted_sets=q.deleted_set.nunique(),control_delta_median=q.delta_rho.median(),control_delta_q025=q.delta_rho.quantile(.025),control_delta_q975=q.delta_rho.quantile(.975),control_delta_min=q.delta_rho.min(),control_delta_max=q.delta_rho.max(),high_delta=h.delta_rho,high_percentile=float((q.delta_rho<=h.delta_rho).mean()),fraction_controls_abs_change_ge_high=float((q.delta_rho.abs()>=abs(h.delta_rho)).mean()),high_outside_control_95=bool(h.delta_rho<q.delta_rho.quantile(.025) or h.delta_rho>q.delta_rho.quantile(.975))))
wr(summ,'matched_control_effect_distribution.tsv');info=dict(status='COMPLETE',signatures=len(models),high_models=int(sel.sum()),matched_controls=int((models.kind=='MatchedControl').sum()),donors=16,cells=1361,original_seeds=[17,29,43],spec_sha256=hashlib.sha256(spec.read_bytes()).hexdigest(),no_new_network=True);(O/'provenance/compute_complete.json').write_text(json.dumps(info,indent=2));print('HIGH_RESULTS',models[sel].to_string(index=False));print('CONTROL_DISTRIBUTIONS',pd.DataFrame(summ).to_string(index=False));print(json.dumps(info),flush=True)
