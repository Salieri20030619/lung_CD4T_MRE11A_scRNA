from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import pandas as pd,numpy as np,json
from scipy import stats
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/human_first_20261004_v1';R=O/'results'
assert (O/'provenance/stage2_complete.txt').exists()
def rd(p):return pd.read_csv(p,sep='\t')
def wr(d,n):pd.DataFrame(d).to_csv(R/n,sep='\t',index=False,na_rep='NA')
def bh(p):
 p=np.asarray(p);i=np.argsort(p);v=p[i]*len(p)/np.arange(1,len(p)+1);o=np.empty(len(p));o[i]=np.minimum(1,np.minimum.accumulate(v[::-1])[::-1]);return o
def correlate(d,x,y,extra=None):
 cols=[x,y,'sender_depth','receiver_depth']+([extra] if extra else []);q=d[np.isfinite(d[cols]).all(axis=1)&(d.sender_cells>=20)&(d.receiver_cells>=20)].copy();n=len(q)
 out=dict(n=n,n_case=int(q.group.eq('CWP').sum()),n_control=int(q.group.eq('Control').sum()),samples=';'.join(q['sample']),rho=np.nan,P=np.nan,df=np.nan,status='INSUFFICIENT')
 if min(out['n_case'],out['n_control'])<2:return out
 X=np.column_stack([np.ones(n),q.group.eq('CWP'),pd.get_dummies(q.series,drop_first=True,dtype=float),stats.rankdata(q.sender_depth),stats.rankdata(q.receiver_depth)])
 if extra:X=np.column_stack([X,stats.rankdata(q[extra])])
 rank=np.linalg.matrix_rank(X);df=n-rank-1;out['df']=df
 if rank<X.shape[1] or df<3:out['status']='NON_IDENTIFIABLE_OR_DF_LT3';return out
 if q[x].nunique()<2 or q[y].nunique()<2:out['status']='CONSTANT';return out
 a=stats.rankdata(q[x]);b=stats.rankdata(q[y]);a-=X@np.linalg.lstsq(X,a,rcond=None)[0];b-=X@np.linalg.lstsq(X,b,rcond=None)[0]
 if min(np.std(a),np.std(b))<1e-10:out['status']='NO_RESIDUAL_VARIATION';return out
 rho=np.clip(np.corrcoef(a,b)[0,1],-1,1);p=2*stats.t.sf(abs(rho)*np.sqrt(df/max(1e-15,1-rho*rho)),df);out.update(rho=rho,P=p,status='EXPLORATORY_SAME_DATA_ASSOCIATION');return out
expr=rd(R/'stage2_LR_expression_all_donors.tsv.gz');prog=rd(R/'stage2_disease_program_donor_scores.tsv.gz');meta=rd(R/'stage2_scope_sample_counts.tsv');cand=rd(R/'stage2_all_candidates_with_negative_results.tsv.gz');design=rd(R/'stage1_donor_design.tsv')
de_status=rd(R/'stage2_DE_status.tsv');set_status=rd(R/'stage2_targetset_status.tsv');extra=[];extra_status=[]
for _,z in de_status[de_status.status.ne('FIT')].iterrows():
 for test in ['FDR05','Nominal05']:
  for direction in ['Up','Down']:
   st=dict(receiver=z.receiver,filter=z['filter'],model=z.model,test=test,direction=direction,background_n=np.nan,geneset_n=np.nan,candidate_n=np.nan,status='DE_'+z.status)
   if len(set_status[(set_status.receiver==z.receiver)&(set_status['filter']==z['filter'])&(set_status.model==z.model)&(set_status.test==test)&(set_status.direction==direction)]):continue
   extra_status.append(st)
   for sen in ['MacAM','MonoMac','MacWide']:
    v=pd.DataFrame({'ligand':sorted(cand.ligand.unique())});v['sender']=sen
    for k,val in st.items():
     if k!='candidate_n':v[k]=val
    v['LR_expression_supported']=False;v['aupr_corrected']=np.nan;v['rank']=np.nan;extra.append(v)
if extra:
 cand=pd.concat([cand]+extra,ignore_index=True);wr(cand,'stage2_all_candidates_with_negative_results.tsv.gz');wr(pd.concat([set_status,pd.DataFrame(extra_status)],ignore_index=True),'stage2_targetset_status.tsv')
tnk=rd(P/'11_reports/missing_links_20261003_v1/results/sender_all_samples.tsv');tnk=tnk[(tnk.species=='Human')&(tnk.sender=='TNK_other')];tnk_cols=[x for x in tnk if x.startswith('L_') and 'complex' not in x]
rows=[];comp=[]
keys=['receiver','filter','model','test','direction']
for key,g in prog.groupby(keys):
 rec,filt,model,test,direction=key;rr=meta[(meta.lineage=='CD4')&(meta.scope==rec)&(meta['filter']==filt)][['sample','cells','mean_depth']].rename(columns={'cells':'receiver_cells','mean_depth':'receiver_depth'})
 base=g[['sample','target_score']].merge(design[['sample','group','series']],on='sample').merge(rr,on='sample')
 cs=cand
 for k,v in zip(keys,key):cs=cs[cs[k]==v]
 for sen,v in cs[cs.LR_expression_supported].groupby('sender'):
  sm=meta[(meta.lineage=='Mac')&(meta.scope==sen)][['sample','cells','mean_depth']].rename(columns={'cells':'sender_cells','mean_depth':'sender_depth'})
  sx=expr[expr.key==sen].pivot(index='sample',columns='gene',values='mean_CP10k')
  d=base.merge(sm,on='sample').merge(sx,left_on='sample',right_index=True)
  for lig in v.ligand.unique():
   r=dict(zip(keys,key));r.update(sender=sen,ligand=lig,**correlate(d,lig,'target_score'));rows.append(r)
   c=dict(zip(keys,key));c.update(sender=sen,ligand=lig)
   col='L_'+lig
   if col not in tnk_cols:c.update(n=0,rho=np.nan,P=np.nan,status='TNK_LIGAND_NOT_IN_EXISTING_EXPORT')
   else:
    tt=tnk[['sample',col,'sender_cells']].rename(columns={col:'TNK_ligand','sender_cells':'TNK_cells'});dd=d.merge(tt,on='sample');dd=dd[dd.TNK_cells>=20];c.update(correlate(dd,lig,'target_score',extra='TNK_ligand'))
   comp.append(c)
for rows_,name in [(rows,'stage2_same_donor_associations.tsv.gz'),(comp,'stage2_TNK_competition_sensitivity.tsv.gz')]:
 d=pd.DataFrame(rows_);d['FDR']=np.nan
 for _,ix in d.groupby(keys).groups.items():
  ix=[i for i in ix if np.isfinite(d.loc[i,'P'])];d.loc[ix,'FDR']=bh(d.loc[ix,'P'])
 wr(d,name);print(name,'rows',len(d),'computed',int(d.P.notna().sum()),'BH05',int((d.FDR<.05).sum()),flush=True)
(O/'provenance/stage2_followup_complete.json').write_text(json.dumps({'associations':len(rows),'competition':len(comp),'inference':'same-data exploratory; no independent validation'}))
