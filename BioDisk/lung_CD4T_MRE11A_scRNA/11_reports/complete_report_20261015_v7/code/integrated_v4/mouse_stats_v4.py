import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
from pathlib import Path
import sys,itertools,json
import numpy as np,pandas as pd
from scipy.stats import rankdata,false_discovery_control
O=Path(sys.argv[1]);MODS=['PyroExec','Inflam','ApoExec','OXPHOS','MitoISR','MitoROS','Mitophagy'];CAND=['IKZF2','STAT1','IL1R2','PENK','TNFRSF9','LAG3','TTC39C','FOXP1','AHNAK','MDFIC','ST6GALNAC3'];SAMPLES=['Saline_1','Saline_2','Saline_3','Silica_1','Silica_2','Silica_3'];GROUP=np.array([0,0,0,1,1,1])
def rd(n):return pd.read_csv(O/n,sep='\t')
def wr(a,n):a.to_csv(O/n,sep='\t',index=False,encoding='utf-8',na_rep='NA')
def bh(p):
 p=np.asarray(p,float);q=np.full(p.shape,np.nan);ok=np.isfinite(p)
 if ok.any():q[ok]=false_discovery_control(p[ok])
 return q
PERM=np.array(list(itertools.permutations(range(6))))
WITHIN=np.array([list(a)+list(b) for a in itertools.permutations(range(3)) for b in itertools.permutations(range(3,6))])
def centered_rank(a,g,adjust):
 z=rankdata(a,axis=0).astype(float)
 if adjust:
  for group in np.unique(g):z[g==group]-=z[g==group].mean(axis=0)
 else:z-=z.mean(axis=0)
 return z

def corr(x,y,g,adjust):
 rx=centered_rank(x,g,adjust);ry=centered_rank(y,g,adjust);den=np.sqrt((rx*rx).sum(axis=0))[:,None]*np.sqrt((ry*ry).sum(axis=0))[None,:];den[den<1e-12]=np.nan
 return rx.T@ry/den,rx,ry,den

def assoc(X,Y,scope,source,seed,kind):
 assert list(X.index)==SAMPLES and list(Y.index)==SAMPLES
 assert np.isfinite(X.to_numpy()).all() and np.isfinite(Y.to_numpy()).all()
 out=[]
 for adjust in [False,True]:
  model='ConditionAdjusted' if adjust else 'Pooled';r,rx,ry,den=corr(X.to_numpy(),Y.to_numpy(),GROUP,adjust);perms=WITHIN if adjust else PERM
  null=np.einsum('psa,sb->pab',rx[perms],ry,optimize=True)/den[None,:,:];p=(np.abs(null)>=np.abs(r)[None,:,:]-1e-12).mean(axis=0);p[~np.isfinite(r)]=np.nan
  loo=[]
  for drop in range(6):
   keep=np.arange(6)!=drop;loo.append(corr(X.to_numpy()[keep],Y.to_numpy()[keep],GROUP[keep],adjust)[0])
  loo=np.asarray(loo);valid=np.isfinite(loo);nv=valid.sum(axis=0);agree=(valid&(np.sign(loo)==np.sign(r)[None,:,:])).sum(axis=0);frac=np.divide(agree,nv,out=np.full_like(r,np.nan),where=nv>0)
  n=len(X.columns)*len(Y.columns)
  frame=pd.DataFrame({'kind':kind,'scope':scope,'sender':source,'seed':seed,'model':model,'x':np.repeat(X.columns,len(Y.columns)),'y':np.tile(Y.columns,len(X.columns)),'rho':r.ravel(),'P':p.ravel(),'n_samples':6,'permutation_allocations':len(perms),'loo_valid':nv.ravel(),'loo_same_sign':frac.ravel()})
  frame['family']=f'{kind}_{scope}_{source}_{seed}_{model}';frame['BH']=bh(frame.P);out.append(frame)
 return pd.concat(out,ignore_index=True)

phen=rd('01_phenotype/sample_program_scores.tsv');phen=phen[phen['filter']=='All'];act=rd('03_receiver/TF_all_regulon_sample_activity.tsv.gz');rec=rd('03_receiver/receptor_candidate_sample_expression.tsv.gz');send=rd('02_sender/all_ligand_sample_expression.tsv.gz');elig=rd('02_sender/ligand_scope_eligibility.tsv');tfresults=[];rnaprogram=[];group=[];pairs=[]
alloc=np.array(list(itertools.combinations(range(6),3)))
for scope in ['Core','Naive','ActCD4','Treg']:
 Y=phen[phen.scope==scope].set_index('sample').loc[SAMPLES,MODS];R=rec[(rec.scope==scope)&rec.gene.isin(CAND)].pivot(index='sample',columns='gene',values='CP10k').reindex(index=SAMPLES,columns=CAND);assert R.notna().all().all()
 rnaprogram.append(assoc(R,Y,scope,'none',0,'RNA_program'))
 for seed in [17,29,43]:
  A=act[(act.scope==scope)&(act.seed==seed)].pivot(index='sample',columns='regulon',values='mean_auc').loc[SAMPLES]
  tfresults.append(assoc(A,Y,scope,'none',seed,'TF_AUC_program'))
  arr=A.to_numpy();delta=arr[3:].mean(axis=0)-arr[:3].mean(axis=0);null=[]
  for ix in alloc:
   mask=np.zeros(6,bool);mask[ix]=True;null.append(arr[~mask].mean(axis=0)-arr[mask].mean(axis=0))
  pp=(np.abs(np.asarray(null))>=np.abs(delta)[None,:]-1e-12).mean(axis=0);pp[np.ptp(arr,axis=0)<1e-12]=np.nan
  z=pd.DataFrame({'scope':scope,'seed':seed,'regulon':A.columns,'delta_mean_auc':delta,'control_mean_auc':arr[:3].mean(axis=0),'case_mean_auc':arr[3:].mean(axis=0),'P':pp,'BH':bh(pp),'n_control':3,'n_case':3,'family':f'AUC_disease_{scope}_{seed}_all_regulons'});group.append(z)
  if seed==17:
   focused=A[[c for c in A.columns if c.rsplit('(',1)[0].upper() in ['IKZF2','STAT1']]].copy();focused.columns=['AUC::'+c for c in focused.columns];rna=R.copy();rna.columns=['RNA::'+c for c in rna.columns];endpoint=pd.concat([focused,rna],axis=1)
   for source in ['AMlike','AllMac']:
    ls=elig[(elig.sender==source)&(elig.scope==scope)&elig.expression_supported].ligand.to_list();X=send[(send.scope==source)&send.gene.isin(ls)].pivot(index='sample',columns='gene',values='CP10k').reindex(index=SAMPLES,columns=ls)
    assert X.notna().all().all()
    pairs.append(assoc(X,endpoint,scope,source,17,'sender_receiver'))
 print('SCOPE_DONE',scope,flush=True)
tfresults=pd.concat(tfresults,ignore_index=True);rnaprogram=pd.concat(rnaprogram,ignore_index=True);pairs=pd.concat(pairs,ignore_index=True);group=pd.concat(group,ignore_index=True)
wr(tfresults,'03_receiver/all_TF_program_associations.tsv.gz');wr(rnaprogram,'03_receiver/candidate_RNA_program_associations.tsv');wr(pairs,'04_chains/all_sender_receiver_associations.tsv.gz');wr(group,'03_receiver/all_TF_activity_group_effects.tsv')
focus=tfresults[tfresults.x.str.split('(').str[0].str.upper().isin(['IKZF2','STAT1'])];wr(focus,'03_receiver/IKZF2_STAT1_program_associations.tsv')
# Fixed-seed exact estimates, condition adjustment, and LOO remain distinct.
primary=focus[(focus.seed==17)&(focus.model=='Pooled')].copy();adj=focus[(focus.seed==17)&(focus.model=='ConditionAdjusted')][['scope','x','y','rho','P','BH','loo_same_sign','loo_valid']].rename(columns={c:'condition_'+c for c in ['rho','P','BH','loo_same_sign','loo_valid']});primary=primary.merge(adj,on=['scope','x','y'],validate='1:1')
for seed in [29,43]:
 b=focus[(focus.seed==seed)&(focus.model=='Pooled')][['scope','x','y','rho','P']].rename(columns={'rho':f'rho_seed{seed}','P':f'P_seed{seed}'});primary=primary.merge(b,on=['scope','x','y'],how='left',validate='1:1')
wr(primary,'03_receiver/IKZF2_STAT1_primary_and_sensitivity.tsv')
# Expression coverage robustness uses the historical10% criterion only as an explicit flag.
lr=rd('inputs/LR_prior.tsv');rob=[]
for source in ['AMlike','AllMac']:
 sm=send[(send.scope==source)&(send.group=='Silica')].groupby('gene').pct.apply(lambda x:int((x>=10).sum()))
 for scope in ['Core','Naive','ActCD4','Treg']:
  rm=rec[(rec.scope==scope)&(rec.group=='Silica')].groupby('gene').pct.apply(lambda x:int((x>=10).sum()))
  a=lr.copy();a['sender']=source;a['scope']=scope;a['source10_libraries']=a.ligand.map(sm).fillna(0).astype(int);a['receiver10_libraries']=a.receptor.map(rm).fillna(0).astype(int);a['legacy10pct_supported']=(a.source10_libraries>=2)&(a.receiver10_libraries>=2);rob.append(a)
wr(pd.concat(rob),'02_sender/LR_legacy10pct_robustness.tsv.gz')
families=[]
for name,a in [('TF_program',tfresults),('RNA_program',rnaprogram),('sender_receiver',pairs),('TF_group',group)]:
 for family,b in a.groupby('family'):families.append(dict(type=name,family=family,records=len(b),valid_tests=int(b.P.notna().sum()),nominal_P05=int((b.P<.05).sum()),BH05=int((b.BH<.05).sum())))
wr(pd.DataFrame(families),'provenance/BH_families.tsv')
summary={'status':'PASS','TF_program_rows':len(tfresults),'RNA_program_rows':len(rnaprogram),'sender_receiver_rows':len(pairs),'TF_group_rows':len(group),'biological_samples':6,'cell_level_P_tests':0,'primary_seed':17,'secondary_seeds_are_not_replicates':True,'matched_mouse_roster':SAMPLES}
(O/'provenance/mouse_stats_complete.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(json.dumps(summary));print('FOCUSED_POOLED_HITS');print(primary[primary.P<.05][['scope','x','y','rho','P','BH','condition_rho','condition_P','loo_same_sign']].to_string(index=False))