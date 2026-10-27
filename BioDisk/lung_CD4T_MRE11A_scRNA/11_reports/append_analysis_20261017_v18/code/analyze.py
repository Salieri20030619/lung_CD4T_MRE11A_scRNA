import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
from pathlib import Path
import pandas as pd,numpy as np,json,itertools,math
from scipy import stats
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/append_analysis_20261017_v18');plan=json.loads((O/'protocol_frozen.json').read_text());score=pd.read_csv(O/'tables/sample_gene_and_program_scores.tsv',sep='\t');auc=pd.read_csv(O/'tables/clean_TF_sample_AUC.tsv',sep='\t');sender=pd.read_csv(O/'tables/sender_ligand_expression.tsv',sep='\t');mods=pd.read_csv(O/'tables/module_member_audit.tsv',sep='\t');pcache={};results=[];effectrows=[];testids=set()
def design(d,model):
 cols=[np.ones(len(d))]
 if model.startswith('Joint'):
  lab=d.dataset.astype(str)+'_'+d.day.astype(str)
  cols.extend(pd.get_dummies(lab,drop_first=True).to_numpy(dtype=float).T)
 if model=='TimeDayAdjusted':cols.extend(pd.get_dummies(d.day.astype(str),drop_first=True).to_numpy(dtype=float).T)
 if 'Condition' in model:cols.append(d.group.eq('Silica').to_numpy(dtype=float))
 return np.column_stack(cols)
def blocks(d,model):
 if model.startswith('Joint'):v=(d.dataset.astype(str)+'_'+d.day.astype(str)).tolist()
 elif model=='TimeDayAdjusted':v=d.day.astype(str).tolist()
 else:v=['all']*len(d)
 if 'Condition' in model:v=[a+'_'+b for a,b in zip(v,d.group.astype(str))]
 return tuple(v)
def permutations(labels):
 if labels in pcache:return pcache[labels]
 n=len(labels);groups=[np.flatnonzero(np.array(labels)==s) for s in sorted(set(labels))];count=math.prod(math.factorial(len(g)) for g in groups)
 if count<=10000:
  out=np.tile(np.arange(n),(count,1));choices=[list(itertools.permutations(g)) for g in groups]
  for i,part in enumerate(itertools.product(*choices)):
   for g,values in zip(groups,part):out[i,g]=values
  exact=True
 else:
  rng=np.random.default_rng(20261017);out=np.tile(np.arange(n),(10000,1))
  for i in range(len(out)):
   for g in groups:out[i,g]=rng.permutation(g)
  exact=False
 pcache[labels]=(out,exact,count);return pcache[labels]
def coefficient(x,y,d,model):
 C=design(d,model);rank=np.linalg.matrix_rank(C)
 if rank!=C.shape[1] or len(d)-rank-1<1:return np.nan,None,None,None,None
 M=np.eye(len(d))-C@np.linalg.pinv(C);rx=M@stats.rankdata(x);ry=M@stats.rankdata(y);den=np.linalg.norm(rx)*np.linalg.norm(ry)
 if den<1e-12:return np.nan,M,rx,ry,rank
 return float(np.clip(np.dot(rx,ry)/den,-1,1)),M,rx,ry,rank

def assoc(d,xcol,ycol,cohort,variant,kind,scope,model,source='none',extra=None):
 key=(cohort,variant,kind,scope,model,source,xcol,ycol)
 if key in testids:return
 testids.add(key);base=dict(cohort=cohort,variant=variant,kind=kind,scope=scope,model=model,sender=source,x=xcol,y=ycol,family='|'.join([cohort,variant,model,kind,source]));base.update(extra or {})
 if xcol not in d or ycol not in d:results.append({**base,'status':'NOT_AVAILABLE','n':0});return
 mask=np.isfinite(d[xcol].to_numpy(float))&np.isfinite(d[ycol].to_numpy(float));z=d.loc[mask].reset_index(drop=True);x=z[xcol].to_numpy(float);y=z[ycol].to_numpy(float);base.update(n=len(z),sample_list=';'.join(z['sample'].astype(str)),dropped_samples=';'.join(d.loc[~mask,'sample'].astype(str)))
 if len(z)<4:results.append({**base,'status':'TOO_FEW_INDEPENDENT_SAMPLES'});return
 r,M,rx,ry,rank=coefficient(x,y,z,model)
 if not np.isfinite(r):results.append({**base,'status':'CONSTANT_OR_NO_RESIDUAL_DF'});return
 ids,exact,total=permutations(blocks(z,model));permy=stats.rankdata(y)[ids]@M;den=np.linalg.norm(permy,axis=1)*np.linalg.norm(rx);valid=den>1e-12;null=(permy[valid]@rx)/den[valid];extreme=int((np.abs(null)>=abs(r)-1e-12).sum());p=extreme/len(null) if exact else (extreme+1)/(len(null)+1);df=len(z)-rank-1;pt=2*stats.t.sf(abs(r)*np.sqrt(df/max(1-r*r,1e-15)),df);loo=[]
 for i in range(len(z)):
  take=np.arange(len(z))!=i;v=coefficient(x[take],y[take],z.loc[take].reset_index(drop=True),model)[0]
  if np.isfinite(v):loo.append(v)
 base.update(status='OK',rho=r,P=p,P_model_t=pt,residual_df=df,permutations_evaluated=len(null),permutations_total=total,permutation_exact_enumeration=exact,loo_valid=len(loo),loo_same_sign=float(np.mean(np.sign(loo)==np.sign(r))) if loo else np.nan,loo_rho_min=min(loo) if loo else np.nan,loo_rho_max=max(loo) if loo else np.nan);results.append(base)

def group_effect(d,col,cohort,variant,scope):
 r=dict(cohort=cohort,variant=variant,scope=scope,program=col[2:],family=cohort+'|'+variant+'|programme_disease_effect');v=d[col].to_numpy(float)
 if not np.isfinite(v).all() or len(np.unique(v))<2:effectrows.append({**r,'status':'CONSTANT_OR_MISSING'});return
 if cohort.startswith('274'):
  C=np.column_stack([np.ones(len(d)),pd.get_dummies(d.day.astype(str),drop_first=True).to_numpy(float)]);b=np.linalg.lstsq(C,v,rcond=None)[0];rss=np.sum((v-C@b)**2);rss0=np.sum((v-v.mean())**2);df1=np.linalg.matrix_rank(C)-1;df2=len(d)-np.linalg.matrix_rank(C);F=max(0,((rss0-rss)/df1)/(rss/df2)) if rss>1e-15 else np.nan;p=stats.f.sf(F,df1,df2) if np.isfinite(F) else np.nan
  effectrows.append({**r,'status':'OK','n':len(d),'P':p,'method':'time_omnibus_F_model','effect':np.nan,'F':F,'n_control':int((d.day==0).sum()),'n_case':int((d.day!=0).sum()),'sample_list':';'.join(d['sample'])});return
 cond=d.group.eq('Silica').to_numpy(float);C=design(d,'JointStratum' if cohort=='Joint' else 'Pooled');M=np.eye(len(d))-C@np.linalg.pinv(C);resid=M@v;xc=M@cond;beta=float(xc@resid/(xc@xc));labels=blocks(d,'JointStratum' if cohort=='Joint' else 'Pooled');groups=[np.flatnonzero(np.array(labels)==g) for g in sorted(set(labels))];choices=[list(itertools.combinations(g,int(cond[g].sum()))) for g in groups];null=[]
 for picked in itertools.product(*choices):
  xx=np.zeros(len(d))
  for ids in picked:xx[list(ids)]=1
  zz=M@xx
  if zz@zz>1e-12:null.append(float(zz@resid/(zz@zz)))
 p=float(np.mean(np.abs(null)>=abs(beta)-1e-12));effectrows.append({**r,'status':'OK','n':len(d),'n_control':int((cond==0).sum()),'n_case':int((cond==1).sum()),'effect':beta,'P':p,'method':'stratified_assignment_permutation','allocations':len(null),'sample_list':';'.join(d['sample'])})

ph=plan['phenotypes'];damage=plan['damage'];aux=plan['auxiliary'];panels=ph+damage+aux
# All planned receiver associations; no result-driven state selection.
for variant in ['Raw','BG1','FixedComp','BG1FixedComp']:
 for cohort in ['314','Joint','274','274_min10']:
  if cohort.startswith('274') and variant!='Raw':continue
  if cohort=='Joint' and variant in ['BG1','BG1FixedComp']:continue
  for scope in (['TLike'] if cohort.startswith('274') else ['Core','Naive','ActCD4','Treg']):
   if variant in ['FixedComp','BG1FixedComp'] and scope!='Core':continue
   ds=['GSE314943'] if cohort=='314' else ['GSE314943','GSE183682'] if cohort=='Joint' else ['GSE274372'];d=score[score.dataset.isin(ds)&score.variant.eq(variant)&score.scope.eq(scope)].sort_values(['dataset','sample']).reset_index(drop=True)
   if cohort=='274_min10':d=d[d.cells>=10].reset_index(drop=True)
   if not len(d):continue
   models=['JointStratum','JointStratumCondition'] if cohort=='Joint' else ['Pooled','Condition'] if cohort=='314' else ['TimePooled','TimeDayAdjusted']
   gg=sorted(set(g for g,sc in plan['focus_pairs'] if sc==scope)) if scope!='TLike' else sorted(set(g for g,sc in plan['focus_pairs']))
   for model in models:
    for g in gg:
     for mods0,kind in [(ph,'candidate_phenotype'),(damage,'candidate_damage'),(aux,'candidate_auxiliary')]:
      for mod in mods0:assoc(d,'G_'+g,'M_'+mod,cohort,variant,kind,scope,model,extra={'gene':g,'program':mod})
     if g=='IL1R2':assoc(d,'G_IL1R2','G_ENO1',cohort,variant,'candidate_auxiliary',scope,model,extra={'gene':g,'program':'ENO1_RNA'})
    for p in ph:
     for m in damage:
      a=set(mods.loc[(mods.module==p)&~mods.removed_candidate,'gene']);b=set(mods.loc[(mods.module==m)&~mods.removed_candidate,'gene']);assert not a&b,('overlapping_pheno_damage',p,m,a&b)
      assoc(d,'M_'+p,'M_'+m,cohort,variant,'phenotype_damage',scope,model,extra={'phenotype':p,'program':m})
   for mod in panels:group_effect(d,'M_'+mod,cohort,variant,scope)

# Source/candidate, source/proximal, and candidate/proximal steps use the same six mice and receiver state.
for variant in ['Raw','BG1','FixedComp','BG1FixedComp']:
 for gene,scope,lig,tfs in plan['routes']:
  if variant in ['FixedComp','BG1FixedComp'] and scope!='Core':continue
  d=score[(score.dataset=='GSE314943')&(score.variant==variant)&(score.scope==scope)].copy();aa=auc[(auc.variant==variant)&(auc.scope==scope)].pivot(index='sample',columns='tf',values='AUC').add_prefix('TF_');d=d.merge(aa,on='sample',how='left')
  for send in ['AllMac','AMlike']:
   ss=sender[(sender.sender==send)&(sender.ligand==lig)][['sample','CP10k']].rename(columns={'CP10k':'L_'+lig});z=d.merge(ss,on='sample',how='left');z['L_'+lig]=np.log1p(z['L_'+lig])
   for model in ['Pooled','Condition']:
    assoc(z,'L_'+lig,'G_'+gene,'314',variant,'ligand_candidate',scope,model,send,{'gene':gene,'ligand':lig})
    if tfs not in ['none','decoy']:
     for tf in tfs.split(';'):
      assoc(z,'L_'+lig,'TF_'+tf,'314',variant,'ligand_proximalTF',scope,model,send,{'ligand':lig,'tf':tf})
      assoc(z,'TF_'+tf,'G_'+gene,'314',variant,'proximalTF_candidate',scope,model,'none',{'gene':gene,'tf':tf})
# Target-specific TF-function compatibility, not asserted as a unique receptor pathway.
for gene,scope,tfs in [('FOXP1','Core',['FOXP1']),('FOXP1','Treg',['FOXP1']),('SKIL','Treg',['TRP53']),('AHNAK','Naive',['NFATC1','NFATC2'])]:
 for var in ['Raw','BG1','FixedComp','BG1FixedComp']:
  if var in ['FixedComp','BG1FixedComp'] and scope!='Core':continue
  d=score[(score.dataset=='GSE314943')&(score.variant==var)&(score.scope==scope)].copy();aa=auc[(auc.variant==var)&(auc.scope==scope)].pivot(index='sample',columns='tf',values='AUC').add_prefix('TF_');d=d.merge(aa,on='sample',how='left')
  for model in ['Pooled','Condition']:
   for tf in tfs:
    assoc(d,'G_'+gene,'TF_'+tf,'314',var,'candidate_focalTF',scope,model,extra={'gene':gene,'tf':tf})
    for mod in ph+damage:assoc(d,'TF_'+tf,'M_'+mod,'314',var,'focalTF_program',scope,model,extra={'tf':tf,'program':mod})

def finalize(rows,name):
 d=pd.DataFrame(rows);d['FDR']=np.nan;d['family_planned_n']=0;d['family_evaluable_n']=0
 for fam,ids in d.groupby('family').groups.items():
  ok=d.loc[ids,'P'].notna();ii=np.array(list(ids))[ok.to_numpy()];d.loc[ids,'family_planned_n']=len(ids);d.loc[ids,'family_evaluable_n']=len(ii)
  if len(ii):d.loc[ii,'FDR']=stats.false_discovery_control(d.loc[ii,'P'].to_numpy())
 d.to_csv(O/'tables'/name,sep='\t',index=False);return d
r=finalize(results,'all_associations.tsv');e=finalize(effectrows,'programme_group_effects.tsv')
# Reference directions for n1+1 timepoints in the second dataset; no p-values invented.
desc=[]
for var in ['Raw','FixedComp']:
 for (scope,day),d in score[(score.dataset=='GSE183682')&(score.variant==var)].groupby(['scope','day']):
  if set(d.group)!=set(['Saline','Silica']) or len(d)!=2:continue
  for col in [c for c in d if c.startswith(('G_','M_'))]:desc.append(dict(variant=var,scope=scope,day=day,variable=col,delta=float(d.loc[d.group=='Silica',col].iloc[0]-d.loc[d.group=='Saline',col].iloc[0]),n_control=1,n_case=1,P=np.nan,status='DESCRIPTIVE_ONLY'))
pd.DataFrame(desc).to_csv(O/'tables/183682_timepoint_descriptive.tsv',sep='\t',index=False)
for name,d in [('association',r),('effect',e)]:print(name,len(d),'families',d.family.nunique(),'evaluable',d.P.notna().sum(),'nominal',(d.P<.05).sum(),'FDR',(d.FDR<.05).sum(),flush=True)
(O/'logs/association_complete.json').write_text(json.dumps({'status':'success','associations':len(r),'association_families':int(r.family.nunique()),'effects':len(e),'effect_families':int(e.family.nunique())},indent=2),encoding='utf-8')