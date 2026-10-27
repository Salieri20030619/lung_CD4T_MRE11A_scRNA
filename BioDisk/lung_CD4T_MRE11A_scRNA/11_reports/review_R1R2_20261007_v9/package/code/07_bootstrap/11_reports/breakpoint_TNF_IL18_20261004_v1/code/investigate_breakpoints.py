from pathlib import Path
import os,sys,json,hashlib,datetime,fcntl,pickle,itertools
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np,pandas as pd
from scipy import stats
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/breakpoint_TNF_IL18_20261004_v1';F=P/'11_reports/front_half_20261004_v1';M=P/'11_reports/missing_links_20261003_v1';A=P/'11_reports/aligned_evidence_20261002_v1'
lock=open(O/'analysis.lock','a+');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
assert not (O/'provenance/complete.json').exists()
manifest={}
def rd(p):
 p=Path(p)
 if str(p) not in manifest:
  h=hashlib.sha256()
  with p.open('rb') as f:
   for b in iter(lambda:f.read(2**20),b''):h.update(b)
  manifest[str(p)]=dict(path=str(p),bytes=p.stat().st_size,sha256=h.hexdigest())
 return pd.read_csv(p,sep='\t',low_memory=False)
def wr(d,n):pd.DataFrame(d).to_csv(O/'results'/n,sep='\t',index=False,na_rep='NA')
pair=rd(F/'results/paired_inputs.tsv')
covorder=['G','S','receiver_library','sender_library','sender_epi','Treg_fraction','TNK_IFNG','IFN_state_NoRegUnion']
def parse(s):
 if not isinstance(s,str) or s in ['None','']:return []
 d={'R':'receiver_library','D':'sender_library','E':'sender_epi','T':'TregNoTF'}
 return [d.get(c,c) for c in s.split('+')]
def cols(cs):return ['group' if c=='G' else 'series' if c=='S' else c for c in cs]
def design(z,cs):
 D=[np.ones(len(z))];names=['Intercept']
 for c in sorted(cs,key=lambda c:covorder.index(c) if c in covorder else 100+sum(map(ord,c))):
  if c=='G':D.append(z.group.eq('CWP').to_numpy(float));names.append('G')
  elif c=='S':
   for v in sorted(z.series.unique())[1:]:D.append(z.series.eq(v).to_numpy(float));names.append('S_'+str(v))
  else:
   v=stats.rankdata(z[c]);sd=v.std();D.append((v-v.mean())/sd if sd>0 else np.zeros(len(v)));names.append(c)
 return np.column_stack(D),names
def signature(z,x,y,cs):
 D,names=design(z,cs)
 return hashlib.sha256(json.dumps([z['sample'].tolist(),stats.rankdata(z[x]).tolist(),stats.rankdata(z[y]).tolist(),np.round(D,10).tolist()]).encode()).hexdigest()
def cache_row(r,origin):
 return dict(rho=r.get('rho',np.nan),P=r.get('P',r.get('p',np.nan)),legacy_FDR=r.get('FDR',np.nan),legacy_id=str(r.get('model_id',r.get('unique_test_id',''))),origin=origin,LOO_valid=r.get('LOO_valid',np.nan),LOO_min=r.get('LOO_min',r.get('LOO_rho_min',np.nan)),LOO_max=r.get('LOO_max',r.get('LOO_rho_max',np.nan)),LOO_same_direction=r.get('LOO_same_direction',np.nan),P_method=r.get('P_method',r.get('p_method','partial_rank_t_exploratory')))
def build_cache():
 cache={};audit=[]
 for root,modelname,inputnames in [(F,'new_link_models.tsv',['new_link_model_inputs.tsv.gz']),(M,'all_new_models_final.tsv',['model_sample_inputs.tsv.gz','matched_modality_sample_inputs.tsv.gz'])]:
  models=rd(root/'results'/modelname);models=models[models.status.eq('FIT')]
  if 'species' in models:models=models[models.species.eq('Human')]
  inputs=pd.concat([rd(root/'results'/n) for n in inputnames],ignore_index=True)
  mi={k:v.sort_values('sample') for k,v in inputs.groupby('model_id')}
  used=0
  for _,r in models.iterrows():
   if r.model_id not in mi:continue
   z=mi[r.model_id];cs=parse(r.covariates)
   if len(z)!=r.n or not all(c in z for c in cols(cs)) or z[['x_value','y_value']+cols(cs)].isna().any().any():continue
   key=signature(z,'x_value','y_value',cs);cache.setdefault(key,cache_row(r,root.name));used+=1
  audit.append(dict(source=root.name,eligible_existing_rows=used,unique_cache_total=len(cache)))
 # Existing aligned MacAM IL18--IKZF2 AUC, using its own saved input values.
 models=rd(A/'results/aligned_all_models.tsv.gz');old=rd(A/'inputs/aligned_receiver_candidates_sender.tsv')
 models=models[(models.species=='Human')&(models['filter']=='All')&(models.evidence=='SENDER_TF')&(models.tf=='IKZF2')&(models.status=='FIT')]
 old=old[(old.species=='Human')&(old['filter']=='All')&(old.tf=='IKZF2')]
 used=0
 for _,r in models.iterrows():
  cs=parse(r.covariates)
  if any(c not in ['G','S','receiver_library','sender_library','sender_epi'] for c in cs):continue
  z=old[(old.scope==r.scope)&(old.seed==r.seed)&old['sample'].isin(str(r.sample_ids).split(';'))].drop_duplicates('sample').sort_values('sample')
  if len(z)!=r.n:continue
  current=pair[(pair.sender=='MacAM')&(pair.scope==r.scope)].set_index('sample').reindex(z['sample'])
  assert np.allclose(z.mean_auc,current['IKZF2_AUC_'+str(r.seed)],rtol=1e-10,atol=1e-12)
  assert np.allclose(z.sender_mean_CP10k,current.L_IL18,rtol=1e-10,atol=1e-12)
  key=signature(z,'sender_mean_CP10k','mean_auc',cs);cache.setdefault(key,cache_row(r,'aligned_evidence'));used+=1
 audit.append(dict(source='aligned_evidence',eligible_existing_rows=used,unique_cache_total=len(cache)));wr(audit,'reuse_cache_inventory.tsv')
 return cache
def generate():
 specs=[];data={};rosters=[];roster_summary=[]
 for (sender,scope),z in pair.groupby(['sender','scope']):
  z=z.sort_values('sample').reset_index(drop=True)
  for axis,tf in [('TNF','STAT1'),('IL18','IKZF2')]:
   chain=f'{sender}/{scope}/{axis}-{tf}'
   prox=['NFkB_proxy_NoRegUnion']+(['AP1_components_NoRegUnion'] if axis=='IL18' else [])
   commoncols=['L_'+axis,tf+'_RNA',tf+'_AUC_17','State_Effector','State_Activation']+prox+['receiver_library','sender_library','TNK_IFNG']
   qr=z.dropna(subset=[c for c in commoncols if c!=tf+'_AUC_17']);qa=z.dropna(subset=commoncols)
   for name,q in [('RNA19',qr),('AUC17',qa),('AUC_cells20',qa[(qa.sender_cells>=20)&(qa.receiver_cells>=20)])]:
    roster_summary.append(dict(chain=chain,roster=name,n=len(q),n_case=int(q.group.eq('CWP').sum()),n_control=int(q.group.eq('Control').sum()),sample_ids=';'.join(q['sample']),excluded_ids=';'.join(sorted(set(z['sample'])-set(q['sample'])))))
    for _,v in z.iterrows():
     reason=[]
     if name!='RNA19' and pd.isna(v[tf+'_AUC_17']):reason.append('frozen_AUC_unavailable_upstream_Core_lt20')
     if name=='AUC_cells20':
      if v.sender_cells<20:reason.append('sender_cells_lt20')
      if v.receiver_cells<20:reason.append('receiver_cells_lt20')
     if v['sample'] not in set(q['sample']) and not reason:reason.append('required_endpoint_missing')
     rosters.append(dict(chain=chain,roster=name,sample=v['sample'],group=v.group,series=v.series,sender_cells=v.sender_cells,receiver_cells=v.receiver_cells,included=v['sample'] in set(q['sample']),reason=';'.join(reason)))
   for modality,q in [('RNA',qr),('AUC17',qa)]:
    tcol=tf+'_RNA' if modality=='RNA' else tf+'_AUC_17'
    edges=[('ligand_proximal','L_'+axis,y) for y in prox]+[('ligand_TF','L_'+axis,tcol)]+[('TF_phenotype',tcol,y) for y in ['State_Effector','State_Activation']]+[('ligand_phenotype','L_'+axis,y) for y in ['State_Effector','State_Activation']]
    base=[[],['G'],['G','S'],['G','S','receiver_library'],['G','S','receiver_library','sender_library']];last=base[-1]
    def add(q0,edge,x,y,cs,policy,seed=17):
     zz=q0.dropna(subset=[x,y]+cols(cs)).sort_values('sample').copy();sid='BP'+str(len(specs)+1).zfill(5)
     fam='primary' if sender=='MacAM' and policy=='aligned_steps' else 'sensitivity'
     meta=dict(spec_id=sid,chain=chain,sender=sender,scope=scope,axis=axis,tf=tf,modality=modality,seed=seed if 'AUC' in x or 'AUC' in y else np.nan,edge=edge,x=x,y=y,covariates='+'.join(cs) or 'None',policy=policy,FDR_family=fam,n=len(zz),sample_ids=';'.join(zz['sample']),signature=signature(zz,x,y,cs))
     specs.append(meta);data[sid]=zz
    for edge,x,y in edges:
     for cs in base:add(q,edge,x,y,cs,'aligned_steps')
     if edge!='ligand_proximal':
      if scope=='Core':add(q,edge,x,y,last+['Treg_fraction'],'Treg_composition')
      add(q,edge,x,y,last+['TNK_IFNG'],'TNK_competition')
     add(q[(q.sender_cells>=20)&(q.receiver_cells>=20)],edge,x,y,last,'cells20')
     if axis=='TNF' and edge in ['ligand_TF','ligand_phenotype']:add(q,edge,x,y,last+['IFN_state_NoRegUnion'],'IFN_environment')
     if modality=='RNA' and edge in ['ligand_TF','TF_phenotype']:
      for cs in [[],last]:add(qa,edge,x,y,cs,'RNA_matched_AUC17')
     if modality=='AUC17' and edge in ['ligand_TF','TF_phenotype']:
      for seed in [29,43]:
       for cs in [[],last]:add(qa,edge,x.replace('_17','_'+str(seed)),y.replace('_17','_'+str(seed)),cs,'seed_sensitivity',seed)
 return pd.DataFrame(specs),data,pd.DataFrame(rosters),pd.DataFrame(roster_summary)
if '--plan' in sys.argv:
 assert not (O/'provenance/frozen_plan.pkl').exists()
 cache=build_cache();specs,data,rosters,summary=generate();specs['already_available']=specs.signature.isin(cache)
 wr(specs,'frozen_specification_inventory.tsv');wr(rosters,'donor_roster_reasons.tsv');wr(summary,'chain_roster_summary.tsv');wr(manifest.values(),'input_manifest.tsv')
 plan=dict(time=datetime.datetime.now().isoformat(),specification_rows=len(specs),unique_signatures=specs.signature.nunique(),reused_specifications=int(specs.already_available.sum()),new_unique_signatures=specs.loc[~specs.already_available,'signature'].nunique(),new_specifications=int((~specs.already_available).sum()),chains=specs.chain.nunique(),no_new_fits_yet=True)
 with (O/'provenance/frozen_plan.pkl').open('wb') as f:pickle.dump((cache,specs,data),f)
 (O/'provenance/plan_inventory.json').write_text(json.dumps(plan,indent=2));print(json.dumps(plan,indent=2));print(summary.to_string(index=False));sys.exit(0)
assert '--run' in sys.argv and not (O/'provenance/models_complete.json').exists()
with (O/'provenance/frozen_plan.pkl').open('rb') as f:cache,specs,data=pickle.load(f)
# Reuse existing aligned per-donor LOO instead of recalculating it.
legacy_models=rd(A/'results/aligned_all_models.tsv.gz').drop_duplicates('unique_test_id').set_index('unique_test_id')
legacy_loo=rd(A/'results/aligned_LOO.tsv.gz');legacy_loo=legacy_loo[(legacy_loo.species=='Human')&(legacy_loo['filter']=='All')&(legacy_loo.evidence=='SENDER_TF')&(legacy_loo.tf=='IKZF2')].copy();legacy_loo['covkey']=legacy_loo.covariates.fillna('None')
loo_lookup={key:v for key,v in legacy_loo.groupby(['scope','seed','covkey','sample_ids'])}
for v in cache.values():
 if v['origin']!='aligned_evidence' or v['legacy_id'] not in legacy_models.index:continue
 r=legacy_models.loc[v['legacy_id']];q=loo_lookup.get((r.scope,r.seed,r.covariates if isinstance(r.covariates,str) else 'None',r.sample_ids))
 if q is not None:
  q=q.drop_duplicates('omitted_sample');values=q.loo_rho.dropna();v.update(LOO_valid=len(values),LOO_min=values.min() if len(values) else np.nan,LOO_max=values.max() if len(values) else np.nan,LOO_same_direction=bool(len(values) and all(np.sign(a)==np.sign(v['rho']) for a in values)))
def diagnostics(z,x,y,cs):
 D,names=design(z,cs);rank=np.linalg.matrix_rank(D);df=len(z)-rank-1;rx=stats.rankdata(z[x]);ry=stats.rankdata(z[y]);a=rx-D@np.linalg.lstsq(D,rx,rcond=None)[0];b=ry-D@np.linalg.lstsq(D,ry,rcond=None)[0]
 vif=(np.sum((rx-rx.mean())**2)/np.sum(a*a)) if np.sum(a*a)>1e-15 else np.inf
 vv=[]
 for j in range(1,D.shape[1]):
  other=np.delete(D,j,axis=1);res=D[:,j]-other@np.linalg.lstsq(other,D[:,j],rcond=None)[0];ss=np.sum((D[:,j]-D[:,j].mean())**2);vv.append(ss/np.sum(res*res) if np.sum(res*res)>1e-15 else np.inf)
 out=dict(n=len(z),df=df,design_rank=rank,design_columns=';'.join(names),condition_number=float(np.linalg.cond(D)),endpoint_VIF=float(vif),max_covariate_VIF=max(vv) if vv else 1.,x_unique=z[x].nunique(),y_unique=z[y].nunique(),x_zero=int(z[x].eq(0).sum()),y_zero=int(z[y].eq(0).sum()),nominal_critical_abs_rho=float(stats.t.ppf(.975,df)/np.sqrt(stats.t.ppf(.975,df)**2+df)) if df>0 else np.nan,status='FIT',small_n_low_df=len(z)<15 or df<=4)
 if len(z)<4:out['status']='TOO_FEW_FOR_DESCRIPTION'
 elif rank<D.shape[1]:out['status']='RANK_DEFICIENT'
 elif df<3:out['status']='RESIDUAL_DF_BELOW3'
 elif min(np.linalg.norm(a),np.linalg.norm(b))<1e-10:out['status']='CONSTANT_ENDPOINT_OR_RESIDUAL'
 elif len(z)<10:out['status']='DESCRIPTIVE_N_LT10'
 return out,a,b
def calc(z,x,y,cs):
 d,a,b=diagnostics(z,x,y,cs);r=dict(rho=np.nan,P=np.nan,**d)
 if d['status'] in ['FIT','DESCRIPTIVE_N_LT10']:
  r['rho']=float(np.clip(a@b/np.linalg.norm(a)/np.linalg.norm(b),-1,1))
  if d['status']=='FIT':r['P']=float(2*stats.t.sf(abs(r['rho'])*np.sqrt(d['df']/max(1-r['rho']**2,1e-15)),d['df']))
 return r
rows=[];loolist=[];inputs=[];newcache={}
for _,s in specs.iterrows():
 z=data[s.spec_id];cs=parse(s.covariates);dg,_,_=diagnostics(z,s.x,s.y,cs)
 if s.signature in cache:
  out=cache[s.signature].copy();out.update(dg);out['calculation']='reused_historical';out['LOO_detail']='historical_summary_only'
  if dg['status']!='FIT':out['P']=np.nan
 elif s.signature in newcache:out=newcache[s.signature].copy();out['calculation']='reused_new_signature'
 else:
  out=calc(z,s.x,s.y,cs);out.update(origin='breakpoint_missing_alignment',legacy_id='',legacy_FDR=np.nan,calculation='new_missing_signature',P_method='partial_rank_t_exploratory',LOO_detail='per_donor_new_diagnostic')
  lv=[]
  if np.isfinite(out['rho']):
   for ss in z['sample']:
    v=calc(z[z['sample'].ne(ss)],s.x,s.y,cs);loolist.append(dict(signature=s.signature,omitted_sample=ss,**v))
    if np.isfinite(v['rho']):lv.append(v['rho'])
  out.update(LOO_valid=len(lv),LOO_min=min(lv) if lv else np.nan,LOO_max=max(lv) if lv else np.nan,LOO_same_direction=bool(lv and all(np.sign(v)==np.sign(out['rho']) for v in lv)))
  newcache[s.signature]=out.copy()
 if np.isfinite(out['rho']) and dg['df']>1:
  se=1/np.sqrt(dg['df']-1);zz=np.arctanh(np.clip(out['rho'],-.999999,.999999));out.update(CI95_low=np.tanh(zz-1.96*se),CI95_high=np.tanh(zz+1.96*se))
 rows.append(dict(**s.to_dict(),**{k:v for k,v in out.items() if k!='n'}))
 for _,v in z.iterrows():
  rec=dict(spec_id=s.spec_id,signature=s.signature,sample=v['sample'],group=v.group,series=v.series,x_value=v[s.x],y_value=v[s.y],sender_cells=v.sender_cells,receiver_cells=v.receiver_cells)
  for c in cols(cs):rec[c]=v[c]
  inputs.append(rec)
m=pd.DataFrame(rows);m['FDR']=np.nan
for family,ix in m.groupby('FDR_family').groups.items():
 q=m.loc[ix];u=q[q.P.notna()].drop_duplicates('signature');values=stats.false_discovery_control(u.P.to_numpy(),method='bh');mapping=dict(zip(u.signature,values));m.loc[ix,'FDR']=m.loc[ix,'signature'].map(mapping)
m['nominal_P_lt05']=m.P.lt(.05);m['FDR_lt05']=m.FDR.lt(.05)
wr(m,'aligned_chain_models.tsv');wr(loolist,'new_signature_LOO.tsv.gz');wr(inputs,'aligned_model_inputs.tsv.gz')
# Within-source diagnostics: only predeclared four focal edges, baseline/group, no small-n P.
within=[]
for _,s in specs[(specs.policy=='aligned_steps')&(specs.covariates=='None')].iterrows():
 if s.edge=='TF_phenotype' and s.y!='State_Effector':continue
 if s.edge=='ligand_phenotype' and s.y!='State_Effector':continue
 if s.edge=='ligand_proximal' and s.y!='NFkB_proxy_NoRegUnion':continue
 for source,z in data[s.spec_id].groupby('series'):
  for cs in ([[],['G']] if z.group.nunique()>1 else [[]]):
   out=calc(z,s.x,s.y,cs);lv=[]
   if np.isfinite(out['rho']):
    for ss in z['sample']:
     v=calc(z[z['sample'].ne(ss)],s.x,s.y,cs)
     if np.isfinite(v['rho']):lv.append(v['rho'])
   within.append(dict(chain=s.chain,modality=s.modality,edge=s.edge,x=s.x,y=s.y,source=source,covariates='+'.join(cs) or 'None',sample_ids=';'.join(z['sample']),n_case=int(z.group.eq('CWP').sum()),n_control=int(z.group.eq('Control').sum()),LOO_valid=len(lv),LOO_min=min(lv) if lv else np.nan,LOO_max=max(lv) if lv else np.nan,**out))
w=pd.DataFrame(within);w['FDR']=np.nan;valid=w.P.notna()
if valid.any():w.loc[valid,'FDR']=stats.false_discovery_control(w.loc[valid,'P'].to_numpy(),method='bh')
wr(w,'within_source_descriptions.tsv')
existing_manifest=pd.read_csv(O/'results/input_manifest.tsv',sep='\t');wr(pd.concat([existing_manifest,pd.DataFrame(manifest.values())]).drop_duplicates('path'),'input_manifest.tsv')
summary=dict(time=datetime.datetime.now().isoformat(),specifications=len(m),chains=m.chain.nunique(),status=m.status.value_counts().to_dict(),calculation=m.calculation.value_counts().to_dict(),unique_new_fits=len(newcache),within_source_rows=len(w),no_old_grid_rerun=True)
(O/'provenance/models_complete.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
