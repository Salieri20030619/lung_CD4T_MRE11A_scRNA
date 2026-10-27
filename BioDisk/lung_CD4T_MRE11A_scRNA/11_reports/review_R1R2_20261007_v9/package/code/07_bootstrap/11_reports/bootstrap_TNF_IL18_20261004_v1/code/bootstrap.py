"""Finite donor-paired rank bootstrap. Run --preflight, then --run; no P tests."""
from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import argparse, json, hashlib, datetime, platform, itertools
import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata

B=Path(__file__).resolve().parents[1]; I=B/'inputs'; O=B/'results'; P=B/'provenance'
for p in [I,O,P]: p.mkdir(exist_ok=True)
N=5000; EDGES=[('ligand_TF',0,1),('TF_phenotype',1,2),('ligand_phenotype',0,2)]
AXES={'IL18':('CD4-IL7R','IKZF2'),'TNF':('Core','STAT1')}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def wr(df,name):pd.DataFrame(df).to_csv(O/name,sep='\t',index=False,na_rep='NA')
def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
pair=pd.read_csv(I/'paired_inputs.tsv',sep='\t');old=pd.read_csv(I/'aligned_chain_models.tsv',sep='\t')
rosters={}
for roster,mod in [('RNA19','RNA'),('MATCH17','AUC17')]:
 ids=[]
 for axis,(scope,tf) in AXES.items():
  q=old[(old.chain==f'MacAM/{scope}/{axis}-{tf}')&(old.modality==mod)&(old.policy=='aligned_steps')]
  ss=q.iloc[0].sample_ids.split(';');ids.append(ss)
 assert ids[0]==ids[1] and len(ids[0])==int(roster[-2:])
 rosters[roster]=ids[0]
def subset(axis,roster):
 scope,tf=AXES[axis]
 return pair[(pair.sender=='MacAM')&(pair.scope==scope)].set_index('sample').loc[rosters[roster]].reset_index()
def fit(z,axis,adjust):
 tf=AXES[axis][1]; names=['L_'+axis,tf+'_RNA']+[tf+'_AUC_'+str(s) for s in [17,29,43]]+['State_Effector']
 if len(z)==19:names=['L_'+axis,tf+'_RNA','State_Effector']
 V=z[names].to_numpy(float); ranks=rankdata(V,axis=0,method='average')
 D=[np.ones(len(z))];dn=['Intercept'];zero=[]
 if adjust!='unadjusted':
  D += [z.group.eq('CWP').to_numpy(float),z.series.eq('G227').to_numpy(float),z.series.eq('G288').to_numpy(float)]
  dn += ['G','S_G227','S_G288']
  for c in ['receiver_library','sender_library']+(['Treg_fraction'] if adjust=='B_Treg' else []):
   a=rankdata(z[c].to_numpy(float));sd=a.std();D.append((a-a.mean())/sd if sd>0 else a*0);dn.append(c)
 D=np.column_stack(D)
 for j in range(1,len(dn)):
  if D[:,j].std()<1e-12:zero.append(dn[j])
 rank=int(np.linalg.matrix_rank(D));res=ranks-D@np.linalg.lstsq(D,ranks,rcond=None)[0]
 norm=np.linalg.norm(res,axis=0);sd=ranks.std(axis=0)
 counts=z['sample'].value_counts().to_numpy();u=len(counts);df=len(z)-rank-1
 dg=dict(n=len(z),unique_donors=u,multiplicity_effective_n=float(len(z)**2/np.sum(counts**2)),design_columns=D.shape[1],design_rank=rank,nominal_df=df,unique_support_df=u-rank-1,low_unique_support=u-rank-1<3,condition_number=float(np.linalg.cond(D)),constant_covariates=';'.join(zero),zero_variance_endpoints=';'.join(np.array(names)[sd<1e-12]),zero_residual_endpoints=';'.join(np.array(names)[norm<1e-10]))
 results=[]
 for tfcol in names[1:-1]:
  modality='RNA19' if len(z)==19 else ('RNA17' if tfcol.endswith('_RNA') else 'AUC'+tfcol.rsplit('_',1)[1])
  for edge,a,b in EDGES:
   ix=[0,names.index(tfcol),len(names)-1];j,k=ix[a],ix[b];reasons=[]
   if not np.isfinite(V[:,[j,k]]).all() or not np.isfinite(D).all():reasons.append('NONFINITE')
   if rank<D.shape[1]:reasons.append('RANK_DEFICIENT')
   if df<3:reasons.append('NOMINAL_DF_BELOW3')
   if min(sd[j],sd[k])<1e-12:reasons.append('ZERO_VARIANCE_ENDPOINT')
   if min(norm[j],norm[k])<1e-10:reasons.append('ZERO_RESIDUAL_ENDPOINT')
   rho=float(np.clip(res[:,j]@res[:,k]/(norm[j]*norm[k]),-1,1)) if not reasons else np.nan
   results.append(dict(axis=axis,modality=modality,adjustment=adjust,edge=edge,x=names[j],y=names[k],rho=rho,status='VALID' if not reasons else ';'.join(reasons)))
 return results,dg
def original():
 rows=[];diagnostics=[]
 for roster in rosters:
  for axis in AXES:
   for adj in ['unadjusted','B']+(['B_Treg'] if axis=='TNF' else []):
    r,d=fit(subset(axis,roster),axis,adj);rows+=r;diagnostics.append(dict(axis=axis,roster=roster,adjustment=adj,**d))
 return pd.DataFrame(rows),pd.DataFrame(diagnostics)
def preflight():
 assert not (P/'FROZEN.json').exists()
 rows,dg=original();assert rows.rho.notna().all();checks=[]
 for r in rows.to_dict('records'):
  chain=f'MacAM/{AXES[r["axis"]][0]}/{r["axis"]}-{AXES[r["axis"]][1]}'
  cov={'unadjusted':'None','B':'G+S+receiver_library+sender_library','B_Treg':'G+S+receiver_library+sender_library+Treg_fraction'}[r['adjustment']]
  rr='RNA19' if r['modality']=='RNA19' else 'MATCH17';ids=';'.join(rosters[rr])
  q=old[(old.chain==chain)&(old.x==r['x'])&(old.y==r['y'])&(old.sample_ids==ids)&(old.covariates.fillna('None')==cov)]
  delta=float(np.max(np.abs(q.rho-r['rho']))) if len(q) else np.nan
  if len(q):assert delta<1e-10,(r,delta)
  checks.append(dict(**r,legacy_matches=len(q),max_abs_difference=delta,legacy_ids=';'.join(q.spec_id)))
 wr(rows,'original_effects.tsv');wr(dg,'original_diagnostics.tsv');wr(checks,'original_legacy_checks.tsv')
 cells=[]
 for roster in rosters:
  z=subset('TNF',roster);wr(z[['sample','group','series']],f'roster_{roster}.tsv')
  for (s,g),v in z.groupby(['series','group']):cells.append(dict(roster=roster,series=s,group=g,n=len(v),samples=';'.join(v['sample'])))
 wr(cells,'design_cells.tsv')
 inputs=[dict(path=str(p.relative_to(B)),sha256=sha(p),bytes=p.stat().st_size) for p in sorted(I.rglob('*')) if p.is_file()]
 wr(inputs,'input_hashes.tsv')
 freeze=dict(time=stamp(),plan_sha256=sha(B/'FROZEN_PLAN.txt'),code_sha256=sha(__file__),inputs=inputs,n_bootstrap=N,rosters=rosters,original_rows=len(rows),legacy_verified_rows=sum(x['legacy_matches']>0 for x in checks),environment=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__),status='FROZEN_BEFORE_RESAMPLING')
 (P/'FROZEN.json').write_text(json.dumps(freeze,indent=2),encoding='utf-8');print(json.dumps(freeze,indent=2));print(dg.to_string(index=False));print(pd.DataFrame(cells).to_string(index=False))
def run():
 frozen=json.loads((P/'FROZEN.json').read_text());assert frozen['plan_sha256']==sha(B/'FROZEN_PLAN.txt');assert frozen['code_sha256']==sha(__file__)
 for x in frozen['inputs']:assert sha(B/x['path'])==x['sha256']
 assert not (P/'COMPLETE.json').exists()
 lock=open(P/'run.lock','x');lock.write(str(os.getpid()));lock.flush()
 indices={};index_records=[]
 for scheme in ['stratified','unstratified']:
  for roster in rosters:
   seed={('stratified','RNA19'):2026100401,('stratified','MATCH17'):2026100402,('unstratified','RNA19'):2026100403,('unstratified','MATCH17'):2026100404}[(scheme,roster)]
   rng=np.random.Generator(np.random.PCG64(seed));z=subset('TNF',roster);n=len(z)
   if scheme=='unstratified':idx=rng.integers(0,n,(N,n))
   else:
    blocks=[]
    for _,v in z.groupby(['series','group'],sort=True):
     ix=v.index.to_numpy();blocks.append(rng.choice(ix,size=(N,len(ix)),replace=True))
    idx=np.concatenate(blocks,axis=1)
   indices[(scheme,roster)]=idx;np.save(O/f'indices_{scheme}_{roster}.npy',idx)
   wr(pd.DataFrame(idx,columns=[f'draw_{j+1}' for j in range(n)]).assign(replicate=np.arange(1,N+1)),f'indices_{scheme}_{roster}.tsv.gz')
   index_records.append(dict(scheme=scheme,roster=roster,seed=seed,sha256=sha(O/f'indices_{scheme}_{roster}.npy')))
 (P/'INDICES_FROZEN.json').write_text(json.dumps(dict(time=stamp(),indices=index_records),indent=2))
 rows=[];diagnostics=[];start=stamp()
 for (scheme,roster),idx in indices.items():
  data={axis:subset(axis,roster) for axis in AXES}
  for i,draw in enumerate(idx,1):
   for axis in AXES:
    z=data[axis].iloc[draw]
    for adj in ['unadjusted','B']+(['B_Treg'] if axis=='TNF' else []):
     r,d=fit(z,axis,adj);rows.extend(dict(scheme=scheme,replicate=i,roster=roster,**a) for a in r);diagnostics.append(dict(scheme=scheme,replicate=i,axis=axis,roster=roster,adjustment=adj,**d))
   if i%1000==0:print(scheme,roster,i,stamp(),flush=True)
  print('finished',scheme,roster,flush=True)
 wr(rows,'bootstrap_effects.tsv.gz');wr(diagnostics,'bootstrap_diagnostics.tsv.gz')
 for x in frozen['inputs']:assert sha(B/x['path'])==x['sha256']
 (P/'COMPLETE.json').write_text(json.dumps(dict(start=start,complete=stamp(),attempted_effects=len(rows),diagnostic_rows=len(diagnostics),inputs_unchanged=True),indent=2));lock.close()
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--preflight',action='store_true');ap.add_argument('--run',action='store_true');a=ap.parse_args()
 if a.preflight:preflight()
 elif a.run:run()
