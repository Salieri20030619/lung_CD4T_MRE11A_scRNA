import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
from pathlib import Path
import json,itertools
import numpy as np,pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/death_priority_20261018_v2');N=O/'joint_network';F=O/'figures'
def rd(n):return pd.read_csv(O/'tables'/n,sep='\t')
checks=[]
def check(n,ok):checks.append({'check':n,'passed':bool(ok)});assert ok,n
def bh(p):
 p=np.array(p,float);n=len(p);i=np.argsort(p);v=np.minimum.accumulate((p[i]*n/np.arange(1,n+1))[::-1])[::-1];o=np.empty(n);o[i]=np.minimum(v,1);return o
for arm in ['All','BG1']:
 audit=pd.read_csv(N/'Mouse'/arm/'training_audit.tsv',sep='\t')
 for seed in [17,29,43]:
  z=audit[audit.seed==seed];check(f'{arm} {seed} training10samples',len(z)==10 and z.cells.between(1,250).all())
  for stage in ['grn','ctx']:check(f'{arm} {seed} {stage} exit',json.loads((N/'Mouse'/arm/f'seed{seed}'/(stage+'_exit.json')).read_text())['returncode']==0)
for name in ['joint_network_TF_associations.tsv','joint_network_TF_disease_effects.tsv','joint_network_upstream_TF_associations.tsv']:
 d=rd(name)
 for fam,z in d.groupby('family'):
  z=z[z.P.notna()]
  if len(z):check(name+' BH '+fam,np.allclose(bh(z.P),z.FDR,atol=1e-12,rtol=0))
 check(name+' n10',d.loc[d.P.notna(),'n'].eq(10).all())
auc=rd('joint_network_TF_sample_AUC.tsv');score=rd('sample_scores.tsv');aa=rd('joint_network_TF_associations.tsv');cache={};err=[]
for r in aa[(aa.variant=='Raw')&aa.kind.eq('jointnet_TF_candidate_primary')&aa.P.notna()].itertuples():
 d=score[(score.variant==r.variant)&(score.scope==r.scope)&score.dataset.ne('GSE274372')].merge(auc[(auc.variant==r.variant)&(auc.scope==r.scope)&(auc.tf==r.tf)][['dataset','sample','AUC']],on=['dataset','sample']).sort_values(['dataset','sample']).reset_index(drop=True)
 x=stats.rankdata(d.AUC);y=stats.rankdata(d['G_'+r.gene]);lab=d.dataset+'_'+d.day.astype(str);C=np.column_stack([np.ones(len(d)),pd.get_dummies(lab,drop_first=True,dtype=float)])
 if 'Condition' in r.model:C=np.column_stack([C,(d.group=='Silica').astype(float)]);lab=lab+'_'+d.group
 Q=np.eye(len(d))-C@np.linalg.pinv(C);rx=Q@x;ry=Q@y;rho=rx@ry/(np.linalg.norm(rx)*np.linalg.norm(ry));key=tuple(lab)
 if key not in cache:
  blocks=[np.flatnonzero(lab.to_numpy()==v) for v in sorted(lab.unique())];ids=[]
  for parts in itertools.product(*[list(itertools.permutations(b)) for b in blocks]):
   ix=np.arange(len(d))
   for b,v in zip(blocks,parts):ix[b]=v
   ids.append(ix)
  cache[key]=np.array(ids)
 vals=y[cache[key]]@Q;den=np.linalg.norm(vals,axis=1)*np.linalg.norm(rx);null=vals[den>1e-12]@rx/den[den>1e-12];p=float(np.mean(abs(null)>=abs(rho)-1e-12));err.append(dict(gene=r.gene,scope=r.scope,tf=r.tf,model=r.model,rho_error=abs(rho-r.rho),P_error=abs(p-r.P)))
num=pd.DataFrame(err,columns=['gene','scope','tf','model','rho_error','P_error']);check('joint primary TF independent rho/P',len(num)==0 or (num.rho_error.max()<1e-10 and num.P_error.max()<1e-12));num.to_csv(O/'tables/joint_network_independent_QA.tsv',sep='\t',index=False)
pd.DataFrame(checks).to_csv(O/'tables/joint_network_QA_checks.tsv',sep='\t',index=False)
cand=rd('frozen_candidate_pairs.tsv');edges=rd('joint_network_candidate_edges.tsv');routes=rd('joint_network_candidate_TF_routes.tsv');old=rd('TF_associations.tsv');rows=[]
for r in cand.itertuples():
 for var in ['Raw','BG1']:
  z=aa[(aa.variant==var)&aa.kind.eq('jointnet_TF_candidate_'+r.branch)&(aa.gene==r.gene)&(aa.scope==r.scope)&(aa.model=='JointStratum')];w=z[z.P<.05].sort_values('P')
  rows.append(dict(gene=r.gene,scope=r.scope,branch=r.branch,variant=var,tested=len(z),nominal=len(w),FDR=int((z.FDR<.05).sum()),signals=';'.join(f'{v.tf}:rho={v.rho:.4g},P={v.P:.4g},FDR={v.FDR:.4g}' for v in w.itertuples())))
summ=pd.DataFrame(rows);summ.to_csv(O/'tables/joint_network_candidate_summary.tsv',sep='\t',index=False)
coverage=rd('joint_network_TF_coverage.tsv');run=rd('joint_network_run_summary.tsv');effects=rd('joint_network_TF_disease_effects.tsv')
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':190})
figs=[]
def save(fid,title,source):
 plt.savefig(F/(fid+'.png'),bbox_inches='tight',facecolor='white');plt.close();figs.append(dict(id=fid,title=title,source=source,file=str(F/(fid+'.png'))))
fig,axs=plt.subplots(1,2,figsize=(11,4))
for ax,arm in zip(axs,['All','BG1']):
 z=run[run.arm==arm].sort_values('seed');ax.bar(z.seed.astype(str),z.regulons,color='#287c8e' if arm=='All' else '#bd7252')
 for i,v in enumerate(z.regulons):ax.text(i,v,str(v),ha='center',va='bottom')
 scored=int(((coverage.arm==arm)&(coverage.status=='SCORED')).sum());ax.set(title=f'{arm}: joint10 donor discovery\nStable, cleaned >=5-target regulons: {scored}',xlabel='Random seed',ylabel='Motif-supported regulons')
plt.tight_layout();save('JN01','联合10样本三种子重新建网与稳定regulon','joint_network_run_summary.tsv; joint_network_TF_coverage.tsv')
focal=['HIF1A','FOXP1','ATF4','TRP53','STAT1','STAT3','STAT5A','STAT5B','NFKB1','NFKB2','RELA','RELB','JUN','JUNB','FOS','FOSB','NFATC1','NFATC2','FOXO1','FOXO3','TFEB']
z=effects[(effects.variant=='Raw')&effects.scope.isin(['Naive','ActCD4','Treg'])];focal=[g for g in focal if g in set(z.program)]
mat=z.pivot(index='scope',columns='program',values='effect').reindex(index=['Naive','ActCD4','Treg'],columns=focal);pv=z.pivot(index='scope',columns='program',values='P').reindex_like(mat)
fig,ax=plt.subplots(figsize=(12,3.5));v=mat.to_numpy();lim=np.nanmax(abs(v));im=ax.imshow(v,cmap='RdBu_r',vmin=-lim,vmax=lim,aspect='auto');ax.set(xticks=range(len(focal)),xticklabels=focal,yticks=range(3),yticklabels=mat.index);plt.setp(ax.get_xticklabels(),rotation=45,ha='right')
for i in range(3):
 for j in range(len(focal)):
  if pv.iloc[i,j]<.05:ax.text(j,i,'*',ha='center',va='center',fontsize=15)
fig.colorbar(im,ax=ax,label='Silica - Saline, stable-regulon AUC');plt.tight_layout();save('JN02','联合新网络的TF活性疾病效应','joint_network_TF_disease_effects.tsv')
# Structural prior is drawn independently of association significance.
chosen=[]
for gene in ['HIF1A','IL2RA','ARL6IP1']:
 z=edges[(edges.arm=='All')&edges.stable_2of3&(edges.target==gene)].sort_values(['seed_support','weight'],ascending=False).head(5)
 chosen.extend(z.to_dict('records'))
sel=pd.DataFrame(chosen);sel.to_csv(O/'tables/joint_network_display_edges.tsv',sep='\t',index=False)
tfs=sorted(set(sel.tf)) if len(sel) else [];genes=['HIF1A','IL2RA','ARL6IP1'];fig,ax=plt.subplots(figsize=(11,6));ty={tf:i for i,tf in enumerate(tfs)};gy={g:v for g,v in zip(genes,np.linspace(0,max(1,len(tfs)-1),3))}
for r in sel.itertuples():
 q=aa[(aa.variant=='Raw')&(aa.model=='JointStratum')&(aa.gene==r.target)&(aa.tf==r.tf)&aa.kind.eq('jointnet_TF_candidate_primary')];p=q.P.iloc[0] if len(q) else np.nan
 ax.plot([0,1],[ty[r.tf],gy[r.target]],color='#bd6042' if p<.05 else '#b9c6ca',lw=1.2 if r.seed_support==2 else 2.2,zorder=1)
for tf,yy in ty.items():ax.scatter(0,yy,s=90,c='#287c8e',zorder=3);ax.text(-.04,yy,tf,ha='right',va='center')
for gene,yy in gy.items():ax.scatter(1,yy,s=180,c='#d6a257',zorder=3);ax.text(1.04,yy,gene,ha='left',va='center')
ax.set(xlim=(-.35,1.4),ylim=(-1,max(1,len(tfs))),title='Joint All: stable motif edges (top 5 by support/weight per candidate)\nRed: matched-state RNA/AUC association raw P < 0.05; grey: otherwise/NA');ax.axis('off');plt.tight_layout();save('JN03','联合网络的候选调控先验与样本关联并列','joint_network_candidate_edges.tsv; joint_network_TF_associations.tsv')
(F/'manifest_joint_network.json').write_text(json.dumps(figs,ensure_ascii=False,indent=2),encoding='utf-8')
status={'status':'success','checks':len(checks),'independent_primary_TF_tests':len(num),'runs':6,'All_scored':int(((coverage.arm=='All')&(coverage.status=='SCORED')).sum()),'BG1_scored':int(((coverage.arm=='BG1')&(coverage.status=='SCORED')).sum()),'primary_raw_nominal':int(summ[(summ.branch=='primary')&(summ.variant=='Raw')].nominal.sum())}
(N/'QA_complete.json').write_text(json.dumps(status,indent=2),encoding='utf-8')
print(json.dumps(status));print(summ.to_string(index=False))
