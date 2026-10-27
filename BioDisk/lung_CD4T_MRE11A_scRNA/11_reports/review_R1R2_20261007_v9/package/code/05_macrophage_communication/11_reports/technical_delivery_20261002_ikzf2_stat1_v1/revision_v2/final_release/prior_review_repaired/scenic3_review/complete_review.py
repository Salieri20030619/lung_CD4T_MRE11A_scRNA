"""Run: python complete_review.py RESULT_ROOT OUTPUT_DIR. Uses returned tables only."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
import pandas as pd
from scipy import stats
ROOT=Path(sys.argv[1]);OUT=Path(sys.argv[2]);OUT.mkdir(exist_ok=True,parents=True)
sys.path.insert(0,str(ROOT/'scenic2_results'))
from stats import bh
manifest=[]
def rd(rel):
 p=ROOT/rel;manifest.append(dict(path=rel,sha256=hashlib.sha256(p.read_bytes()).hexdigest()));return pd.read_csv(p,sep='\t',low_memory=False)
def wr(x,name):pd.DataFrame(x).to_csv(OUT/name,sep='\t',index=False)
labels={0:('C1qAM','MacAM','C1QA/C1QB high, FABP4/MARCO positive'),1:('Mono','Excluded','FCN1/S100A8 high, C1Q/FABP4 low'),2:('AM_Dom','MacAM','C1Q/MARCO/PPARG high; CTL068 dominates'),3:('MonoMac','MonoMac','C1Q/MARCO with FCN1/SPP1; mixed differentiation state'),4:('Mono','Excluded','FCN1/VCAN/S100A8'),5:('Mono','Excluded','FCN1/S100A8/S100A9'),6:('DC2like','Excluded','CD1C/FCER1A/HLA-II high'),7:('MonoMac','MonoMac','SPP1/GPNMB/TREM2/MARCO with FCN1'),8:('EpiMix','Excluded','SFTPC/SCGB1A1/EPCAM; epithelial/myeloid mixed signal'),9:('TMix','Excluded','CD3D/CD3E/CCL5 with myeloid signal'),10:('VascMix','Excluded','VWF/CLDN5/RNASE1 with epithelial/myeloid signal'),11:('AM_Dom','MacAM','PPARG/MARCO/C1Q high; CTL070 dominates'),12:('DC1like','Excluded','CLEC9A/HLA-II high'),13:('AM_MARCO','MacAM','MARCO/FABP4/C1Q/PPARG high'),14:('Mono','Excluded','S100A8/S100A9/FCN1')}
wr([dict(cluster=k,label=v[0],scope=v[1],evidence=v[2]) for k,v in labels.items()],'reviewed_annotation.tsv')
c=rd('scenic3_results/tables/human_refined_cells.tsv.gz');assert set(c.r4)==set(labels) and c.cell.is_unique
c['label']=c.r4.map({k:v[0] for k,v in labels.items()});c['review_scope']=c.r4.map({k:v[1] for k,v in labels.items()})
assert not c[['label','review_scope']].isna().any().any()
wr(c[['cell','sample','group','series','parent','r2','r4','r8','label','review_scope']],'reviewed_cells.tsv.gz')
wr(c.groupby(['r4','label','sample','group','series'],observed=True).size().reset_index(name='cells'),'annotation_sample_counts.tsv')
ex=rd('scenic3_results/tables/human_sender_parent_sensitivity.tsv')
ligands=['IL15','TNF','IL18','CXCL16']
assert not ex.duplicated(['r4','sample','gene','sender_filter']).any()
sendrows=[];covrows=[]
scopes={'MacAM':[0,2,11,13],'MonoMac':[3,7],'MacWide':[0,2,3,7,11,13]}
for scope,clusters in scopes.items():
 for arm in ['All','NoParent16']:
  cc=c[c.r4.isin(clusters)&c['sample'].ne('Unknown')]
  if arm=='NoParent16':cc=cc[cc.parent!=16]
  for ss,z in cc.groupby('sample',observed=True):
   covrows.append(dict(sender=scope,sender_filter=arm,sample=ss,group=z.group.iloc[0],series=z.series.iloc[0],sender_cells=len(z),sender_library=z.library.mean(),sender_epi=z.epi_fraction.mean(),IL15_IL15RA_pct=100*z.IL15_IL15RA_det.mean()))
  ee=ex[ex.r4.isin(clusters)&ex.sender_filter.eq(arm)&ex['sample'].ne('Unknown')&ex.gene.isin(ligands)]
  for (ss,gene),z in ee.groupby(['sample','gene'],observed=True):
   n=z.cells.sum();assert n==len(cc[cc['sample']==ss])
   sendrows.append(dict(sender=scope,sender_filter=arm,sample=ss,ligand=gene,sender_mean_CP10k=np.average(z.mean_CP10k,weights=z.cells),sender_pct=np.average(z.pct,weights=z.cells)))
send=pd.DataFrame(sendrows).merge(pd.DataFrame(covrows),on=['sender','sender_filter','sample'],validate='many_to_one')
wr(send,'reviewed_sender_expression.tsv')
act=rd('scenic2_results/tables/TF_sample_activity.tsv.gz')
act=act[act.species.eq('Human')&act.variant.eq('DropOverlap')&act.regulon.isin(['STAT1(+)','IKZF2(+)','ETV7(+)'])]

def assoc(x,y,d,model):
 x=np.asarray(x,float);y=np.asarray(y,float);n=len(d)
 result=dict(n=n,rho=np.nan,p=np.nan,df=np.nan,reason='')
 if n<3:result['reason']='TOO_FEW';return result
 if len(np.unique(x))<2 or len(np.unique(y))<2:result['reason']='CONSTANT';return result
 cols=[np.ones(n)]
 if model!='Spearman':
  cols.append(d.group.eq('CWP').astype(float).to_numpy())
  cols.extend([(d.series==s).astype(float).to_numpy() for s in sorted(d.series.unique())[1:]])
  cols.append(stats.rankdata(d.receiver_library))
 if model=='BothEpi':cols.extend([stats.rankdata(d.sender_library),stats.rankdata(d.sender_epi)])
 X=np.column_stack(cols);rank=np.linalg.matrix_rank(X);df=n-rank-1;result['df']=df
 if rank<X.shape[1] or df<=0:result['reason']='UNIDENTIFIABLE';return result
 rx=stats.rankdata(x);ry=stats.rankdata(y);rx-=X@np.linalg.lstsq(X,rx,rcond=None)[0];ry-=X@np.linalg.lstsq(X,ry,rcond=None)[0]
 if np.linalg.norm(rx)<1e-10 or np.linalg.norm(ry)<1e-10:result['reason']='CONSTANT_RESIDUAL';return result
 rho=float(np.clip(np.corrcoef(rx,ry)[0,1],-1,1));result['rho']=rho
 if n>=10:result['p']=float(2*stats.t.sf(abs(rho)*np.sqrt(df/max(1-rho*rho,1e-15)),df))
 return result
res=[];paired=[]
for arm in ['All','BG1']:
 m=rd(f'scenic1_results/Human/{arm}/cells.tsv')
 for scope in ['Core','CD4-IL7R','CD4-GZMA','CD4-CCR7','Treg']:
  z=m if scope=='Core' else m[m.subtype.eq(scope)]
  md=z.groupby('sample',observed=True).agg(group=('group','first'),series=('series','first'),receiver_library=('library','mean'),receiver_cells=('cell','size')).reset_index()
  aa=act[act['filter'].eq(arm)&act.scope.eq(scope)]
  for (seed,tf),a in aa.groupby(['seed','regulon'],observed=True):
   d=a[['sample','mean_auc','cells']].merge(md,on='sample',validate='one_to_one');assert (d.cells==d.receiver_cells).all()
   for (sender,sarm,ligand),se in send.groupby(['sender','sender_filter','ligand'],observed=True):
    dd=d.merge(se,on='sample',validate='one_to_one',suffixes=('','_sender'))
    assert (dd.group==dd.group_sender).all() and (dd.series==dd.series_sender).all()
    dd=dd[(dd.receiver_cells>=20)&(dd.sender_cells>=20)].reset_index(drop=True)
    key=dict(filter=arm,scope=scope,seed=seed,tf=tf,sender=sender,sender_filter=sarm,ligand=ligand)
    paired.append(dd.assign(**key))
    for sample_set in ['Matched','NoDom']:
     dx=dd if sample_set=='Matched' else dd[~dd['sample'].isin(['CTL068','CTL070'])]
     dx=dx.reset_index(drop=True);x=dx.sender_mean_CP10k.to_numpy();y=dx.mean_auc.to_numpy()
     for model in ['Spearman','SourceDep','BothEpi']:
      r=assoc(x,y,dx,model);vals=[]
      for j in range(len(dx)):
       k=np.arange(len(dx))!=j;q=assoc(x[k],y[k],dx.iloc[k],model)
       if np.isfinite(q['rho']):vals.append(q['rho'])
      res.append(dict(**key,sample_set=sample_set,model=model,**r,loo_valid=len(vals),loo_min=min(vals) if vals else np.nan,loo_max=max(vals) if vals else np.nan))
r=pd.DataFrame(res);r['FDR']=np.nan
for _,idx in r.groupby(['filter','seed','sender_filter','sample_set','model'],observed=True).groups.items():
 valid=r.loc[idx].p.dropna().index
 if len(valid):r.loc[valid,'FDR']=bh(r.loc[valid,'p'])
wr(r,'reviewed_sender_TF_models.tsv');wr(pd.concat(paired),'aligned_samples.tsv.gz')
keys=['filter','scope','tf','sender','sender_filter','ligand','sample_set','model']
summary=r.groupby(keys,observed=True).agg(seeds=('seed','nunique'),n_min=('n','min'),rho_min=('rho','min'),rho_max=('rho','max'),p_min=('p','min'),p_max=('p','max'),FDR_min=('FDR','min'),FDR_max=('FDR','max'),loo_min=('loo_min','min'),loo_max=('loo_max','max')).reset_index();wr(summary,'seed_ranges.tsv')
rec=rd('scenic3_results/tables/CD4_receptor_detection.tsv');co=rd('scenic3_results/tables/CD4_receptor_codetection.tsv')
for tab,gene,name in [(rec,'gene','receptor_summary.tsv'),(co,'complex','receptor_complex_summary.tsv')]:
 tab=tab[tab.eligible20]
 sm=tab.groupby(['species','filter','scope',gene],observed=True).agg(samples=('sample','nunique'),pct_min=('pct','min'),pct_median=('pct','median'),pct_max=('pct','max')).reset_index();wr(sm,name)
# Descriptive sensitivity excludes FOXP3-detected ActCD4 cells without defining them as Tregs.
mo=rd('scenic3_results/tables/mouse_identity_cells.tsv.gz');effects=[]
for (arm,seed),d in mo[mo.subtype.eq('ActCD4')&mo.FOXP3_det.eq(0)].groupby(['filter','seed'],observed=True):
 z=d.groupby('sample',observed=True).agg(group=('group','first'),cells=('cell','size'),Ikzf2=('Ikzf2_auc','mean'),Stat1=('Stat1_auc','mean')).reset_index();z=z[z.cells>=20]
 for tf in ['Ikzf2','Stat1']:
  y=z[tf].to_numpy();X=np.column_stack([np.ones(len(z)),z.group.eq('Silica').astype(float)]);df=len(z)-2;b=np.linalg.lstsq(X,y,rcond=None)[0];se=np.sqrt(np.sum((y-X@b)**2)/df*np.linalg.inv(X.T@X)[1,1])
  effects.append(dict(filter=arm,seed=seed,tf=tf,n=len(z),coef=b[1],p=2*stats.t.sf(abs(b[1]/se),df),df=df,description='FOXP3_not_detected_description_only'))
wr(effects,'mouse_FOXP3_nondetect_sensitivity.tsv')
(OUT/'input_manifest.json').write_text(json.dumps(manifest,indent=2))
(OUT/'methods.txt').write_text('Marker-reviewed working annotation: MacAM clusters0/2/11/13, MonoMac3/7, MacWide union. Other clusters retained in annotation but excluded from macrophage sender models. AM clusters2/11 are subject-dominated and not new validated cell types. Sender aggregation uses cell-number-weighted original meanCP10k, verified against cell metadata. Matched >=20 cells at both ends. NoParent16 retains the previous exclusion sensitivity. NoDom omits CTL068/CTL070 only if present in paired subjects. Human receiver All/BG1, all five pre-existing CD4 scopes, three focused TFs, four fixed ligands, three seeds; never pool seeds as subjects. Spearman; SourceDep partial rank adjusts group/source/receiver library; BothEpi additionally sender library and epithelial fraction. No correlation P for n<10. BH within receiver filter/seed/sender filter/sample set/model across all sender/ligand/TF/receiver scope pairs. LOO same model; singular fits NA. Annotation models are exploratory, not cell-level doublet correction or causality. Mouse FOXP3 nondetection sensitivity is a descriptive subset, not a Treg-negative identity definition; nominal group model P, no confirmatory multiplicity claim.\n')
(OUT/'status.json').write_text(json.dumps(dict(status='PASS',stage='local_scenic3_reviewed_sender_models',tests=len(r),final_human_labels='marker_reviewed_working_annotations',causal_mechanism_proven=False),indent=2))
print('SCENIC3_REVIEW_PASS',len(r),'models')
