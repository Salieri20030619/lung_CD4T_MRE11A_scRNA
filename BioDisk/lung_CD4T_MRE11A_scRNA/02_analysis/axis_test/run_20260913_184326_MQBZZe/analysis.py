import os, sys, json, itertools, platform
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R=Path(os.environ['PROJECT_ROOT']); O=Path(sys.argv[1]); A=R/'02_analysis'
for sub in ['tables','figures']: (O/sub).mkdir()
genes=['IL18','CXCL16','IL18R1','IL18RAP','CXCR6']
specs=[
 ('GSE314943','Mac','Core',A/'GSE314943/stage11_20260906_100314_pALwAB/objects/macrophages_core.h5ad',1796),
 ('GSE314943','CD4','Core',A/'GSE314943/stage5_20260905_222457_a2dApT/objects/CD4_core.h5ad',2311),
 ('GSE183682','CD4','Core',A/'screen2_CD4_target_wz1lfx/GSE183682/CD4_core.h5ad',1413),
 ('GSE174725','CD4','Core',A/'screen2_Trefine_nuzVQu/GSE174725/T_refined.h5ad',10286),
 ('GSE274372','TLike','TLike',A/'GSE274372/stage9_20260908_140721_Yj56oq/objects/CD4_author_reclustered.h5ad',888)]
assert all(p.is_file() for _,_,_,p,_ in specs), 'Missing fixed input: '+str([str(p) for _,_,_,p,_ in specs if not p.is_file()])
def wr(df,name): df.to_csv(O/'tables'/name,sep='\t',index=False,na_rep='NA')
rows=[]; cellframes=[]; audits=[]; mouse_ids={}
for dataset,lineage,scope,path,expected in specs:
 print('DATASET_START',dataset,lineage,flush=True)
 a=ad.read_h5ad(path)
 assert a.n_obs==expected and a.obs_names.is_unique
 assert 'counts' in a.layers and 'sample' in a.obs
 x=sparse.csr_matrix(a.layers['counts']);x.eliminate_zeros()
 assert np.isfinite(x.data).all() and (x.data>=0).all()
 assert np.allclose(x.data,np.rint(x.data),atol=1e-6,rtol=0)
 lib=np.asarray(x.sum(axis=1)).ravel();assert (lib>0).all()
 ids=pd.Index(a.var_names.astype(str)).str.replace(r'\.\d+$','',regex=True)
 if 'symbol' in a.var:
  symbols=a.var['symbol'].astype('string').fillna('').str.strip().str.upper().to_numpy(dtype=str)
 else: symbols=ids.str.upper().to_numpy()
 if dataset=='GSE314943':
  for gene in genes:
   jj=np.flatnonzero(symbols==gene);assert len(jj)>0, gene
   mouse_ids[gene]=sorted(set(ids[jj]))
 if dataset=='GSE274372':
  for gene in genes:
   jj=np.flatnonzero(ids.isin(mouse_ids[gene]));assert len(jj)>0, 'Missing mapped gene '+gene
   symbols[jj]=gene
 vals={}
 for gene in genes:
  jj=np.flatnonzero(symbols==gene);assert len(jj)>0, dataset+' missing '+gene
  vals[gene]=np.asarray(x[:,jj].sum(axis=1)).ravel()
  audits.append(dict(dataset=dataset,lineage=lineage,path=str(path),gene=gene,features=';'.join(ids[jj]),n_features=len(jj),full_features=x.shape[1],input_cells=a.n_obs))
 m=pd.DataFrame({'cell':a.obs_names.astype(str),'sample':a.obs['sample'].astype(str).to_numpy(),'library_UMI':lib,**vals})
 m['day']=-1
 if dataset=='GSE314943':
  m['subtype']=a.obs['MacType' if lineage=='Mac' else 'Subtype'].astype(str).to_numpy()
  m['group']=a.obs['group'].astype(str).to_numpy()
  assert set(m.group)=={'Saline','Silica'} and m['sample'].nunique()==6
 elif dataset=='GSE183682':
  mp={'GSM5567501':('NS',7),'GSM5567502':('Silica',7),'GSM5567503':('NS',56),'GSM5567504':('Silica',56)}
  m['group']=m['sample'].map({k:v[0] for k,v in mp.items()});m['day']=m['sample'].map({k:v[1] for k,v in mp.items()})
  m['subtype']=a.obs.F08.astype(str).map({'0':'NaiveCD4','3':'ActCD4','5':'Treg'}).to_numpy()
 elif dataset=='GSE174725':
  m['subtype']=a.obs.F08.astype(str).map({'0':'MemCD4','1':'NaiveCD4','12':'Treg','2':'Mix2','3':'Mix3'}).to_numpy()
  m=m[m.subtype.notna()].copy()
  m['group']=m['sample'].map({'GSM5324866':'Exposed','GSM5324867':'Exposed','GSM5324868':'Silicosis','GSM5324869':'Silicosis','GSM5324870':'Silicosis'})
  assert len(m)==5481 and (~m.subtype.isin(['Mix2','Mix3'])).sum()==2958
 else:
  m['subtype']=a.obs.R08.astype(str).map({'0':'ActReg','1':'TLike1'}).to_numpy()
  m['day']=pd.to_numeric(a.obs.day).to_numpy();m['group']='D'+m.day.astype(int).astype(str)
  m=m[m.subtype.notna()].copy();assert len(m)==494 and m['sample'].nunique()==12
 assert m[['sample','subtype','group','day']].notna().all().all()
 m['dataset']=dataset;m['lineage']=lineage
 m['IL18pair']=(m.IL18R1>0)&(m.IL18RAP>0)
 cellframes.append(m)
 scopes={'Core':m} if dataset!='GSE274372' else {'TLike':m}
 if dataset=='GSE174725': scopes={'Core':m[~m.subtype.isin(['Mix2','Mix3'])],'Expanded':m}
 for scp,z in scopes.items():
  for subtype in ['All']+sorted(z.subtype.unique()):
   for sample in sorted(z['sample'].unique()):
    base=z[z['sample']==sample];v=base if subtype=='All' else base[base.subtype==subtype]
    n=len(v)
    for gene in genes+['IL18pair']:
     val=v[gene].to_numpy();den=v.library_UMI.to_numpy()
     rows.append(dict(dataset=dataset,lineage=lineage,scope=scp,subtype=subtype,sample=sample,group=base.group.iloc[0],day=base.day.iloc[0],gene=gene,cells=n,positive=int((val>0).sum()),pct=100*np.mean(val>0) if n else np.nan,CPM=1e6*val.sum()/den.sum() if n and gene!='IL18pair' else np.nan,CP10k=np.mean(val/den*1e4) if n and gene!='IL18pair' else np.nan,mean_UMI=den.mean() if n else np.nan,eligible20=n>=20))
 del a,x
 print('DATASET_DONE',dataset,lineage,flush=True)
tab=pd.DataFrame(rows);cells=pd.concat(cellframes,ignore_index=True)
assert not tab.duplicated(['dataset','lineage','scope','subtype','sample','gene']).any()
assert ((tab.pct.dropna()>=0)&(tab.pct.dropna()<=100)).all()
print('ANALYSIS_LOGIC_PASS RAW_COUNTS_FIXED_SCOPE_SAMPLE_UNIT',flush=True)
wr(tab,'sample_expression.tsv');wr(pd.DataFrame(audits),'input_gene_audit.tsv')
cells.to_csv(O/'tables/target_cells.tsv.gz',sep='\t',index=False)

# Group means weight independent samples equally. No cell-level p-values.
results=[]
for keys,v in tab.groupby(['dataset','lineage','scope','subtype','gene'],observed=True):
 dataset,lineage,scp,subtype,gene=keys
 pairs=[('Saline','Silica',-1)] if dataset=='GSE314943' else [('NS','Silica',7),('NS','Silica',56)] if dataset=='GSE183682' else [('Exposed','Silicosis',-1)] if dataset=='GSE174725' else [('D0','D7',-1),('D0','D28',-1),('D0','D56',-1)]
 for g0,g1,day in pairs:
  vv=v if day==-1 else v[v.day==day]
  for threshold in [1,20]:
   z=vv[vv.cells>=threshold]
   for metric in ['pct','CPM','CP10k']:
    b=z.loc[z.group==g0,metric].dropna().to_numpy();c=z.loc[z.group==g1,metric].dropna().to_numpy()
    effect=c.mean()-b.mean() if len(b) and len(c) else np.nan
    p=np.nan
    if dataset=='GSE314943' and len(b)>=2 and len(c)>=2:
     y=np.r_[b,c];null=[]
     for ix in itertools.combinations(range(len(y)),len(c)):
      take=np.zeros(len(y),bool);take[list(ix)]=True
      null.append(y[take].mean()-y[~take].mean())
     p=np.mean(np.abs(null)>=abs(effect)-1e-12)
    results.append(dict(dataset=dataset,lineage=lineage,scope=scp,subtype=subtype,gene=gene,day=day,control=g0,case=g1,min_cells=threshold,metric=metric,n_control=len(b),n_case=len(c),control_mean=b.mean() if len(b) else np.nan,case_mean=c.mean() if len(c) else np.nan,difference=effect,p=p,method='exact_sample_permutation' if np.isfinite(p) else 'descriptive',pairwise_higher_fraction=np.mean(c[:,None]>b[None,:]) if len(b) and len(c) else np.nan))
wr(pd.DataFrame(results),'group_directions.tsv')

# Co-detection is an expression prerequisite; no communication probability constructed.
co=[];q=tab[(tab.dataset=='GSE314943')&(tab.scope=='Core')]
for sample in sorted(q['sample'].unique()):
 z=q[q['sample']==sample]
 for axis,lig,rec in [('IL18','IL18','IL18pair'),('CXCL16','CXCL16','CXCR6')]:
  for mac in ['All','AMlike','C1qMac','Folr2Mac','Spp1Mac']:
   for tt in ['All','Naive','ActCD4','Treg']:
    b=z[(z.lineage=='Mac')&(z.subtype==mac)&(z.gene==lig)].iloc[0]
    c=z[(z.lineage=='CD4')&(z.subtype==tt)&(z.gene==rec)].iloc[0]
    co.append(dict(sample=sample,group=b.group,axis=axis,mac=mac,CD4=tt,mac_cells=b.cells,CD4_cells=c.cells,ligand_positive=b.positive,receptor_positive=c.positive,ligand_pct=b.pct,receptor_pct=c.pct,both_detected=bool(b.positive>0 and c.positive>0),both_eligible20=bool(b.cells>=20 and c.cells>=20)))
wr(pd.DataFrame(co),'sample_codetection.tsv')

# Fixed subtype weights separate composition from expression direction in main cohort.
standard=[]
for (lineage,gene),v in q[q.subtype!='All'].groupby(['lineage','gene'],observed=True):
 ct=v.pivot(index='sample',columns='subtype',values='cells');w=ct.div(ct.sum(axis=1),axis=0).mean(axis=0)
 for metric in ['pct','CP10k']:
  if gene=='IL18pair' and metric=='CP10k': continue
  mat=v.pivot(index='sample',columns='subtype',values=metric).reindex(columns=w.index)
  for sample,row in mat.iterrows():
   complete=row.notna().all()
   standard.append(dict(lineage=lineage,gene=gene,sample=sample,group=v.loc[v['sample']==sample,'group'].iloc[0],metric=metric,value=float(row@w) if complete else np.nan,status='complete' if complete else 'missing_subtype',weights=json.dumps(w.to_dict())))
wr(pd.DataFrame(standard),'main_standardized.tsv')
for dataset in tab.dataset.unique():
 z=tab[(tab.dataset==dataset)&(tab.scope!='Expanded')&(tab.subtype!='All')]
 panels=[('Mac',g) for g in ['IL18','CXCL16']]+[('CD4',g) for g in ['IL18R1','IL18RAP','IL18pair','CXCR6']]
 if dataset=='GSE274372':panels=[('TLike',g) for g in ['IL18R1','IL18RAP','IL18pair','CXCR6']]
 panels=[p for p in panels if len(z[(z.lineage==p[0])&(z.gene==p[1])])]
 fig,axs=plt.subplots(2,3,figsize=(12,7),squeeze=False)
 for ax,(lin,gene) in zip(axs.flat,panels):
  v=z[(z.lineage==lin)&(z.gene==gene)];types=sorted(v.subtype.unique());groups=sorted(v.group.unique())
  for j,gr in enumerate(groups):
   for i,ty in enumerate(types):
    rr=v[(v.group==gr)&(v.subtype==ty)&(v.cells>0)].sort_values('sample')
    xx=i+(j-(len(groups)-1)/2)*.16+np.linspace(-.04,.04,len(rr))
    color=plt.get_cmap('tab10')(j)
    for xp,(_,r) in zip(xx,rr.iterrows()):ax.scatter(xp,r.pct,s=30,edgecolor=color,facecolor=color if r.cells>=20 else 'none')
   ax.scatter([],[],color=plt.get_cmap('tab10')(j),label=gr)
  ax.set(xticks=range(len(types)),xticklabels=types,ylabel='Detect (%)',title=lin+' / '+gene);ax.tick_params(axis='x',rotation=25);ax.legend(frameon=False,fontsize=7);ax.spines[['top','right']].set_visible(False)
 for ax in list(axs.flat)[len(panels):]:ax.set_visible(False)
 fig.tight_layout()
 for ext in ['png','pdf']:fig.savefig(O/'figures'/f'{dataset}_detect.{ext}',dpi=350,bbox_inches='tight')
 plt.close(fig)

# Inventory supplementary full-lung objects for the next source-specific check.
inventory=[]
for dataset in ['GSE183682','GSE274372','GSE174725']:
 for p in A.rglob('*.h5ad'):
  if dataset in str(p) and p not in [s[3] for s in specs]:
   inventory.append(dict(dataset=dataset,path=str(p),size_bytes=p.stat().st_size))
wr(pd.DataFrame(inventory,columns=['dataset','path','size_bytes']),'supplement_object_inventory.tsv')
(O/'captions.txt').write_text('Detection panels: x=pre-existing cell subtype, y=percentage with positive raw counts. Each point is an independent sample/library. Filled points >=20 cells; hollow points 1-19 cells. IL18pair=both IL18R1 and IL18RAP detected in the same cell; double dropout affects this measurement. AMlike=alveolar-like macrophage; C1qMac=C1q-positive macrophage; Folr2Mac=Folr2-positive macrophage; Spp1Mac=Spp1-positive macrophage; Naive/NaiveCD4=naive-like CD4; ActCD4=activated-like CD4; Treg=regulatory T; MemCD4=memory-like CD4; ActReg=activated/regulatory-like T candidate; TLike1=R08 cluster1 T-like candidate. GSE274372 is a lower-depth single-nucleus T-like scope, not equivalent to high-confidence CD4. NS/Saline=saline; Silica=silica-exposed; Exposed=human exposed comparator; Silicosis=human disease; D0/D7/D28/D56=sampling day. Plots show detection, not protein activation or communication significance.\n')
(O/'methods.txt').write_text('Fixed five target genes, original raw counts, complete measured count denominator. Duplicate symbols summed. Ensembl mapping for GSE274372 taken from actual mouse reference object and recorded. Absent subtypes have n=0 and missing expression. Main sample-level exact two-sided permutations are exploratory; all tested outcomes retained. Other cohorts descriptive: GSE183682 one library per time/condition; GSE174725 2 exposed/3 disease with Core and Expanded separate; GSE274372 day-associated direction without later contemporaneous controls. Threshold1 descriptive, threshold20 sensitivity. CPM is pseudobulk; CP10k is mean per-cell normalized counts. Fixed equal-sample subtype weights used only where all subtypes present. No new CellChat fit, no cross-model probability comparison. Main macrophage and CD4 endpoints tested; supplementary CD4/T-like receptor endpoints tested here. Supplementary macrophage-source validation and receptor-associated functional programs remain pending.\n')
(O/'status.json').write_text(json.dumps({'status':'PASS','completed':'main_Mac_CD4_expression_and_supplement_receptor_screen','pending':['supplement_macrophage_source_validation','receptor_associated_functional_programs'],'python':platform.python_version(),'anndata':ad.__version__,'inputs':[str(s[3]) for s in specs]},indent=2))
print('AXIS_TEST_ANALYSIS_PASS',flush=True)
