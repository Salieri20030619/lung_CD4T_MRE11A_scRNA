import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',MPLBACKEND='Agg')
from pathlib import Path
import sys,json,itertools,hashlib,re
import numpy as np,pandas as pd
from scipy.io import mmread
from scipy import sparse
from scipy.stats import false_discovery_control,permutation_test
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=Path(sys.argv[1]);M=P/'02_analysis/mechanism_step1_20260927_152218_s4V5eF';T=P/'02_analysis/GSE314943/stage9_20260906_080538_wN3sjO/tables';D=P/'02_analysis/GSE314943/stage6_20260905_231807_L44XdT/tables'
manifest=[]
def record(p,role):
 p=Path(p);manifest.append(dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),role=role));return p
def read(p,role,**kw):return pd.read_csv(record(p,role),sep='\t',**kw)
def wr(x,folder,name):x.to_csv(O/folder/name,sep='\t',index=False,encoding='utf-8',na_rep='NA')
def bh(v):
 a=np.asarray(v,float);q=np.full(a.shape,np.nan);ok=np.isfinite(a)
 if ok.any():q[ok]=false_discovery_control(a[ok],method='bh')
 return q
def canon(x):return {'MRE11A':'MRE11','ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','PARK2':'PRKN'}.get(str(x).upper(),str(x).upper())
mods=['PyroExec','Inflam','ApoExec','OXPHOS','MitoISR','MitoROS','Mitophagy','Autophagy_GO']
cells=read(M/'comm_input/CD4_cells.tsv','canonical counted CD4 cells');genes=read(M/'comm_input/CD4_genes.tsv','canonical mouse symbols')
x=mmread(record(M/'comm_input/CD4_counts.mtx.gz','canonical raw counts; gene by cell')).T.tocsr();assert x.shape==(len(cells),len(genes))
assert (np.asarray(x.sum(axis=1)).ravel()<=cells.library.values+1e-6).all()
y=x.multiply((1e4/cells.library.values)[:,None]).tocsr();y.data=np.log1p(y.data)
cellscore=read(M/'objects/GSE314943_cell_scores.tsv.gz','frozen matched-expression scores').set_index('cell').loc[cells.cell].reset_index()
assert (cellscore['sample'].values==cells['sample'].values).all()
auto=read(O/'inputs/GO0006914_autophagy_annotations.tsv','installed GO autophagy annotation');autogenes=set(auto.SYMBOL.map(canon));gidx=pd.Index(genes.canonical)
ix=np.array([j for j,g in enumerate(gidx) if g in autogenes]);assert len(ix)>=3
oldset=read(M/'tables/module_definition.tsv','frozen project panels');exclude=set(oldset.gene)|autogenes|set('IL1R2 PENK TNFRSF9 LAG3 TTC39C FOXP1 MRE11 IKZF2 STAT1'.split())
ref=np.mean([np.asarray(y[cells['sample'].values==s].mean(axis=0)).ravel() for s in sorted(cells['sample'].unique())],axis=0)
pool=np.array([g not in exclude and not re.match(r'^(MT-|RPL|RPS|MRPL|MRPS|IGH|IGK|IGL)',g) for g in gidx])&(ref>0)
# Background identity exclusions match the fixed old project panels, never a target-dependent threshold.
for g in 'SFTPA1 SFTPA2 SFTPB SFTPC SFTPD SCGB1A1 SCGB3A1 SCGB3A2 ABCA3 SLC34A2 NAPSA LYZ LYZ1 LYZ2 S100A8 S100A9 CSF1R LST1 TYROBP FCER1G'.split():
 if g in gidx:pool[gidx.get_loc(g)]=False
bins=pd.qcut(pd.Series(ref).rank(method='first'),25,labels=False).to_numpy();controls=set()
for j in ix:
 c=np.flatnonzero(pool&(bins==bins[j]));order=np.lexsort((c,np.abs(ref[c]-ref[j])));controls.update(c[order[:20]])
controls=np.array(sorted(controls));assert len(controls)>=20
score=np.asarray(y[:,ix].mean(axis=1)).ravel()-np.asarray(y[:,controls].mean(axis=1)).ravel();cellscore['score_Autophagy_GO']=score
cellscore['detect_Autophagy_GO']=np.asarray((x[:,ix]>0).sum(axis=1)).ravel()
wr(pd.DataFrame({'gene':gidx[ix],'role':'GO0006914_member'}),'02_mito_autophagy','autophagy_measured_members.tsv')
wr(pd.DataFrame({'gene':gidx[controls],'role':'matched_control'}),'02_mito_autophagy','autophagy_control_members.tsv')
wr(pd.DataFrame([dict(module='Autophagy_GO',annotated_symbols=len(autogenes),measured_symbols=len(ix),expressed_symbols=int((np.asarray(x[:,ix].sum(axis=0)).ravel()>0).sum()),controls=len(controls),coverage=len(ix)/len(autogenes),source='installed org.Mm.eg.db GOALL GO:0006914')]),'02_mito_autophagy','autophagy_coverage.tsv')
wr(cellscore[['cell','sample','subtype','group','library','bg']+['score_'+m for m in mods]],'02_mito_autophagy','cell_program_inputs.tsv.gz')
# Aggregate every nonempty library/state, equally weighted across libraries.
samples=[]
for filt in ['All','BG1']:
 z=cellscore if filt=='All' else cellscore[cellscore.bg<=.01]
 for scope in ['Core','Naive','ActCD4','Treg']:
  a=z if scope=='Core' else z[z.subtype==scope]
  for sample,b in a.groupby('sample'):
   r=dict(filter=filt,scope=scope,sample=sample,group=b.group.iloc[0],cells=len(b))
   for mod in mods:r[mod]=b['score_'+mod].mean()
   samples.append(r)
samples=pd.DataFrame(samples);wr(samples,'01_CD4_death','all_nonempty_sample_program_scores.tsv')
def compare(c,s):
 c=np.asarray(c,float);s=np.asarray(s,float);delta=s.mean()-c.mean();p=np.nan
 if min(len(c),len(s))>=2 and np.ptp(np.r_[c,s])>0:
  p=permutation_test((c,s),lambda a,b:abs(b.mean()-a.mean()),permutation_type='independent',n_resamples=np.inf,alternative='greater',vectorized=False).pvalue
 loo=[s.mean()-np.delete(c,j).mean() for j in range(len(c)) if len(c)>1]+[np.delete(s,j).mean()-c.mean() for j in range(len(s)) if len(s)>1]
 return dict(n_control=len(c),n_case=len(s),control_mean=c.mean(),case_mean=s.mean(),difference=delta,P=p,loo_same_sign=np.mean(np.sign(loo)==np.sign(delta)) if loo else np.nan)
results=[]
for (filt,scope),a in samples.groupby(['filter','scope']):
 for mod in mods:
  c=a[a.group=='Saline'][mod].dropna();s=a[a.group=='Silica'][mod].dropna();r=compare(c,s)
  r.update(filter=filt,scope=scope,module=mod,family=f'{filt}_'+('Core' if scope=='Core' else 'Subtypes')+'_'+('Death' if mod in mods[:3] else 'MitoAuto'),min_cells_observed=a.cells.min());results.append(r)
res=pd.DataFrame(results);res['BH']=res.groupby('family').P.transform(bh)
wr(res[res.module.isin(mods[:3])],'01_CD4_death','disease_program_comparisons.tsv');wr(res[res.module.isin(mods[3:])],'02_mito_autophagy','disease_program_comparisons.tsv')
# Direct within-library state contrasts; descriptive absolute level ranking is separate.
contrasts=[]
for filt in ['All','BG1']:
 for s1,s2 in itertools.combinations(['Naive','ActCD4','Treg'],2):
  a=samples[(samples['filter']==filt)&(samples.scope==s1)].set_index('sample');b=samples[(samples['filter']==filt)&(samples.scope==s2)].set_index('sample');common=a.index.intersection(b.index)
  for mod in mods:
   vals=a.loc[common,mod]-b.loc[common,mod];gr=a.loc[common,'group'];r=compare(vals[gr=='Saline'],vals[gr=='Silica']);r.update(filter=filt,contrast=f'{s1}-{s2}',module=mod,family=f'{filt}_state_response_'+('Death' if mod in mods[:3] else 'MitoAuto'));contrasts.append(r)
contrasts=pd.DataFrame(contrasts);contrasts['BH']=contrasts.groupby('family').P.transform(bh);wr(contrasts,'01_CD4_death','between_state_differential_response.tsv')
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
for panel,name in [(mods[:3],'01_death_programs'),(mods[3:],'02_mito_autophagy_programs')]:
 fig,axs=plt.subplots(4,len(panel),figsize=(3.4*len(panel),10.5),squeeze=False)
 for i,scope in enumerate(['Core','Naive','ActCD4','Treg']):
  for j,mod in enumerate(panel):
   ax=axs[i,j];a=samples[(samples['filter']=='All')&(samples.scope==scope)]
   for k,group in enumerate(['Saline','Silica']):
    b=a[a.group==group].sort_values('sample');ax.scatter(k+np.linspace(-.08,.08,len(b)),b[mod],c=['#2c6f9f','#c84d3e'][k],s=30);ax.hlines(b[mod].mean(),k-.15,k+.15,colors='black')
   r=res[(res['filter']=='All')&(res.scope==scope)&(res.module==mod)].iloc[0]
   ax.set(title=f'{scope}: {mod}\ndelta {r.difference:+.4f}; P {r.P:.2g}; BH {r.BH:.2g}',xticks=[0,1],xticklabels=['Saline','Silica'],ylabel='Sample mean program score' if j==0 else '')
 fig.suptitle('All six nonempty libraries retained; no minimum cell-count exclusion',fontsize=15)
 fig.tight_layout(rect=[0,0,1,.96]);fig.savefig(O/'figures'/(name+'.png'),dpi=170,bbox_inches='tight');fig.savefig(O/'figures'/(name+'.svg'),bbox_inches='tight');plt.close(fig)
print('STAGES 1/2 COMPLETE',res.shape,flush=True)
# Pseudotime inferred upstream is reused unchanged; inference now uses independent-library slopes.
pt=read(T/'pseudotime_all_runs.tsv.gz','four fixed existing CD4 DPT orderings',index_col=0)
assert len(pt)==1734 and set(pt.index)==set(cells.loc[cells.subtype.isin(['Naive','ActCD4']),'cell'])
ptmeta=cells.set_index('cell').loc[pt.index].join(pt)
wr(ptmeta.reset_index(),'03A_state_pseudotime','pseudotime_cell_membership.tsv.gz')
wr(read(T/'trajectory_stability.tsv','existing root/graph stability'),'03A_state_pseudotime','existing_ordering_stability.tsv')
sloperows=[];generes=[];programres=[];qc=[];binrows=[]
for state in ['Naive','ActCD4']:
 de=read(D/f'{state}_Full_DEG.tsv','frozen all-library discovery DE');de['canonical']=de.symbol.map(canon);assert not de.canonical.duplicated().any()
 selected=de[de.canonical.isin(gidx)].copy();idx=gidx.get_indexer(selected.canonical);names=selected.canonical.to_list()
 wr(de,'04_candidate_intersection',f'{state}_discovery_DE_frozen.tsv')
 base=(cells.subtype==state).to_numpy();global_rows=np.flatnonzero(base);ss=cells.loc[base,'sample'].to_numpy();cc=cells.loc[base,'cell'].to_numpy()
 yy=y[base][:,idx].toarray();pp=cellscore.loc[base,['score_'+m for m in mods]].to_numpy();lib=np.log1p(cells.loc[base,'library'].to_numpy());allvalues=np.column_stack([yy,pp]);features=names+mods
 for run in pt.columns:
  tx=pt.loc[cc,run].to_numpy();tx=(tx-tx.min())/(tx.max()-tx.min());slopes=[]
  for sample in sorted(set(ss)):
   take=ss==sample;z=np.column_stack([np.ones(take.sum()),lib[take]-lib[take].mean()]);u=tx[take];u=u-z@np.linalg.lstsq(z,u,rcond=None)[0]
   if take.sum()<=np.linalg.matrix_rank(z) or np.dot(u,u)<1e-10:beta=np.full(len(features),np.nan)
   else:beta=u@allvalues[take]/np.dot(u,u)
   slopes.append(beta)
   qc.append(dict(state=state,run=run,sample=sample,cells=int(take.sum()),pt_min=tx[take].min(),pt_max=tx[take].max(),depth_adjusted_pt_ss=float(u@u),slope_estimable=bool(np.isfinite(beta).all())))
   for f,v in zip(features,beta):sloperows.append(dict(state=state,run=run,sample=sample,feature=f,type='gene' if f in names else 'program',slope=v))
   if run=='K15_R1':
    # Within-state sample quantiles for visualization only; bins are never replicates.
    bins=pd.qcut(pd.Series(tx[take]).rank(method='first'),4,labels=False).to_numpy()
    for q in range(4):
     for m,score in zip(mods,pp[take].T):binrows.append(dict(state=state,sample=sample,bin=q,cells=int((bins==q).sum()),mean_pt=tx[take][bins==q].mean(),module=m,score=score[bins==q].mean()))
  slopes=np.asarray(slopes);n=slopes.shape[0];assert n==6 and np.isfinite(slopes).all()
  means=slopes.mean(axis=0);signs=np.array(list(itertools.product([-1,1],repeat=n)));null=signs@slopes/n
  pv=(np.abs(null)>=np.abs(means)[None,:]-1e-12).mean(axis=0);pv[np.max(np.abs(slopes),axis=0)<1e-12]=np.nan
  for j,f in enumerate(features):
   row=dict(state=state,run=run,feature=f,mean_slope=means[j],P=pv[j],n_libraries=n,positive_libraries=int((slopes[:,j]>0).sum()),negative_libraries=int((slopes[:,j]<0).sum()),same_sign_fraction=float((np.sign(slopes[:,j])==np.sign(means[j])).mean()),min_slope=slopes[:,j].min(),max_slope=slopes[:,j].max(),method='exact 64 sample-level sign flips; depth-adjusted expression slope')
   (generes if j<len(names) else programres).append(row)
generes=pd.DataFrame(generes);generes['family']=generes.state+'_'+generes.run+'_all_tested_genes';generes['BH']=generes.groupby('family').P.transform(bh)
programres=pd.DataFrame(programres);programres['family']=programres.run+'_all_state_programs';programres['BH']=programres.groupby('family').P.transform(bh)
wr(generes,'03A_state_pseudotime','genomewide_trajectory_tests.tsv.gz');wr(pd.DataFrame(sloperows),'03A_state_pseudotime','sample_level_slopes.tsv.gz');wr(programres,'03A_state_pseudotime','program_trajectory_tests.tsv');wr(pd.DataFrame(qc),'03A_state_pseudotime','sample_trajectory_coverage.tsv');wr(pd.DataFrame(binrows),'03A_state_pseudotime','program_sample_quantile_curves.tsv')
fig,axs=plt.subplots(2,4,figsize=(16,7))
binrows=pd.DataFrame(binrows)
for i,state in enumerate(['Naive','ActCD4']):
 for j,mod in enumerate(['PyroExec','ApoExec','Mitophagy','Autophagy_GO']):
  ax=axs[i,j];a=binrows[(binrows.state==state)&(binrows.module==mod)]
  for sample,b in a.groupby('sample'):ax.plot(b.mean_pt,b.score,'o-',lw=.8,alpha=.6,c='#2c6f9f' if sample.startswith('Saline') else '#c84d3e',ms=3)
  r=programres[(programres.state==state)&(programres.run=='K15_R1')&(programres.feature==mod)].iloc[0]
  ax.set(title=f'{state}: {mod}\nslope {r.mean_slope:+.3f}; P {r.P:.3g}; BH {r.BH:.3g}',xlabel='Within-state rescaled DPT',ylabel='Program score' if j==0 else '')
fig.suptitle('Within-state trajectories: six sample curves; P uses six library slopes',fontsize=15);fig.tight_layout(rect=[0,0,1,.95]);fig.savefig(O/'figures/03A_program_trajectories.png',dpi=170,bbox_inches='tight');plt.close(fig)
wr(pd.DataFrame(manifest).drop_duplicates('path'),'provenance','input_manifest_phase123A.tsv')
status={'status':'PASS','cell_count_exclusion':'none','samples_per_main_scope':samples[samples['filter']=='All'].groupby('scope').size().to_dict(),'autophagy_genes':len(ix),'gene_trajectory_tests':len(generes),'primary_nominal_gene_counts':generes[generes.run=='K15_R1'].groupby('state').P.apply(lambda v:int((v<.05).sum())).to_dict(),'program_trajectory_tests':len(programres),'state_cells':cells.groupby('subtype').size().to_dict()}
(O/'provenance/phase123A_complete.json').write_text(json.dumps(status,indent=2),encoding='utf-8');print(json.dumps(status),flush=True)