import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',MPLBACKEND='Agg')
from pathlib import Path
import sys,json,hashlib,itertools
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(sys.argv[1]);P=O.parents[1];V=P/'11_reports/customer_sequence_20261015_v2'
G=['AHNAK','MDFIC','TTC39C','ST6GALNAC3']
PAN={'NaiveIdentity':'TCF7 LEF1 CCR7 SELL IL7R','TregIdentity':'FOXP3 IL2RA CTLA4 IKZF2 IKZF4 TIGIT','Activation':'CD69 ICOS NR4A1 BATF IRF4 TNFRSF4','Cytotoxic':'NKG7 PRF1 GZMB CCL5 CTSW CD8A CD8B1'}
PAN={k:v.split() for k,v in PAN.items()};BG='SFTPA1 SFTPA2 SFTPB SFTPC SFTPD SCGB1A1 SCGB3A1 SCGB3A2 ABCA3 SLC34A2 NAPSA LYZ LYZ1 LYZ2 S100A8 S100A9 CSF1R LST1 TYROBP FCER1G'.split()
manifest=[];cellframes=[];coverage=[];mapping=[]
def rec(p):
 p=Path(p);manifest.append(dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()));return p
def wr(a,n):a.to_csv(O/'direction'/n,sep='\t',index=False,encoding='utf-8',na_rep='NA')
paths={'GSE314943':P/'02_analysis/GSE314943/stage5_20260905_222457_a2dApT/objects/CD4_core.h5ad','GSE183682':P/'02_analysis/screen2_CD4_target_wz1lfx/GSE183682/CD4_core.h5ad'}
for ds,path in paths.items():
 a=ad.read_h5ad(rec(path));x=sparse.csr_matrix(a.layers['counts'] if 'counts' in a.layers else a.X);assert np.allclose(x.data,np.rint(x.data));lib=np.asarray(x.sum(axis=1)).ravel();assert (lib>0).all()
 sym=a.var['symbol'].astype(str).str.upper().to_numpy() if 'symbol' in a.var else a.var_names.str.upper().to_numpy()
 md=pd.DataFrame({'cell':a.obs_names,'sample':a.obs['sample'].astype(str).to_numpy(),'library':lib,'dataset':ds})
 if ds=='GSE314943':md['scope']=a.obs.Subtype.astype(str).to_numpy();md['group']=a.obs.group.astype(str).to_numpy();md['day']=-1
 else:
  mp={'GSM5567501':('Saline',7),'GSM5567502':('Silica',7),'GSM5567503':('Saline',56),'GSM5567504':('Silica',56)}
  md['scope']=a.obs.F08.astype(str).map({'0':'Naive','3':'ActCD4','5':'Treg'}).to_numpy();md['group']=md['sample'].map({k:v[0] for k,v in mp.items()});md['day']=md['sample'].map({k:v[1] for k,v in mp.items()})
 assert md.scope.notna().all()
 for gene in G:
  ix=np.flatnonzero(sym==gene);assert len(ix)==1,(ds,gene,len(ix));v=np.asarray(x[:,ix].sum(axis=1)).ravel();md[gene+'_count']=v;md[gene+'_CP10k']=v/lib*1e4;md[gene+'_logCP10k']=np.log1p(v/lib*1e4)
  mapping.append(dict(dataset=ds,gene=gene,feature_id=';'.join(a.var_names[ix]),symbol_matches=len(ix),total_UMI=v.sum()))
 bg=np.isin(sym,BG)|np.array([s.startswith(('IGH','IGK','IGL')) for s in sym]);md['bg_fraction']=np.asarray(x[:,bg].sum(axis=1)).ravel()/lib
 for name,gs in PAN.items():
  ix=np.flatnonzero(np.isin(sym,gs));z=x[:,ix].multiply((1e4/lib)[:,None]).tocsr();z.data=np.log1p(z.data);md[name]=np.asarray(z.mean(axis=1)).ravel();coverage.append(dict(dataset=ds,panel=name,members=';'.join(sym[ix]),measured=len(ix),expected=len(gs)))
 cellframes.append(md);del a,x
cells=pd.concat(cellframes,ignore_index=True);wr(pd.DataFrame(mapping),'candidate_raw_feature_mapping.tsv');wr(pd.DataFrame(coverage),'identity_panel_coverage.tsv');wr(cells,'mouse_candidate_cell_inputs.tsv.gz')
metrics=[];identity=[]
for filt in ['All','BG1']:
 z=cells if filt=='All' else cells[cells.bg_fraction<=.01]
 for (ds,day,sample,group),a in z.groupby(['dataset','day','sample','group']):
  for scope in ['Core','Naive','ActCD4','Treg']:
   b=a if scope=='Core' else a[a.scope==scope]
   if not len(b):continue
   identity.append(dict(dataset=ds,day=day,sample=sample,group=group,filter=filt,scope=scope,cells=len(b),median_library=b.library.median(),mean_bg=b.bg_fraction.mean(),**{k:b[k].mean() for k in PAN}))
   for gene in G:
    c=b[gene+'_count'].to_numpy();v=b[gene+'_CP10k'].to_numpy();pos=c>0
    metrics.append(dict(dataset=ds,day=day,sample=sample,group=group,filter=filt,scope=scope,gene=gene,cells=len(b),positive_cells=int(pos.sum()),target_UMI=float(c.sum()),full_library_UMI=float(b.library.sum()),detection_fraction=pos.mean(),mean_CP10k=v.mean(),pseudobulk_CPM=c.sum()/b.library.sum()*1e6,mean_logCP10k=np.log1p(v).mean(),positive_mean_CP10k=v[pos].mean() if pos.any() else np.nan,max_cell_expression_share=v.max()/v.sum() if v.sum()>0 else np.nan))
met=pd.DataFrame(metrics);wr(met,'sample_expression_estimands.tsv');wr(pd.DataFrame(identity),'sample_identity_context.tsv')
contrasts=[];decomp=[]
for (ds,day,filt,scope,gene),a in met.groupby(['dataset','day','filter','scope','gene']):
 c=a[a.group=='Saline'];s=a[a.group=='Silica']
 if not len(c) or not len(s):continue
 base=dict(dataset=ds,day=day,filter=filt,scope=scope,gene=gene,n_control=len(c),n_case=len(s),method='equal-library group means; diagnostic only')
 for metric in ['mean_CP10k','pseudobulk_CPM','mean_logCP10k','detection_fraction','positive_mean_CP10k']:
  cm=c[metric].mean();sm=s[metric].mean();contrasts.append(dict(**base,metric=metric,control=cm,case=sm,difference=sm-cm,log2_case_control=np.log2(sm/cm) if cm>0 and sm>0 else np.nan))
 pc=c.detection_fraction.mean();ps=s.detection_fraction.mean();mc=c.mean_CP10k.mean();ms=s.mean_CP10k.mean();uc=mc/pc if pc>0 else np.nan;us=ms/ps if ps>0 else np.nan
 det=(ps-pc)*(uc+us)/2;intensity=(us-uc)*(pc+ps)/2
 if np.isfinite(det+intensity):assert abs(det+intensity-(ms-mc))<1e-9
 decomp.append(dict(**base,delta=ms-mc,control_detection=pc,case_detection=ps,control_positive_intensity=uc,case_positive_intensity=us,detection_component=det,positive_intensity_component=intensity))
con=pd.DataFrame(contrasts);wr(con,'direction_by_metric_and_filter.tsv');wr(pd.DataFrame(decomp),'detection_intensity_decomposition.tsv')
# Core composition decomposition using equal-library mixtures, exact reconstruction.
comp=[]
for (ds,day,filt,gene),a in met.groupby(['dataset','day','filter','gene']):
 sides={}
 for group in ['Saline','Silica']:
  b=a[a.group==group];core=b[b.scope=='Core'].set_index('sample');sub=b[b.scope!='Core'].copy();sub['w']=sub.cells/sub['sample'].map(core.cells);sub['contribution']=sub.w*sub.mean_CP10k
  vals={}
  for state,v in sub.groupby('scope'):
   w=v.w.sum()/len(core);mu=v.contribution.sum()/len(core)/w if w>0 else np.nan;vals[state]=(w,mu)
  sides[group]=vals
 for state in ['Naive','ActCD4','Treg']:
  wc,mc=sides['Saline'][state];ws,ms=sides['Silica'][state];comp.append(dict(dataset=ds,day=day,filter=filt,gene=gene,state=state,control_fraction=wc,case_fraction=ws,control_conditional_expression=mc,case_conditional_expression=ms,composition_component=(ws-wc)*(mc+ms)/2,within_state_component=(ms-mc)*(wc+ws)/2))
comp=pd.DataFrame(comp);wr(comp,'Core_state_composition_decomposition.tsv')
for key,a in comp.groupby(['dataset','day','filter','gene']):
 b=con[(con.dataset==key[0])&(con.day==key[1])&(con['filter']==key[2])&(con.gene==key[3])&(con.scope=='Core')&(con.metric=='mean_CP10k')].iloc[0]
 assert abs(a.composition_component.sum()+a.within_state_component.sum()-b.difference)<1e-8
# Marker-stratified standardization in ActCD4; all prespecified diagnostics retained.
# Fixed tertiles use equal-library weights, are blind to candidate expression.
def wquant(v,w,q):
 i=np.argsort(v,kind='stable');return np.interp(q,np.cumsum(w[i])/w.sum(),v[i])
standard=[];strata=[]
for (ds,day),a in cells[cells.scope=='ActCD4'].groupby(['dataset','day']):
 a=a.copy();sample_n=a.groupby('sample').size();w=a['sample'].map(1/sample_n).to_numpy()
 for panel in PAN:
  cut=np.unique(wquant(a[panel].to_numpy(),w,np.array([1/3,2/3])));bins=np.searchsorted(cut,a[panel].to_numpy(),side='left');a['bin']=bins
  ct=pd.crosstab(a['sample'],a.bin);common=ct.columns[(ct>0).all(axis=0)].to_list()
  for sample,b in a.groupby('sample'):
   for bin in ct.columns:
    v=b[b.bin==bin];strata.append(dict(dataset=ds,day=day,panel=panel,sample=sample,group=b.group.iloc[0],stratum=bin,cells=len(v),total_cells=len(b),common_support=bin in common,cutpoints=';'.join(map(str,cut))))
  if not common:continue
  ctab=ct[common];prop=ctab.div(ctab.sum(axis=1),axis=0);groups=a.groupby('sample').group.first();weights=(prop.loc[groups=='Saline'].mean()+prop.loc[groups=='Silica'].mean())/2
  for gene in G:
   vals=[]
   for sample,b in a.groupby('sample'):
    means=b.groupby('bin')[gene+'_CP10k'].mean();z=sum(weights[bin]*means[bin] for bin in common);vals.append((sample,b.group.iloc[0],z,b[b.bin.isin(common)][gene+'_CP10k'].mean(),len(b[b.bin.isin(common)])/len(b)))
   v=pd.DataFrame(vals,columns=['sample','group','std','restricted_raw','coverage']);delta=v[v.group=='Silica']['std'].mean()-v[v.group=='Saline']['std'].mean();raw=v[v.group=='Silica'].restricted_raw.mean()-v[v.group=='Saline'].restricted_raw.mean()
   standard.append(dict(dataset=ds,day=day,scope='ActCD4',gene=gene,panel=panel,standardized_difference=delta,restricted_raw_difference=raw,min_cell_coverage=v.coverage.min(),common_bins=';'.join(map(str,common)),sign_changed=np.sign(delta)!=np.sign(raw),interpretation='descriptive marker-stratified standardization; not reannotation or causal adjustment'))
wr(pd.DataFrame(standard),'ActCD4_marker_standardization.tsv');wr(pd.DataFrame(strata),'ActCD4_marker_stratum_coverage.tsv')
# Human contradiction: keep all donors, explicit state and frozen sensitivity masks.
h=pd.read_csv(rec(V/'05_supporting_evidence/GSE174725_candidate_donor_expression.tsv'),sep='\t');hr=pd.read_csv(rec(V/'05_supporting_evidence/GSE174725_candidate_group_comparisons.tsv'),sep='\t');wr(hr,'human_scope_group_directions.tsv');loo=[]
for (filt,scope,gene),a in h.groupby(['filter','scope','gene']):
 a=a.dropna(subset=['CP10k']);base=a[a.group=='Silicosis'].CP10k.mean()-a[a.group=='Exposed'].CP10k.mean()
 for sample in a['sample']:
  b=a[a['sample']!=sample];c=b[b.group=='Exposed'];s=b[b.group=='Silicosis'];delta=s.CP10k.mean()-c.CP10k.mean() if len(c) and len(s) else np.nan
  loo.append(dict(filter=filt,scope=scope,gene=gene,omitted_donor=sample,n_control=len(c),n_case=len(s),full_difference=base,loo_difference=delta,same_sign=np.sign(delta)==np.sign(base)))
wr(pd.DataFrame(loo),'human_leave_one_donor_out.tsv')
# Human Core fixed-state composition: observed state proportions preserved for both groups.
parts=[]
for (filt,gene),a in h.groupby(['filter','gene']):
 sides={}
 for group in ['Exposed','Silicosis']:
  b=a[a.group==group];core=b[b.scope=='Core'].set_index('sample');sub=b[b.scope!='Core'].copy();sub['w']=sub.cells/sub['sample'].map(core.cells);sub['contribution']=sub.w*sub.CP10k;vv={}
  for state,z in sub.groupby('scope'):
   w=z.w.sum()/len(core);v=z.contribution.sum()/len(core)/w if w>0 else np.nan;vv[state]=(w,v)
  sides[group]=vv
 for state in ['NaiveCD4','MemCD4','Treg']:
  wc,mc=sides['Exposed'][state];ws,ms=sides['Silicosis'][state];parts.append(dict(filter=filt,gene=gene,state=state,control_fraction=wc,case_fraction=ws,control_conditional_expression=mc,case_conditional_expression=ms,composition_component=(ws-wc)*(mc+ms)/2,within_state_component=(ms-mc)*(wc+ws)/2))
wr(pd.DataFrame(parts),'human_Core_composition_decomposition.tsv')
# Reconcile direct raw expression with the published v2 sample table.
old=pd.read_csv(rec(V/'05_supporting_evidence/GSE183682_candidate_sample_expression.tsv'),sep='\t');new=met[(met.dataset=='GSE183682')&(met['filter']=='All')];j=old.merge(new,on=['scope','sample','gene'],suffixes=('_old','_new'));assert len(j)==len(old);err=float((j.CP10k-j.mean_CP10k).abs().max());assert err<1e-5,err
wr(pd.DataFrame(manifest).drop_duplicates('path'),'input_manifest.tsv')
# Direction grid distinguishes metrics and time points.
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axs=plt.subplots(1,2,figsize=(12,5))
for ax,gene in zip(axs,['MDFIC','TTC39C']):
 a=con[(con.gene==gene)&(con.scope=='ActCD4')&(con['filter']=='All')];metrics=['mean_CP10k','pseudobulk_CPM','mean_logCP10k','detection_fraction','positive_mean_CP10k'];matrix=np.array([[a[(a.dataset==ds)&(a.day==day)&(a.metric==m)].log2_case_control.iloc[0] for m in metrics] for ds,day in [('GSE314943',-1),('GSE183682',7),('GSE183682',56)]])
 im=ax.imshow(matrix,cmap='RdBu_r',vmin=-2,vmax=2,aspect='auto')
 for i in range(3):
  for k in range(5):ax.text(k,i,f'{matrix[i,k]:+.2f}',ha='center',va='center',color='white' if abs(matrix[i,k])>1.3 else 'black')
 ax.set(title=gene+' / ActCD4',xticks=range(5),xticklabels=['Cell-mean\nCP10k','Pseudobulk\nCPM','Cell-mean\nlogCP10k','Detection\nfraction','Positive-cell\nCP10k'],yticks=range(3),yticklabels=['Discovery','External D7','External D56'])
fig.colorbar(im,ax=axs.tolist(),label='log2(case/control) of the specified sample-mean estimand',shrink=.7);fig.suptitle('Direction audit: normalization and detection components',fontsize=15);fig.subplots_adjust(top=.87,bottom=.2,wspace=.35,right=.85);fig.savefig(O/'figures/direction_estimand_grid.png',dpi=170,bbox_inches='tight');plt.close(fig)
summary={'status':'PASS','raw_feature_mapping_unique':True,'external_v2_max_expression_difference':err,'mouse_sample_metric_rows':len(met),'direction_records':len(con),'marker_standardization_records':len(standard),'cell_minimum':'none','external_inferential_tests':0,'old_results_modified':False}
(O/'provenance/direction_complete.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary));print(con[(con.gene.isin(['MDFIC','TTC39C']))&(con.scope=='ActCD4')&(con['filter']=='All')][['dataset','day','gene','metric','control','case','difference']].to_string(index=False));print('MARKER_STD',pd.DataFrame(standard).query("gene in ['MDFIC','TTC39C']")[['dataset','day','gene','panel','standardized_difference','sign_changed','min_cell_coverage']].to_string(index=False))